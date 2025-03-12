use std::collections::HashMap;
use std::io::{BufRead, BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use crossbeam_channel::{Sender, Receiver, unbounded, bounded};
use std::thread;


use anyhow::Result;
use clap::Parser;
use log::*;
use noodles_bam as bam;
use noodles_sam as sam;
use noodles_bam::record::Record as BamRecord;
use noodles_sam::Header;
use noodles_sam::alignment::io::Write as SamWrite;
use noodles_sam::alignment::RecordBuf;
use simple_logger::SimpleLogger;

#[derive(Debug, Parser)]
struct Config {
    #[arg(long = "reference")]
    reference: PathBuf,
    #[arg(long = "bam")]
    bam: PathBuf,
    #[arg(long = "out")]
    out: PathBuf,
    #[arg(long = "bwa", default_value = "bwa")]
    bwa: PathBuf,
    #[arg(long = "args", default_value="-t8 -Y -h30,200 -W13 -k6 -xont2d -T20 -u")]
    args: String,
}


enum SamResult {
    Header(Header),
    Record(RecordBuf),
}

fn open_bwa_process(bwa_path: &Path, reference: &Path, args: &str) -> Result<(Sender<String>, Receiver<SamResult>)> {
    let (tx_fq, rx_fq) = bounded::<String>(1000);

    let bwa_args: Vec<_> = args.split_whitespace().collect();
    debug!("bwa args: {:?}", &bwa_args);
    let mut child = Command::new(bwa_path)
        .arg("mem")
        .args(bwa_args)
        .arg(reference)
        .arg("-")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap_or_else(|e| panic!("Error spawning aligner process:\n\t{}", e));

    let mut childin = child.stdin.take().unwrap();
    let childout = child.stdout.take().unwrap();
    let childerr = child.stderr.take().unwrap();

    // the stderr logger
    thread::spawn(move || {
        let er = BufReader::new(childerr);
        for line in er.lines() {
            let bl = line.expect("Error reading BWA stderr stream");
            info!("BWA: '{bl}'");
        }
    });

    thread::spawn(move || {
        for s in rx_fq.iter() {
            childin
                .write_all(s.as_bytes())
                .expect("Error write record to aligner");
            }
        info!("align_thread: Sent all records to aligner");
        drop(childin);
        if child
            .wait()
                .expect("Error in aligner command result")
                .success()
        {
            info!("BWA finished mapping");
        } else {
            panic!("Aligner return with non-zero exit code");
        }
    });


    let (tx_sam, rx_sam) = unbounded();
    // from another thread read the sam out put from stdout and parse it
    thread::spawn(move || {
        // read sam file from child process
        let mut reader = sam::io::reader::Builder::default()
            .build_from_reader(childout)
            .expect("Error opening BWA output as SAM");
        let header = reader.read_header()
            .expect("Error reading BWA sam output header");
        tx_sam.send(SamResult::Header(header.clone()))
            .expect("Error sending SAM header to main");

        for (i, record) in reader.record_bufs(&header).enumerate() {
            tx_sam.send(SamResult::Record(record.expect("Error parsing SAM record from BWA output")))
                .unwrap_or_else(|_| panic!("Error sending BWA record to main {i}"));
        }
        drop(tx_sam);
        info!("all reads mapped");
    });


    Ok((tx_fq, rx_sam))
}

fn record_fastq(r: &BamRecord) -> String {
    let seq = r.sequence();
    let quals = r.quality_scores();
    let name = r.name().unwrap();

    let mut fastq_record = String::with_capacity(seq.len() * 2 + name.len() + 2 + 4);
    fastq_record.push('@');
    fastq_record.push_str(&String::from_utf8_lossy(name.as_ref()));
    fastq_record.push('\n');
    fastq_record.extend(seq.iter().map(|c| c as char));
    fastq_record.push('\n');

    fastq_record.push('+');
    fastq_record.push('\n');
    fastq_record.extend(quals.as_ref().iter().map(|&v| (v + 33) as char));
    fastq_record.push('\n');

    fastq_record
}

fn merge_data(unmapped: &BamRecord, mapped: &mut RecordBuf) -> Result<()> {
    //merge the aux tags
    let samdata = mapped.data_mut();

    for d in unmapped.data().iter() {
        let (tag, value) = d?;
        if samdata.get_index_of(&tag).is_none() {
            samdata.insert(tag, value.try_into()?);
        }
    }

    Ok(())
}

fn main() -> Result<()> {
    SimpleLogger::new().init().unwrap();

    let config = Config::parse();

    // check for bwa at the provided arg or in path
    // Start the aligner thread
    let (to_bwa, bwa_output) = open_bwa_process(&config.bwa, &config.reference, &config.args)?;

    // read bam file or stdin
    let mut reader = bam::io::reader::Builder::default().build_from_path(&config.bam)?;
    let _header = reader.read_header()?;

    let mut writer = bam::io::writer::Builder::default().
        build_from_path(&config.out)?;

    let mut record_map = HashMap::new();

    let mut header = None;

    // read from bam in a thread
    for (i, result) in reader.records().enumerate() {
        let record = result?;

        let fq = record_fastq(&record);
        let name = record.name().unwrap();
        // remove \0 string terminaator
        let rid: Vec<u8> = name.to_vec();
        record_map.insert(rid, record);

        to_bwa.send(fq)?;

        // try to read from bwa output
        if i % 1000 == 0 {
            match bwa_output.try_recv() {
                Ok(SamResult::Header(h)) => {
                    debug!("Received header from bwa");
                    writer.write_header(&h)?;
                    header = Some(h);
                },
                Ok(SamResult::Record(mut mapped)) => {
                    if mapped.flags().is_supplementary() {
                        continue;
                    }
                    let unmapped = record_map.remove(mapped.name().unwrap())
                        .expect("Got unknown record from BWA");
                    merge_data(&unmapped, &mut mapped)?;
                    writer.write_alignment_record(header.as_ref().unwrap(), &mapped)?;
                },
                Err(_) => {},
            }
        }
    }

    drop(to_bwa);

    //receive the remainder of the mapped records
    loop {
        match bwa_output.recv() {
            Ok(SamResult::Header(h)) => {
                debug!("Received header from bwa");
                writer.write_header(&h)?;
                header = Some(h);
            },
            Ok(SamResult::Record(mut mapped)) => {
                if mapped.flags().is_supplementary() {
                    continue;
                }
                let unmapped = record_map.remove(mapped.name().unwrap())
                    .expect("Got unknown record from BWA");
                merge_data(&unmapped, &mut mapped)?;
                writer.write_alignment_record(header.as_ref().unwrap(), &mapped)?;
            },
            Err(_) => break,
        }
    }
    info!("All done");

    Ok(())
}
