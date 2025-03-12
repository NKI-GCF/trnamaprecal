/// Program that adjusts the BWA mapping quality to somehting more usable for aligning tRNA
/// sequences.
///
/// The tRNA reference sequences differ as little as 1bp between genes, but many reads do map with
/// a higher score to these sequences than any other sequences as reported in the XB tag
/// (use bwa -u). If this is reflected in the mapping score it is easier to filter the reads
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use anyhow::{anyhow, Error, Result};
use atoi::atoi;
use clap::Parser;
use indexmap::IndexSet;
use log::*;
use noodles_bam as bam;
use noodles_fasta::io::Reader as FastaReader;
use noodles_sam::alignment::io::Write as WriteSam;
use noodles_sam::alignment::record::{data::field::Tag, MappingQuality};
use noodles_sam::alignment::record_buf::{data::field::Value, RecordBuf};
use simple_logger::SimpleLogger;
use triple_accel::levenshtein;

#[derive(Debug, Default)]
struct Reference {
    names: IndexSet<Vec<u8>>,
    sequences: IndexSet<Vec<u8>>,
    distances: DistanceMatrix,
    codon_map: IndexSet<Vec<u8>>,
    sequence_codon: Vec<usize>,
}

#[derive(Debug, Default)]
struct DistanceMatrix(HashMap<(usize, usize), u32>);

impl Reference {
    pub fn from_fasta<P: AsRef<Path>>(p: P) -> Result<Self> {
        let mut reader = File::open(p.as_ref())
            .map(BufReader::new)
            .map(FastaReader::new)?;

        let mut reference = Self::default();
        for r in reader.records() {
            let record = r?;
            let name = record.definition().name().to_vec();
            let seq = record.sequence().as_ref().to_vec();

            //extract the target codon name
            let codon_end = name.windows(5).position(|w| w == b"-trna")
                .ok_or(anyhow!("refname codon not extractable"))?;

            let (idx, _) = reference.codon_map.insert_full(name[0..codon_end].to_vec());
            reference.sequence_codon.push(idx);
            reference.names.insert(name);
            reference.sequences.insert(seq);
        }

        reference.calc_dist_matrix();

        Ok(reference)
    }
    
    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn find_chr(&self, s: &[u8]) -> Option<usize> {
        self.names.get_index_of(s)
    }

    pub fn calc_dist_matrix(&mut self) {
        let l = self.sequences.len();
        for a in 0..(l - 1) {
            for b in (a + 1)..l {
                let dist = levenshtein(
                    self.sequences.get_index(a).unwrap(),
                    self.sequences.get_index(b).unwrap(),
                );
                self.distances.0.insert((a, b), dist);
            }
        }
    }

    pub fn get_dist(&self, a: usize, b: usize) -> Option<u32> {
        if a == b {
            Some(0)
        } else if a < b {
            self.distances.0.get(&(a, b)).copied()
        } else {
            self.distances.0.get(&(b, a)).copied()
        }
    }
}

type MapQCounts = [u64; 61];

#[derive(Debug, Default)]
struct ReferenceCounts {
    refcounts: Vec<MapQCounts>,
    multi_same_target: Vec<usize>,
    multi_same_codon: Vec<usize>,
    unmapped: usize,
}

impl ReferenceCounts {

    pub fn new(n: usize) -> ReferenceCounts {
        let refcounts = (0..n).map(|_| [0; 61]).collect();
        ReferenceCounts {
            refcounts,
            multi_same_target: vec![0; n],
            multi_same_codon: vec![0; n],
            unmapped: 0
        }
    }

    pub fn count_mapq(&mut self, tid: usize, mapq: u8) {
        self.refcounts[tid][mapq as usize] += 1;
    }

    pub fn count_same_target(&mut self, tid: usize) {
        self.multi_same_target[tid] += 1;
    }

    pub fn count_same_codon(&mut self, tid: usize) {
        self.multi_same_codon[tid] += 1;
    }

    pub fn count_unmapped(&mut self) {
        self.unmapped += 1;
    }

    pub fn write<P: AsRef<Path>>(&self, path: P, reference: &Reference) -> Result<()> {
        let file = File::create(path.as_ref())?;
        let mut writer = BufWriter::new(file);

        for (refidx, mapqcounts) in self.refcounts.iter().enumerate() {
            let refname = String::from_utf8_lossy(&reference.names[refidx]);
            writeln!(&mut writer, "{refname}\tQ0_sametarget\t{}", self.multi_same_target[refidx])?;
            writeln!(&mut writer, "{refname}\tQ0_samecodon\t{}", self.multi_same_codon[refidx])?;
            for (qval, qscore) in mapqcounts.iter().enumerate() {
                writeln!(&mut writer, "{}\t{}\t{}", refname, qval, qscore)?;
            }
            writeln!(&mut writer, "*\t*\t{}", self.unmapped)?;
        }

        Ok(())
    }
}

struct AltMapping {
    target_id: usize,
    target_codon_id: usize,
    alignment_score: u32,
    #[allow(dead_code)]
    mismatches: u32,
    #[allow(dead_code)]
    mapq: u8,
}

// not used yet
#[allow(dead_code)]
struct AltMappingComp {
    same_target: bool,
    same_codon: bool,
    map_len_diff: isize,
    as_diff: i32,
    ref_dist: i64,
}

impl AltMapping {
    //not used yet
    #[allow(dead_code)]
    pub fn from_record(r: &RecordBuf, reference: &Reference) -> AltMapping {
        let mapq = r.mapping_quality().unwrap().get();
        let target_id = r.reference_sequence_id().unwrap();
        let alignment_score = r
            .data()
            .get(&Tag::ALIGNMENT_SCORE)
            .unwrap()
            .as_int()
            .unwrap() as u32;
        let mismatches = r
            .data()
            .get(&Tag::MISMATCHED_POSITIONS)
            .unwrap()
            .as_int()
            .unwrap() as u32;

        AltMapping {
            target_id,
            target_codon_id: reference.sequence_codon[target_id],
            alignment_score,
            mismatches,
            mapq,
        }
    }

    pub fn from_xb(xb: &Xb, reference: &Reference) -> AltMapping {
        let target_id = reference.find_chr(xb.refname).expect("Alt chr not in ref");

        AltMapping {
            target_id,
            target_codon_id: reference.sequence_codon[target_id],
            alignment_score: xb.alignment_score,
            mismatches: xb.mismatches,
            mapq: xb.mapq,
        }
    }

    //not used yet
    #[allow(dead_code)]
    pub fn compare_primary(&self, primary_mapping: AltMapping, reference: &Reference) -> AltMappingComp {
        AltMappingComp {
            same_target: self.target_id == primary_mapping.target_id,
            same_codon: self.target_codon_id == primary_mapping.target_codon_id,
            map_len_diff: 0,
            as_diff: primary_mapping.alignment_score as i32 - self.alignment_score as i32,
            ref_dist: reference.get_dist(primary_mapping.target_id, self.target_id).unwrap() as i64
        }
    }
}

#[derive(Default)]
struct AltMappings(Vec<AltMapping>);
impl AltMappings {
    pub fn all_same_codon(&self, codon_id: usize) -> bool {
        self.0.is_empty() ||
            self.0.iter().all(|a| a.target_codon_id == codon_id)
    }

    pub fn all_same_target(&self, target_id: usize) -> bool {
        self.0.is_empty() ||
            self.0.iter().all(|a| a.target_id == target_id)
    }

    pub fn iter(&self) -> impl Iterator<Item=&AltMapping> {
        self.0.iter()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }
}

#[allow(dead_code)]
struct Xb<'a> {
    refname: &'a [u8],
    strand_fw: bool,
    pos: i32,
    cigar: &'a [u8],
    alignment_score: u32,
    mismatches: u32,
    mapq: u8,
}

impl<'a> TryFrom<&'a [u8]> for Xb<'a> {
    type Error = Error;

    fn try_from(s: &'a [u8]) -> Result<Self, Self::Error> {
        let mut it = s.split(|&b| b == b',');

        let refname = it.next().ok_or(anyhow!("No refname in XB"))?;
        let strandpos = it.next().ok_or(anyhow!("No strand/pos in XB"))?;
        let strand_fw = strandpos[0] == b'+';
        let pos = atoi(&strandpos[1..]).ok_or(anyhow!("Cannot parse XB mapq"))?;
        let cigar = it.next().ok_or(anyhow!("No cigar in XB"))?;
        let mismatches = it
            .next()
            .ok_or(anyhow!("No MM in XB"))
            .and_then(|v| atoi::<u32>(v).ok_or(anyhow!("Cannot parse XB MM")))?;
        let alignment_score = it
            .next()
            .ok_or(anyhow!("No AS in XB"))
            .and_then(|v| atoi::<u32>(v).ok_or(anyhow!("Cannot parse XB AS")))?;
        let mapq = it
            .next()
            .ok_or(anyhow!("No mapq in XB"))
            .and_then(|v| atoi::<u8>(v).ok_or(anyhow!("Cannot parse XB mapq")))?;

        Ok(Xb {
            refname,
            strand_fw,
            pos,
            cigar,
            mismatches,
            alignment_score,
            mapq,
        })
    }
}

fn recalc_mapq(r: &RecordBuf, reference: &Reference, alts: &AltMappings) -> u8 {
    let mq = r.mapping_quality().unwrap().get();
    let primary_tid = r.reference_sequence_id().unwrap();
    let primarty_as = r
        .data()
        .get(&Tag::ALIGNMENT_SCORE)
        .unwrap()
        .as_int()
        .unwrap();
    //eprintln!("primary as {primarty_as} ");
    let md = alts
        .iter()
        .map(|alt| {
            //get the difference in ref seqs
            if primary_tid == alt.target_id {
                return 0;
            }
            let _dist = reference.get_dist(primary_tid, alt.target_id).unwrap() as i64;
            //eprintln!("dist {dist} alt score {} ", alt.alignment_score);
            let score = primarty_as - alt.alignment_score as i64;
            if score < 0 {
                0
            } else {
                score as u8
            }
        })
        .min()
        .unwrap_or(0);

    //eprintln!("mapq was {mq} altmin {md}");
    std::cmp::max(mq, md)
}


#[derive(Parser)]
struct Config {
    /// The fasta reference file that was used for the alignment
    #[arg(long = "reference")]
    reference: PathBuf,
    /// The BAM file containing the BWA alignments to the tRNA reference
    #[arg(long = "bam")]
    bam: PathBuf,
    #[arg(long = "out")]
    /// The output BAM file
    out: PathBuf,
    /// The outcount quantification file (grouped by reference sequence and mapping quality)
    #[arg(long = "counts")]
    counts: Option<PathBuf>,
    /// The minimum read quality to include in the count files (this does not affect the BAM output
    #[arg(long = "qs", default_value_t = 0.0)]
    qs: f32,
}

fn main() -> Result<()> {
    SimpleLogger::new().init().unwrap();

    let conf = Config::parse();
    info!("Processing file {}", conf.bam.display());

    let reference = Reference::from_fasta(&conf.reference)?;

    let mut reader = bam::io::reader::Builder::default().build_from_path(&conf.bam)?;
    let header = reader.read_header()?;

    // check that the header chr order is indentical to the reference .fa
    header
        .reference_sequences()
        .keys()
        .enumerate()
        .all(|(i, k)| {
            k == reference
                .names
                .get_index(i)
                .expect("Bam index out of range for reference sequence")
        });

    //open the output bam
    let out = File::create(&conf.out)?;
    let mut writer = bam::io::Writer::new(out);
    writer.write_header(&header)?;

    let xb_tag: Tag = (*b"XB").into();
    let qs_tag: Tag = (*b"qs").into();

    let mut refcounts = ReferenceCounts::new(reference.len());
    let mut lowqs = 0;
    let mut warnqs = false;
    let mut nrecords = 0;

    for result in reader.record_bufs(&header) {
        let mut record = result?;
        nrecords += 1;

        // Collect the alternative mappings from the XB tag
        let alt_mappings = if let Some(xb) = record.data().get(&xb_tag) {
            if let Value::String(ts) = xb {
                let am = ts.split_inclusive(|&b| b == b';')
                    .map(|s| Xb::try_from(s).map(|xb| AltMapping::from_xb(&xb, &reference)))
                    .collect::<Result<Vec<_>, _>>()?;
                AltMappings(am)
            } else {
                warn!("XB tag value present but not string type. Ignoring");
                AltMappings::default()
            }
        } else {
            AltMappings::default()
        };

        // Recalculate the mapping quality by comparing the Alternative alignments to the primary
        // one

        if !alt_mappings.is_empty() {
            let new_mapq = recalc_mapq(&record, &reference, &alt_mappings);
            let mq = MappingQuality::new(new_mapq);
            if let Some(old) = record.mapping_quality_mut().replace(mq.unwrap()) {
                record.data_mut().insert(Tag::from(*b"om"), Value::UInt8(old.into()));
            }
        }

        writer.write_alignment_record(&header, &record)?;

        // Record the read mapping if the read has mapq > 0 or 
        // all alternative alignments fall on the same tRNA transcript
        if let Some(refid) = record.reference_sequence_id() {
            if let Some(Value::Float(qs)) = record.data().get(&qs_tag) {
                if *qs < conf.qs {
                    lowqs += 1;
                    continue;
                } 
            } else {
                if !warnqs {
                    warn!("qs tag not present in aligned data, unable to filter on mean read quality");
                    warnqs = true;
                }
            }

            if !record.flags().is_reverse_complemented() {
                let mq: u8 = record.mapping_quality().expect("Mapped read wihout a mapping quality").into();
                if mq > 0 {
                    refcounts.count_mapq(refid, mq);
                } else {
                    let all_same_codon = alt_mappings.all_same_codon(reference.sequence_codon[refid]);
                    let all_same_target = alt_mappings.all_same_target(refid);

                    if all_same_codon {
                        refcounts.count_same_codon(refid);
                    }

                    if all_same_target {
                        refcounts.count_same_target(refid);
                    }

                    if !(all_same_codon || all_same_target) {
                        refcounts.count_mapq(refid, 0);
                    }
                }
            } else {
                refcounts.count_unmapped();
            }
        }
    }

    info!("Processed {} records into outfile: {}", nrecords, conf.out.display());

    if conf.qs > 0.0 {
        info!("Skipped {} records with a mean quality score below {}, ({:.2}%)", lowqs, conf.qs,
            100.0 * lowqs as f64  / nrecords as f64)
    }

    if let Some(counts_file) = &conf.counts {
        refcounts.write(counts_file, &reference)?;
        info!("Counts written to {}", counts_file.display());
    }

    info!("All done.");

    Ok(())
}
