use crate::checker::{FileReport, Stats};
use crate::progress::DualProgressReader;
use crate::sha256::SharedHashingReader;
use indicatif::ProgressBar;
use sha2::{Digest, Sha256};
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

#[derive(Debug, Copy, Clone)]
pub enum ReadLengthCheck {
    Fixed(usize),
    Auto,
    Skip,
}

pub fn check_single_fastq(
    path: &Path,
    length_check: ReadLengthCheck,
    file_pb: &ProgressBar,
    global_pb: &ProgressBar,
) -> FileReport {
    file_pb.set_message(
        path.file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string(),
    );

    let file = match fs::File::open(path) {
        Ok(file) => file,
        Err(e) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                checksum: None,
                errors: vec![format!("Failed to open file for reading: {}", e)],
                warnings: vec![],
            };
        }
    };
    let hasher = Arc::new(Mutex::new(Sha256::new()));
    let hashing_reader = SharedHashingReader::new(BufReader::new(file), hasher.clone());
    let progress_reader =
        DualProgressReader::new(hashing_reader, file_pb.clone(), global_pb.clone());
    let decompressed_reader = match niffler::get_reader(Box::new(progress_reader)) {
        Ok((reader, _)) => reader,
        Err(e) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                checksum: None,
                errors: vec![format!("Failed to decompress file: {}", e)],
                warnings: vec![],
            };
        }
    };

    let mut reader = noodles::fastq::io::Reader::new(BufReader::new(decompressed_reader));
    let mut records = reader.records();

    let first_record = match records.next() {
        Some(Ok(rec)) => rec,
        Some(Err(e)) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                checksum: None,
                errors: vec![format!("Failed to parse first record: {}", e)],
                warnings: vec![],
            };
        }
        None => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                checksum: None,
                errors: vec!["File is empty. Expected at least one record.".to_string()],
                warnings: vec![],
            };
        }
    };

    let file_read_length = first_record.sequence().len();
    let mut num_reads: u64 = 1;

    if let ReadLengthCheck::Fixed(expected) = length_check {
        if expected != file_read_length {
            return FileReport {
                path: path.to_path_buf(),
                stats: Some(Stats {
                    num_records: num_reads,
                    read_length: Some(file_read_length),
                }),
                checksum: None,
                errors: vec![format!(
                    "Provided read length ({}) does not match first read's length ({}).",
                    expected, file_read_length
                )],
                warnings: vec![],
            };
        }
    }

    if !matches!(length_check, ReadLengthCheck::Skip) {
        for (i, record) in records.enumerate() {
            num_reads += 1;
            let record = match record {
                Ok(rec) => rec,
                Err(e) => {
                    return FileReport {
                        path: path.to_path_buf(),
                        stats: Some(Stats {
                            num_records: num_reads,
                            read_length: Some(file_read_length),
                        }),
                        checksum: None,
                        errors: vec![format!("Failed to parse record #{}: {}", i + 2, e)],
                        warnings: vec![],
                    };
                }
            };

            if record.sequence().len() != file_read_length {
                return FileReport {
                    path: path.to_path_buf(),
                    stats: Some(Stats {
                        num_records: num_reads,
                        read_length: Some(file_read_length),
                    }),
                    checksum: None,
                    errors: vec![format!(
                        "Found inconsistent read length at record #{}. Expected {}, but got {}.",
                        i + 2,
                        file_read_length,
                        record.sequence().len()
                    )],
                    warnings: vec![],
                };
            }
        }
    } else {
        num_reads += records.count() as u64;
    }

    drop(reader);

    let checksum = match Arc::try_unwrap(hasher) {
        Ok(mutex) => {
            let final_hasher = mutex.into_inner().unwrap();
            Some(format!("{:x}", final_hasher.finalize()))
        }
        Err(_) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: Some(Stats {
                    num_records: num_reads,
                    read_length: Some(file_read_length),
                }),
                checksum: None,
                errors: vec!["Failed to finalize checksum: hasher is still in use.".to_string()],
                warnings: vec![],
            };
        }
    };

    FileReport {
        path: path.to_path_buf(),
        stats: Some(Stats {
            num_records: num_reads,
            read_length: Some(file_read_length),
        }),
        checksum,
        errors: vec![],
        warnings: vec![],
    }
}

#[derive(Debug)]
pub struct FastqCheckJob {
    pub fq1: PathBuf,
    pub fq2: Option<PathBuf>,
    pub fq1_length_check: ReadLengthCheck,
    pub fq2_length_check: Option<ReadLengthCheck>,
    pub fq1_size: u64,
    pub fq2_size: Option<u64>,
}
