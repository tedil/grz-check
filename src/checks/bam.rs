use crate::checker::{FileReport, Stats};
use crate::progress::DualProgressReader;
use crate::sha256::SharedHashingReader;
use indicatif::ProgressBar;
use noodles::bam;
use noodles::sam::alignment::record::cigar::op::Kind;
use sha2::{Digest, Sha256};
use std::fs;
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

pub fn check_single_bam(path: &Path, file_pb: &ProgressBar, global_pb: &ProgressBar) -> FileReport {
    file_pb.set_message(
        path.file_name()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string(),
    );

    let file = match fs::File::open(path) {
        Ok(f) => f,
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
    let mut reader = bam::io::Reader::new(progress_reader);

    let header = match reader.read_header() {
        Ok(h) => h,
        Err(e) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                checksum: None,
                errors: vec![format!("Failed to read BAM header: {}", e)],
                warnings: vec![],
            };
        }
    };

    let mut warnings = Vec::new();
    if !header.reference_sequences().is_empty()
        || !header.read_groups().is_empty()
        || (header.programs().roots().count() != 0)
        || !header.comments().is_empty()
    {
        warnings.push(
            "Detected a header in BAM file, ensure it contains no private information!".to_string(),
        );
    }

    let mut num_records = 0;
    for (i, result) in reader.records().enumerate() {
        let record = match result {
            Ok(rec) => rec,
            Err(e) => {
                return FileReport {
                    path: path.to_path_buf(),
                    stats: Some(Stats {
                        num_records,
                        read_length: None,
                    }),
                    checksum: None,
                    errors: vec![format!("Failed to parse record #{}: {}", i + 1, e)],
                    warnings,
                };
            }
        };
        num_records += 1;

        if record.flags().is_secondary() {
            warnings.push(format!(
                "Record #{} ('{}') is a secondary alignment.",
                num_records,
                record.name().map(|n| n.to_string()).unwrap_or_default()
            ));
        }

        if !record.flags().is_secondary()
            && record
                .cigar()
                .iter()
                .any(|op| op.is_ok_and(|op| op.kind() == Kind::HardClip))
        {
            warnings.push(format!(
                "Record #{} ('{}') is a primary alignment and contains hard-clipped bases.",
                num_records,
                record.name().map(|n| n.to_string()).unwrap_or_default()
            ));
        }
    }

    if num_records == 0 {
        return FileReport {
            path: path.to_path_buf(),
            stats: None,
            checksum: None,
            errors: vec!["File is empty. Expected at least one record.".to_string()],
            warnings: vec![],
        };
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
                    num_records,
                    read_length: None,
                }),
                checksum: None,
                errors: vec!["Failed to finalize checksum: hasher is still in use.".to_string()],
                warnings,
            };
        }
    };

    FileReport {
        path: path.to_path_buf(),
        stats: Some(Stats {
            num_records,
            read_length: None,
        }),
        checksum,
        errors: vec![],
        warnings,
    }
}

#[derive(Debug)]
pub struct BamCheckJob {
    pub path: PathBuf,
    pub size: u64,
}
