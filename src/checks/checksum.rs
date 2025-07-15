use crate::checker::FileReport;
use crate::progress::DualProgressReader;
use crate::sha256::SharedHashingReader;
use indicatif::ProgressBar;
use sha2::{Digest, Sha256};
use std::io::BufReader;
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::{fs, io};

pub fn check_checksum_only(
    path: &Path,
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
    let mut progress_reader =
        DualProgressReader::new(hashing_reader, file_pb.clone(), global_pb.clone());

    if let Err(e) = io::copy(&mut progress_reader, &mut io::sink()) {
        return FileReport {
            path: path.to_path_buf(),
            stats: None,
            checksum: None,
            errors: vec![format!("Failed to read file: {}", e)],
            warnings: vec![],
        };
    }

    drop(progress_reader);

    let checksum = match Arc::try_unwrap(hasher) {
        Ok(mutex) => {
            let final_hasher = mutex.into_inner().unwrap();
            Some(format!("{:x}", final_hasher.finalize()))
        }
        Err(_) => {
            return FileReport {
                path: path.to_path_buf(),
                stats: None,
                checksum: None,
                errors: vec!["Failed to finalize checksum: hasher is still in use.".to_string()],
                warnings: vec![],
            };
        }
    };

    FileReport {
        path: path.to_path_buf(),
        stats: None,
        checksum,
        errors: vec![],
        warnings: vec![],
    }
}

#[derive(Debug)]
pub struct ChecksumOnlyJob {
    pub path: PathBuf,
    pub size: u64,
}
