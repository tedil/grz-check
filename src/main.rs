use anyhow::{Context, Result};
use clap::Parser;
use fastq::PairReport;
use std::error::Error as StdError;
use std::fmt;
use std::path::PathBuf;

mod fastq;

#[derive(Debug, clap::Parser)]
#[command(author, version, about)]
struct Args {
    #[command(subcommand)]
    cmd: Commands,

    /// Flag to show progress bars during processing.
    #[arg(long)]
    show_progress: Option<bool>,
}

#[derive(Debug, Parser)]
enum Commands {
    /// Checks integrity and consistency of FASTQ files.
    ///
    /// Use --paired or --single for each sample. These flags can be used multiple times.
    ///
    /// By default, the tool will exit immediately after the first error is found.
    /// Use --continue-on-error to check all files regardless of errors.
    Fastq {
        /// A paired-end sample. Provide FQ1, FQ2, FQ1 read length, and FQ2 read length.
        ///
        /// Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check.
        /// Example: --paired fq1.fastq.gz fq2.fastq.gz 150 151
        #[arg(
            long,
            action = clap::ArgAction::Append,
            num_args = 4,
            value_names = ["FQ1_PATH", "FQ2_PATH", "FQ1_READ_LEN", "FQ2_READ_LEN"]
        )]
        paired: Vec<String>,

        /// A single-end sample. Provide the file path and read length.
        ///
        /// Read Length: >0 for fixed, 0 for auto-detect, <0 to skip length check.
        /// Example: --single sample.fastq.gz 0
        #[arg(
            long,
            action = clap::ArgAction::Append,
            num_args = 2,
            value_names = ["FQ_PATH", "READ_LEN"]
        )]
        single: Vec<String>,

        /// Path to write the output TSV report.
        #[arg(long, required = true)]
        output: PathBuf,

        /// Continue processing all files even if an error is found.
        #[arg(long, action = clap::ArgAction::SetTrue)]
        continue_on_error: bool,

        /// Number of threads to use for processing.
        #[arg(long)]
        threads: Option<usize>,
    },
}

#[derive(Debug)]
struct EarlyExitError(PairReport);

impl fmt::Display for EarlyExitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "A validation error occurred, exiting.")
    }
}

impl StdError for EarlyExitError {}

fn main() -> Result<()> {
    let args = Args::parse();

    match args.cmd {
        Commands::Fastq {
            paired,
            single,
            output,
            threads,
            continue_on_error,
        } => {
            if let Some(num_threads) = threads {
                rayon::ThreadPoolBuilder::new()
                    .num_threads(num_threads)
                    .build_global()
                    .context("Failed to set up Rayon thread pool")?;
            }
            fastq::run_check(
                paired,
                single,
                &output,
                continue_on_error,
                args.show_progress,
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fastq::run_check;
    use anyhow::Result;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::fs;
    use std::io::Write;
    use std::path::Path;
    use tempfile::tempdir;

    fn create_gzipped_fastq(path: &Path, content: &str) -> Result<()> {
        let file = fs::File::create(path)?;
        let mut writer = GzEncoder::new(file, Compression::default());
        writer.write_all(content.as_bytes())?;
        writer.finish()?;
        Ok(())
    }

    struct TestFiles {
        _tempdir: tempfile::TempDir,
        pub dir: PathBuf,
    }

    impl TestFiles {
        fn new() -> Result<Self> {
            let tempdir = tempdir()?;
            let dir = tempdir.path().to_path_buf();

            // Case 1: Valid pair
            create_gzipped_fastq(
                &dir.join("ok_r1.fastq.gz"),
                "@SEQ1\nACGT\n+\nFFFF\n@SEQ2\nTGCA\n+\nFFFF\n",
            )?;
            create_gzipped_fastq(
                &dir.join("ok_r2.fastq.gz"),
                "@SEQ1\nAAAA\n+\nFFFF\n@SEQ2\nTTTT\n+\nFFFF\n",
            )?;

            // Case 2: Inconsistent length in a single file
            create_gzipped_fastq(
                &dir.join("badlen.fastq.gz"),
                "@SEQ1\nACGT\n+\nFFFF\n@SEQ2\nTCG\n+\nFFF\n",
            )?;

            // Case 3: Mismatched read counts in a pair
            create_gzipped_fastq(
                &dir.join("counts1.fastq.gz"),
                "@SEQ1\nACGT\n+\nFFFF\n@SEQ2\nTGCA\n+\nFFFF\n",
            )?;
            create_gzipped_fastq(&dir.join("counts2.fastq.gz"), "@SEQ1\nACGT\n+\nFFFF\n")?;

            // Case 4: Malformed file (seq len != qual len)
            create_gzipped_fastq(&dir.join("malformed.fastq.gz"), "@SEQ1\nACGT\n+\nFF\n")?;

            // Case 5: Empty file
            create_gzipped_fastq(&dir.join("empty.fastq.gz"), "")?;

            // Case 6: Read length 5 for R2 (for checking different read lengths in the same pair)
            create_gzipped_fastq(
                &dir.join("ok_r2_len5.fastq.gz"),
                "@SEQ1\nAAAAA\n+\nFFFFF\n@SEQ2\nTTTTA\n+\nFFFFF\n",
            )?;

            Ok(Self {
                _tempdir: tempdir,
                dir,
            })
        }

        fn path_str(&self, filename: &str) -> String {
            self.dir.join(filename).to_string_lossy().to_string()
        }
    }

    fn read_report_status(report_path: &Path) -> Result<Vec<(String, String, String)>> {
        let mut reader = csv::ReaderBuilder::default()
            .delimiter(b'\t')
            .from_path(&report_path)?;
        reader
            .records()
            .map(|r| {
                let record = r?;
                Ok((
                    record[0].to_string(),
                    record[1].to_string(),
                    record[4].to_string(),
                ))
            })
            .collect()
    }

    #[test]
    fn test_valid_pair_with_different_lengths() -> Result<()> {
        let fixture = TestFiles::new()?;
        let output = fixture.dir.join("report.tsv");
        let paired = vec![
            fixture.path_str("ok_r1.fastq.gz"),
            fixture.path_str("ok_r2_len5.fastq.gz"),
            "4".to_string(),
            "5".to_string(),
        ];

        run_check(paired, vec![], &output, false, Some(false))?;

        let statuses = read_report_status(&output)?;
        assert_eq!(statuses.len(), 2);
        assert_eq!(statuses[0].1, "OK");
        assert_eq!(statuses[1].1, "OK");
        Ok(())
    }

    #[test]
    fn test_multiple_inputs_with_continue_on_error() -> Result<()> {
        let fixture = TestFiles::new()?;
        let output = fixture.dir.join("report.tsv");

        let all_paired = vec![
            fixture.path_str("counts1.fastq.gz"),
            fixture.path_str("counts2.fastq.gz"),
            "0".to_string(),
            "0".to_string(),
            fixture.path_str("ok_r1.fastq.gz"),
            fixture.path_str("ok_r2.fastq.gz"),
            "4".to_string(),
            "4".to_string(),
        ];

        let all_single = vec![fixture.path_str("badlen.fastq.gz"), "0".to_string()];

        run_check(all_paired, all_single, &output, true, Some(false))?;

        let statuses = read_report_status(&output)?;
        assert_eq!(statuses.len(), 5);

        assert_eq!(statuses[0].0, fixture.path_str("counts1.fastq.gz"));
        assert_eq!(statuses[0].1, "ERROR");
        assert!(statuses[0].2.contains("Mismatched read counts"));

        assert_eq!(statuses[2].0, fixture.path_str("ok_r1.fastq.gz"));
        assert_eq!(statuses[2].1, "OK");

        assert_eq!(statuses[4].0, fixture.path_str("badlen.fastq.gz"));
        assert_eq!(statuses[4].1, "ERROR");
        Ok(())
    }
}
