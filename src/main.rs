use anyhow::{Context, Result};
use clap::Parser;
use fastq::PairReport;
use std::error::Error as StdError;
use std::fmt;
use std::path::PathBuf;

mod fastq;
mod progress;
mod sha256;

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
