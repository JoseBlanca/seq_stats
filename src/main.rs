use bio::io::fastq::{Error as FastqError, FastqRead, Reader, Record};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, ErrorKind, Read};

/// Opens a file or stdin, and decompresses if gzipped.
/// `"-"` means read from stdin.
fn open_maybe_gzipped(path: &str) -> Result<Box<dyn Read>, Box<dyn std::error::Error>> {
    let raw_input: Box<dyn Read> = if path == "-" {
        Box::new(io::stdin())
    } else {
        Box::new(File::open(path)?)
    };

    let mut buf_reader = BufReader::new(raw_input);

    // Use fill_buf to peek at the buffer without consuming bytes
    let buffer = buf_reader.fill_buf()?;
    let is_gzipped = buffer.starts_with(&[0x1f, 0x8b]);

    // Now wrap in gzip decoder or not, preserving buffer
    if is_gzipped {
        Ok(Box::new(MultiGzDecoder::new(buf_reader)))
    } else {
        Ok(Box::new(buf_reader))
    }
}

fn calc_read_stats(in_seq_fpath: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("Calculating read stats for file: {}", in_seq_fpath);

    let input = open_maybe_gzipped(in_seq_fpath)?;

    let bufreader = BufReader::new(input);
    let mut reader = Reader::new(bufreader);

    let mut record = Record::new();
    let mut total_records: i32 = 0;
    let mut gc_count: usize;
    let mut gc_percent: u8;
    let mut len: usize;
    let mut seq: &[u8];

    loop {
        // Try to read the next record
        let read_result = reader.read(&mut record);

        // Handle possible errors
        if let Err(e) = read_result {
            match e {
                FastqError::ReadError(ref io_err) if io_err.kind() == ErrorKind::UnexpectedEof => {
                    break
                }
                FastqError::IncompleteRecord => {
                    return Err("Encountered incomplete FASTQ record".into())
                }
                other => return Err(Box::new(other)),
            }
        }

        // Check for implicit EOF case (blank trailing lines)
        if record.seq().is_empty() {
            break;
        }

        seq = record.seq();
        gc_count = seq
            .iter()
            .filter(|&&b| b == b'G' || b == b'g' || b == b'C' || b == b'c')
            .count();
        len = seq.len();
        gc_percent = ((gc_count as f64 / len as f64) * 100.0).round() as u8;
        println!("{len} {gc_percent}");
        total_records += 1;
    }

    println!("Total records: {}", total_records);
    Ok(())
}

#[derive(Parser)]
#[command(
    name = "seq_stats",
    version = "0.1",
    about = "Calculate GC content and length of sequences in a FASTQ file"
)]
struct Cli {
    /// Input file path (default: "-" (stdin))
    #[arg(default_value = "-")]
    input_seq: String,

    /// output seq file path (default: "-" (stdout))
    #[arg(default_value = "-")]
    output_seq: String,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Cli::parse();
    let in_seq_fpath = &args.input_seq;

    let _ = calc_read_stats(in_seq_fpath);

    Ok(())
}
