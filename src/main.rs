use bio::io::fastq::{Error as FastqError, FastqRead, Reader, Record};
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

fn calc_read_stats(in_seq_file: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("Calculating read stats for file: {}", in_seq_file);

    let input = open_maybe_gzipped(in_seq_file)?;

    let bufreader = BufReader::new(input);
    let mut reader = Reader::new(bufreader);

    let mut record = Record::new();
    let mut total_records: i32 = 0;
    let mut gc_count: usize;
    let mut gc_percent: u8;
    let mut len: usize;
    let mut seq: &[u8];
    loop {
        match reader.read(&mut record) {
            Ok(()) => {
                if record.seq().is_empty() {
                    // End of valid content (e.g., trailing blank lines)
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
            Err(FastqError::ReadError(ref e)) if e.kind() == ErrorKind::UnexpectedEof => break,
            Err(FastqError::IncompleteRecord) => {
                return Err("Encountered incomplete FASTQ record".into());
            }
            Err(e) => return Err(Box::new(e)),
        }
    }

    println!("Total records: {}", total_records);
    Ok(())
}

fn main() {
    let in_seq_file = "-";
    //let in_seq_file = "seq.fastq.gz";
    //let in_seq_file = "seq.fastq";
    let _ = calc_read_stats(in_seq_file);
}
