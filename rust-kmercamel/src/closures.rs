// closures.rs
// Functions for FASTA parsing and invoking all necessary functions for mapping and alignment.

use std::io::{self};
use std::error::Error;
use std::io::{BufRead, BufReader, Write, BufWriter};
use std::path::Path;
use crate::BufReadDecompressor;
use std::fs::{File};
use std::sync::{Arc};
use seq_io::BaseRecord;
use seq_io::parallel::{read_process_fasta_records, read_process_fastq_records};
use dashmap::{DashMap, DashSet};
use std::path::PathBuf;
use super::Params;
use crate::get_reader;
use crate::parse::*;
use crate::kmers::*;
use std::time::Instant;
//use crate::align::{get_slices, align_slices, AlignStats};
use std::sync::atomic::{AtomicUsize, Ordering};
use rust_parallelfastx::parallel_fastx;
use std::str::FromStr;
use std::sync::mpsc;

// Main function for all FASTA parsing + mapping / alignment functions.
pub fn run_mers(filename: &PathBuf, params: &Params, threads: usize, queue_len: usize, fasta_reads: bool, output_prefix: &PathBuf) {
    let kmers = DashSet::new();
    let ref_i = AtomicUsize::new(0);
    let ref_map : DashMap<usize, (String, usize)> = DashMap::new(); // Sequence lengths per reference

    // PAF file generation
    let paf_filename = format!("{}{}", output_prefix.to_str().unwrap(), ".paf");
    let mut paf_file = match File::create(&paf_filename) {
        Err(why) => panic!("Couldn't create {}: {}", paf_filename, why.description()),
        Ok(paf_file) => BufWriter::new(paf_file),
    };
    // Closure for indexing reference k-min-mers
    let index_mers = |seq_id: &str, seq: &[u8], params: &Params| -> usize {
        let origin = ref_i.fetch_add(1, Ordering::Relaxed);
        add_kmers(&kmers, seq, params);
        //let nb_mers = mers::ref_extract(ref_idx, seq, params, &mers_index);
        ref_map.insert(origin, (seq_id.to_string(), seq.len()));
        1
    };

    // Closures for obtaining k-min-mers from references

    let ref_process_read_aux_mer = |ref_str: &[u8], ref_id: &str| -> Option<u64> {
        let nb_mers = index_mers(ref_id, ref_str, params);
        return Some(1)
    };

    let ref_process_read_fasta_mer = |record: seq_io::fasta::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap();
        *found = ref_process_read_aux_mer(&ref_str, &ref_id);

    };
    let ref_process_read_fastq_mer = |record: seq_io::fastq::RefRecord, found: &mut Option<u64>| {
        let ref_str = record.seq(); 
        let ref_id = record.id().unwrap();
        *found = ref_process_read_aux_mer(&ref_str, &ref_id);
    };
    let ref_main_thread_mer = |_found: &mut Option<u64>| { // runs in main thread
        None::<()>
    };

    //
    
    // Start processing references

    let start = Instant::now();
    let (buf,_dontcare) = get_reader(&filename);
    if fasta_reads {
        let reader = seq_io::fasta::Reader::new(buf);
        read_process_fasta_records(reader, threads as u32, queue_len, ref_process_read_fasta_mer, |_record, found| {ref_main_thread_mer(found)}).ok();
    }
    else {
        let reader = seq_io::fastq::Reader::new(buf);
        read_process_fastq_records(reader, threads as u32, queue_len, ref_process_read_fastq_mer, |_record, found| {ref_main_thread_mer(found)}).ok();
    }
    eprintln!("Added {} unique k-mers.", kmers.len());
}
