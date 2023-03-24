#![allow(unused_variables)]
#![allow(non_upper_case_globals)]
#![allow(warnings)]
#![feature(iter_advance_by)]

// kmercamel v0.1.0
// Copyright 2023 Baris Ekim.
// Licensed under the MIT license (http://opensource.org/licenses/MIT).
// This file may not be copied, modified, or distributed except according to those terms.
use std::io::stderr;
use std::error::Error;
use std::io::Write;
use std::collections::HashMap;
use std::collections::HashSet;
#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: tikv_jemallocator::Jemalloc = tikv_jemallocator::Jemalloc;
use std::fs::{self, File};
use std::io::{BufRead, BufReader};
use std::mem::{MaybeUninit};
use std::path::PathBuf;
use std::time::{Instant};
use chrono::{Utc};
use flate2::read::GzDecoder;
use lzzzz::lz4f::{BufReadDecompressor};
use rust_seq2kminmers::{FH, KH};
use structopt::StructOpt;
use std::cell::UnsafeCell;
use xx_bloomfilter::Bloom;
mod closures;
mod kmers;
mod parse;

pub type H = u64;

pub struct RacyBloom(UnsafeCell<Bloom>); // intentionnally allowing data races as a tradeoff for bloom speed
unsafe impl Sync for RacyBloom {}
// follows https://sodocumentation.net/rust/topic/6018/unsafe-guidelines
// don't try this at home
impl RacyBloom {
    fn new(v: Bloom) -> RacyBloom {
        RacyBloom(UnsafeCell::new(v))
    }
    fn get(&self) -> &mut Bloom {
        // UnsafeCell::get() returns a raw pointer to the value it contains
        // Dereferencing a raw pointer is also "unsafe"
        unsafe {&mut *self.0.get()}
    }
}

pub struct Params {
    k: usize,
    mask: i64,
    use_hpc: bool,
    use_simd: bool,
    debug: bool,
    complement: bool,
}

/// Try to get memory usage (resident set size) in bytes using the `getrusage()` function from libc.
// from https://github.com/digama0/mm0/blob/bebd670c5a77a1400913ebddec2c6248e76f90fe/mm0-rs/src/util.rs
fn get_memory_rusage() -> usize {
    let usage = unsafe {
        let mut usage = MaybeUninit::uninit();
        assert_eq!(libc::getrusage(libc::RUSAGE_SELF, usage.as_mut_ptr()), 0);
        usage.assume_init()
    };
    usage.ru_maxrss as usize * 1024
}

fn get_reader(path: &PathBuf) -> (Box<dyn BufRead + Send>, bool) {
    let mut filetype = "unzip";
    let filename_str = path.to_str().unwrap();
    let file = match File::open(path) {
            Ok(file) => file,
            Err(error) => panic!("Error opening compressed file: {:?}.", error),
        };
    if filename_str.ends_with(".gz")  {filetype = "zip";}
    if filename_str.ends_with(".lz4") {filetype = "lz4";}
    let reader :(Box<dyn BufRead + Send>,bool) = match filetype { 
        "zip" => (Box::new(BufReader::new(GzDecoder::new(file))),true), 
        "lz4" => (Box::new(BufReadDecompressor::new(BufReader::new(file)).unwrap()),true),
        _ =>     (Box::new(BufReader::new(file)),false), 
    }; 
    reader
}

#[derive(Debug, StructOpt)]
#[structopt(name = "kmercamel")]
/// Rust implementation of kmercamel.
struct Opt {
    /// Activate debug mode
    ///
    #[structopt(long)]
    debug: bool,
    /// Input file (raw or gzip-/lz4-compressed FASTX)
    ///
    /// Input file can be FASTA/FASTQ, as well as gzip-compressed (.gz) or
    /// lz4-compressed (.lz4). Lowercase bases are currently not supported;
    /// see documentation for formatting.
    #[structopt(parse(from_os_str))]
    reads: Option<PathBuf>,
    /// Output prefix for PAF file
    /// 
    #[structopt(parse(from_os_str), short, long)]
    prefix: Option<PathBuf>,
    /// k-min-mer length
    ///
    /// The length of each k-min-mer. If
    /// fewer l-mers than this value are obtained
    /// from a read, they will be ignored.
    #[structopt(short, long)]
    k: Option<usize>,
    /// Number of threads
    /// 
    #[structopt(long)]
    threads: Option<usize>,
    /// Enable low-memory reference FASTA parsing
    /// 
    #[structopt(long)]
    low_memory: bool,
    /// Deactivate SIMD (AVX2,AVX512) functions (for old processors)
    #[structopt(long)]
    nosimd: bool,
    /// Deactivate HomoPolymer Compression
    #[structopt(long)]
    nohpc: bool,
}

fn main() {
    let start = Instant::now();
    let opt = Opt::from_args();      
    let mut filename = PathBuf::new();
    let mut output_prefix;
    let mut k : usize = 5;
    let low_memory = opt.low_memory;
    let reference : bool = false;
    let mut use_hpc : bool = true; 
    let mut use_simd : bool = true; 
    let mut threads : usize = 8;
    let mut complement = true;
    if opt.reads.is_some() {filename = opt.reads.unwrap();} 

    if filename.as_os_str().is_empty() {panic!("Please specify an input file.");}
    let mut reads_are_fasta : bool = false;
    let filename_str = filename.to_str().unwrap();
    if filename_str.contains(".fasta.") || filename_str.ends_with(".fna") || filename_str.contains(".fna.") || filename_str.contains(".fa.") || filename_str.ends_with(".fa") || filename_str.ends_with(".fasta") { // not so robust but will have to do for now
        reads_are_fasta = true;
        eprintln!("Input file: {}", filename_str);
        eprintln!("Format: FASTA");
    }
    if opt.k.is_some() {k = opt.k.unwrap()} else {eprintln!("Warning: Using default k value ({}).", k);} 
    if opt.threads.is_some() {threads = opt.threads.unwrap();} else {eprintln!("Warning: Using default number of threads (8).");}
    output_prefix = PathBuf::from(format!("kmercamel-{}-k{}", filename_str, k));
    if opt.prefix.is_some() {output_prefix = opt.prefix.unwrap();} else {eprintln!("Warning: Using default output prefix ({}).", output_prefix.to_str().unwrap());}
    let debug = opt.debug;
    if opt.nohpc  { use_hpc = false; }
    if opt.nosimd { use_simd = false; }
    if ! std::is_x86_feature_detected!("avx512f") { 
        eprintln!("Warning: No AVX-512 CPU found, falling back to scalar implementation");
        use_simd = false; 
    }
    if use_hpc {
        if use_simd {
            eprintln!("Using HPC ntHash, with SIMD");
        }
            else {
            eprintln!("Using HPC ntHash, scalar");
        }
    } else {
        if use_simd {
            eprintln!("Using regular ntHash (not HPC), with SIMD");
        }
            else {
            eprintln!("Using regular ntHash (not HPC), scalar");
        }
    }
    let mask : i64 = (1 << (2 * k)) - 1;
    let params = Params { 
        k,
        mask,
        use_hpc,
        use_simd,
        debug,
        complement,
    };
    // init some useful objects
    // get file size for progress bar
    let _metadata = fs::metadata(&filename).expect("Error opening input file.");
    let mut queue_len = 1000; // https://doc.rust-lang.org/std/sync/mpsc/fn.sync_channel.html
                             // also: controls how many reads objects are buffered during fasta/fastq
                             // parsing
    closures::run_mers(&filename, &params, threads, queue_len, reads_are_fasta, &output_prefix);
    //println!("current time after exiting closures {:?}",Utc::now());
    let duration = start.elapsed();
    eprintln!("Total execution time: {:?}", duration);
    eprintln!("Maximum RSS: {:?}GB", (get_memory_rusage() as f32) / 1024.0 / 1024.0 / 1024.0);
}
