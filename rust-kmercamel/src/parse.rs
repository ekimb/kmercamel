use dashmap::DashSet;
use crate::kmers::{NVAL, rc};
use crate::Params;

pub fn add_kmers(kmers: &DashSet<i64>, seq: &[u8], params: &Params) {
    let k = params.k;
    let mask = params.mask;
    let mut curr_kmer = 0;
    let mut possible_kmer_end = k;
    for i in 1..seq.len() {
        if seq[i - 1] != b'A' && seq[i - 1] != b'C' && seq[i - 1] != b'G' && seq[i - 1] != b'T' {
            possible_kmer_end = i + k;
            curr_kmer = 0;
        }
        else {
            curr_kmer <<= 2;
            curr_kmer &= mask;
            curr_kmer |= *NVAL.get(&seq[i - 1]).unwrap() as i64;
        }
        if i >= possible_kmer_end && (!params.complement || !kmers.contains(&rc(curr_kmer, k))) {
            kmers.insert(curr_kmer);
        }
    }
}