use crate::kmers;
use std::collections::{HashMap, LinkedList};

pub struct OverlapEdge {
    pub u: usize,
    pub v: usize,
    pub length: usize,
}

pub fn overlap_hamiltonian(kmers: &Vec<i64>, params: &Params) -> Vec<OverlapEdge> {
    let mut path = Vec::<OverlapEdge>::new();
    let k = params.k;
    let complement = params.complement;
    let size = kmers.len();
    let n = size / (1 + (complement as usize));
    let mut suffix_forbidden = vec![false; size];
    let mut prefix_forbidden = vec![false; size];
    let mut first = (0..size).collect::<Vec<usize>>();
    let mut last = first.clone();
    for dp in 0..k {
        let mut d = k - dp - 1;
        let mut prefixes = HashMap::<(i64, LinkedList<usize>)>::new();
        for i in 0..size {
            if prefix_forbidden[i] {continue;}
            let mut prefix = kmers::bit_prefix(kmers[i], k, d);
            if prefixes[prefix].is_empty() {prefixes[prefix].push_back(i)}
        }
        for i in 0..size {
            if suffix_forbidden[i] {continue;}
            let suffix = kmers::bit_suffix(kmers[i], d);
            if prefixes[suffix].is_empty() {continue;}
            for j in prefixes[suffix].iter_mut() {
                if (first[i] % n == j % n || first[i] % n == last[j] % n || prefix_forbidden[j]) {
                    let jp = j + 1;
                    if prefix_forbidden[j] {prefixes[suffix].remove(j);}
                    else {j += 1;}
                }
            }
            let mut new_edges : Vec<(usize, usize)> = vec![(i, j)];
            for (x, y) in new_edges.iter() {
                path.push(OverlapEdge{x, y, d});
                prefix_forbidden[y] = true;
                first[last[y]] = first[x];
                last[first[x]] = last[y];
                suffix_forbidden[x] = true;
            }
            prefixes[suffix].remove(j);
        }
        return path;
    }


}