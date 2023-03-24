pub struct OverlapEdge {
    pub u: usize,
    pub v: usize,
    pub length: usize,
}

pub fn overlap_hamiltonian(kmers: &DashSet, params: &Params) -> Vec<OverlapEdge> {
    let k = params.k;
    let complement = params.complement;
    let size = kmers.len();
    let n = size / (1 + (complement as usize));
    let mut suffix_forbidden = vec![false; size];
    let mut prefix_forbidden = vec![false; size];
    let mut first = (0..size).collect::<Vec<usize>>();
    let mut last = first.clone();
    for d in 0..k {

    }


}