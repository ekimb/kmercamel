
use phf::{phf_map, Map};

pub static ALPH : [u8; 4] = [b'A', b'C', b'G', b'T'];

pub static COMPLEMENT: Map<u8, u8> = phf_map! {
    b'A' => b'T',
    b'C' => b'G',
    b'G' => b'C',
    b'T' => b'A',
};

pub static NVAL: Map<u8, u8> = phf_map! {
    b'A' => 0,
    b'C' => 1,
    b'G' => 2,
    b'T' => 3,
};

pub fn complement(n: u8) -> u8 {
    match n {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

pub fn bit_prefix(kmer: i64, k: usize, d: usize) -> i64 {
    kmer >> ((k - d) << 1)
}

pub fn bit_suffix(kmer: i64, k: usize, d: usize) -> i64 {
    kmer & ((1 << (d << 1)) - 1)
}

pub fn rc(kmer: i64, k: usize) -> i64 {
    let mut kmer_m = kmer;
    let mut ans : i64 = 0;
    for _ in 0..k {
        ans <<= 2;
        ans |= 3 ^ (kmer_m & 3);
        kmer_m >>= 2;
    }
    ans
}

pub fn to_kmer(enc: i64, len: usize) -> Vec<u8> {
    let mut ret = vec![b'N'; len];
    let mut enc_m = enc;
    for i in 0..len {
        ret[len - i - 1] = ALPH[(enc_m & 3) as usize];
        enc_m >>= 2;
    }
    ret
}

pub struct Kmer {
    value: Vec<u8>,
    length: usize,
}
impl Kmer {
    pub fn new(value: &[u8]) -> Self {
        Kmer { value: value.to_vec(), length: value.len() }
    }
    pub fn value(&self) -> Vec<u8> {
        self.value.to_vec()
    }
    pub fn length(&self) -> usize {
        self.length
    }
    pub fn to_number(&self) -> i64 {
        let mut ret : i64 = 0;
        for c in self.value().iter() {
            ret <<= 2;
            ret |= *NVAL.get(c).unwrap() as i64;
        }
        ret
    }
    pub fn rc(&self) -> Kmer {
        let mut rc_v = self.value().iter().map(|n| *COMPLEMENT.get(n).unwrap()).collect::<Vec<u8>>();
        rc_v.reverse();
        Kmer::new(&rc_v)
    }
}