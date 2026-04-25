// Cargo.toml:
//   rustc-hash = "2"
//   rayon = "1"          # optional, for parallel seed-extend

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};

use indicatif::{ProgressBar, ProgressStyle};
use rustc_hash::FxHashMap;

// ============================================================
// 2-bit encoding
// ============================================================
// A = 00, C = 01, G = 10, T = 11
// Stored little-endian in u64 limbs: position 0 in bits 0..2 of limb 0, etc.
// Up to 64bp per limb. 80bp -> [u64; 2] with 48 bits unused in the second limb.
//
// We use a const generic over LIMBS so the same code handles different lengths.
// For your 80bp pool, LIMBS = 2.

const NUC_BAD: u8 = 0xFF;

#[inline]
fn nuc_to_2bit(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => NUC_BAD, // Ns, ambiguity codes -> caller decides
    }
}

#[inline]
fn revcomp_u64(key: u64, len: u32) -> u64 {
    let mut out = 0u64;
    for i in 0..len {
        let base = (key >> (2 * i)) & 0b11;
        out |= (base ^ 0b11) << (2 * (len - 1 - i));
    }
    out
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Packed<const LIMBS: usize> {
    limbs: [u64; LIMBS],
    len: u32, // number of bases actually stored
}

impl<const LIMBS: usize> Packed<LIMBS> {
    fn from_bytes(seq: &[u8]) -> Option<Self> {
        assert!(seq.len() <= 32 * LIMBS);
        let mut limbs = [0u64; LIMBS];
        for (i, &b) in seq.iter().enumerate() {
            let n = nuc_to_2bit(b);
            if n == NUC_BAD {
                return None;
            }
            limbs[i / 32] |= (n as u64) << ((i % 32) * 2);
        }
        Some(Packed {
            limbs,
            len: seq.len() as u32,
        })
    }

    fn revcomp(&self) -> Self {
        // Complement = XOR with all 11s. Reverse = reverse 2-bit groups.
        let mut out = [0u64; LIMBS];
        for i in 0..(self.len as usize) {
            let src = i / 32;
            let src_shift = (i % 32) * 2;
            let n = ((self.limbs[src] >> src_shift) & 0b11) ^ 0b11; // complement
            let dst_pos = self.len as usize - 1 - i;
            let dst = dst_pos / 32;
            let dst_shift = (dst_pos % 32) * 2;
            out[dst] |= n << dst_shift;
        }
        Packed {
            limbs: out,
            len: self.len,
        }
    }

    /// 2-bit value of the base at position i (0..len)
    #[inline]
    fn base_at(&self, i: usize) -> u64 {
        (self.limbs[i / 32] >> ((i % 32) * 2)) & 0b11
    }

    /// Compare base ranges: self[a_start..a_start+n]  ==  other[b_start..b_start+n]
    /// SIMD-able; the simple version below is already fast for 80bp.
    fn match_range<const M: usize>(
        &self,
        a_start: usize,
        other: &Packed<M>,
        b_start: usize,
        n: usize,
    ) -> bool {
        // Word-aligned fast path is possible; keep simple here.
        for k in 0..n {
            if self.base_at(a_start + k) != other.base_at(b_start + k) {
                return false;
            }
        }
        true
    }
}

// ============================================================
// Edge type and strand
// ============================================================

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Hash)]
pub enum Strand {
    Fwd,
    Rev,
}

#[inline]
fn flip(s: Strand) -> Strand {
    match s {
        Strand::Fwd => Strand::Rev,
        Strand::Rev => Strand::Fwd,
    }
}

#[derive(Debug, Clone, Copy)]
pub struct Edge {
    pub from_id: u32,
    pub from_strand: Strand,
    pub to_id: u32,
    pub to_strand: Strand,
    pub overlap_len: u32,
}

// ============================================================
// Build overlap graph
// ============================================================
//
// LIMBS: choose so that 64*LIMBS >= max sequence length.
// L_MIN: minimum overlap length (must be <= 32 for the seed-key fits-in-u64 invariant).

pub fn build_overlap_graph<const LIMBS: usize>(seqs: &[&[u8]], l_min: u32) -> Vec<Edge> {
    assert!(
        l_min >= 1 && l_min <= 32,
        "l_min must fit in a u64 seed (<=32)"
    );

    let n_seqs = seqs.len();

    // ---- 1. Pack forward + RC ---------------------------------------------
    // segment_id encoding: 0..n_seqs are forward, n_seqs..2*n_seqs are RC
    let mut packed: Vec<Packed<LIMBS>> = Vec::with_capacity(2 * n_seqs);
    for s in seqs {
        // Skip sequences with N if you want strict matching; for now error out.
        let p = Packed::<LIMBS>::from_bytes(s).expect("non-ACGT base");
        packed.push(p);
    }
    for i in 0..n_seqs {
        let rc = packed[i].revcomp();
        packed.push(rc);
    }
    let n_segs = packed.len(); // == 2 * n_seqs

    // ---- 2. Index every prefix of length L_min -----------------------------
    //
    // Key: u64 holding the L_min 2-bit values (low 2*L_min bits).
    // Value: Vec<(seq_id, strand)> — tagged so we know which orientation matched.

    let mask: u64 = if l_min == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * l_min)) - 1
    };

    let mut prefix_index: FxHashMap<u64, Vec<(u32, Strand)>> =
        FxHashMap::with_capacity_and_hasher(2 * n_seqs, Default::default());

    for i in 0..n_seqs {
        let fwd = &packed[i];
        if fwd.len < l_min {
            continue;
        }
        let fwd_key = fwd.limbs[0] & mask;
        prefix_index
            .entry(fwd_key)
            .or_default()
            .push((i as u32, Strand::Fwd));

        let rc = &packed[n_seqs + i];
        let rc_key = rc.limbs[0] & mask;
        prefix_index
            .entry(rc_key)
            .or_default()
            .push((i as u32, Strand::Rev));
    }

    // ---- 3. Seed-and-extend (single forward scan, dual rolling seeds) ------
    //
    // For each forward sequence A, at each suffix position p:
    //   fwd_key = A_fwd[p..p+l_min]          -> types 1, 2
    //   rc_key  = revcomp(A_fwd[p..p+l_min])  -> type 3 (skip type 4)

    let mut raw: Vec<Edge> = Vec::new();

    let pb = ProgressBar::new(n_seqs as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} seqs ({eta})",
        )
        .unwrap()
        .progress_chars("#>-"),
    );

    for seq_a in 0..n_seqs {
        pb.inc(1);
        let pa = &packed[seq_a];
        if pa.len < l_min {
            continue;
        }

        let pa_rc = &packed[n_seqs + seq_a];
        let len_a = pa.len;
        let last_p = (len_a - l_min) as usize;

        let mut fwd_key: u64 = pa.limbs[0] & mask;
        let mut rc_key: u64 = revcomp_u64(fwd_key, l_min);

        for p in 0..=last_p {
            // --- fwd_key lookup: types 1 and 2 ---
            if let Some(candidates) = prefix_index.get(&fwd_key) {
                let overlap = len_a - p as u32;
                for &(seq_b, strand_b) in candidates {
                    let target = if strand_b == Strand::Fwd {
                        &packed[seq_b as usize]
                    } else {
                        &packed[n_seqs + seq_b as usize]
                    };
                    if target.len < overlap {
                        continue;
                    }
                    // Filter trivial full-length self-loops (any strand combo)
                    if seq_a as u32 == seq_b && overlap == len_a {
                        continue;
                    }
                    if pa.match_range(
                        p + l_min as usize,
                        target,
                        l_min as usize,
                        (overlap - l_min) as usize,
                    ) {
                        raw.push(Edge {
                            from_id: seq_a as u32,
                            from_strand: Strand::Fwd,
                            to_id: seq_b,
                            to_strand: strand_b,
                            overlap_len: overlap,
                        });
                    }
                }
            }

            // --- rc_key lookup: type 3 only (skip Rev-tagged = type 4) ---
            if let Some(candidates) = prefix_index.get(&rc_key) {
                let overlap = p as u32 + l_min;
                for &(seq_b, strand_b) in candidates {
                    if strand_b == Strand::Rev {
                        continue; // type 4, redundant
                    }
                    let pb_fwd = &packed[seq_b as usize];
                    if pb_fwd.len < overlap {
                        continue;
                    }
                    // Filter trivial full-length self-loops
                    if seq_a as u32 == seq_b && overlap == len_a {
                        continue;
                    }
                    // Verify: pa_rc[len_a - p .. len_a] == pb_fwd[l_min .. l_min + p]
                    // When p == 0, this is n=0 (no-op, seed already verified)
                    if pa_rc.match_range(
                        (len_a - p as u32) as usize,
                        pb_fwd,
                        l_min as usize,
                        p,
                    ) {
                        raw.push(Edge {
                            from_id: seq_a as u32,
                            from_strand: Strand::Rev,
                            to_id: seq_b,
                            to_strand: Strand::Fwd,
                            overlap_len: overlap,
                        });
                    }
                }
            }

            // Roll both seeds forward
            if p < last_p {
                let new_base = pa.base_at(p + l_min as usize);
                fwd_key = (fwd_key >> 2) | (new_base << (2 * (l_min - 1)));
                rc_key = ((rc_key << 2) | (new_base ^ 0b11)) & mask;
            }
        }
    }

    pb.finish_with_message("done");

    // ---- 4. Dedup, keep max overlap per canonical edge ---------------------
    canonicalize_and_dedup(raw)
}

fn canonicalize_and_dedup(edges: Vec<Edge>) -> Vec<Edge> {
    type Key = (u32, Strand, u32, Strand);
    let mut best: FxHashMap<Key, u32> = FxHashMap::default();
    for e in edges {
        let a = (e.from_id, e.from_strand, e.to_id, e.to_strand);
        let b = (e.to_id, flip(e.to_strand), e.from_id, flip(e.from_strand));
        let key = if a <= b { a } else { b };
        let entry = best.entry(key).or_insert(0);
        if e.overlap_len > *entry {
            *entry = e.overlap_len;
        }
    }
    best.into_iter()
        .map(
            |((from_id, from_strand, to_id, to_strand), overlap_len)| Edge {
                from_id,
                from_strand,
                to_id,
                to_strand,
                overlap_len,
            },
        )
        .collect()
}

// ============================================================
// Output: GFA writer (optional convenience)
// ============================================================

pub fn write_gfa<W: std::io::Write>(
    seqs: &[&[u8]],
    edges: &[Edge],
    mut w: W,
) -> std::io::Result<()> {
    writeln!(w, "H\tVN:Z:1.0")?;
    for (i, s) in seqs.iter().enumerate() {
        writeln!(w, "S\t{}\t{}", i, std::str::from_utf8(s).unwrap())?;
    }
    for e in edges {
        let f = match e.from_strand {
            Strand::Fwd => '+',
            Strand::Rev => '-',
        };
        let t = match e.to_strand {
            Strand::Fwd => '+',
            Strand::Rev => '-',
        };
        writeln!(
            w,
            "L\t{}\t{}\t{}\t{}\t{}M",
            e.from_id, f, e.to_id, t, e.overlap_len
        )?;
    }
    Ok(())
}

const LIMBS: usize = 10;

fn parse_fasta(path: &str) -> std::io::Result<Vec<Vec<u8>>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut current: Option<Vec<u8>> = None;
    let mut skipped = 0u32;

    for line in reader.lines() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            if let Some(seq) = current.take() {
                let has_non_acgt = seq.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T'));
                if has_non_acgt {
                    skipped += 1;
                } else {
                    seqs.push(seq);
                }
            }
            current = Some(Vec::new());
        } else if let Some(ref mut seq) = current {
            seq.extend(line.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    if let Some(seq) = current {
        let has_non_acgt = seq.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T'));
        if has_non_acgt {
            skipped += 1;
        } else {
            seqs.push(seq);
        }
    }
    if skipped > 0 {
        eprintln!("skipped {} sequences with non-ACGT bases", skipped);
    }
    Ok(seqs)
}

fn usage() -> ! {
    eprintln!("Usage: oligraph-rs <input.fasta> [output.gfa] [-l <min_overlap>]");
    eprintln!("  input.fasta    Input FASTA file of sequences");
    eprintln!("  output.gfa     Output GFA file (default: stdout)");
    eprintln!("  -l <min>       Minimum overlap length (default: 20, max: 32)");
    std::process::exit(1);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        usage();
    }

    let mut fasta_path: Option<&str> = None;
    let mut gfa_path: Option<&str> = None;
    let mut l_min: u32 = 20;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-l" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("error: -l requires a value");
                    usage();
                }
                l_min = args[i].parse().unwrap_or_else(|_| {
                    eprintln!("error: invalid value for -l: {}", args[i]);
                    usage();
                });
            }
            "-" | "--" => {
                eprintln!("error: unknown flag: {}", args[i]);
                usage();
            }
            arg if arg.starts_with('-') => {
                eprintln!("error: unknown flag: {}", arg);
                usage();
            }
            _ => {
                if fasta_path.is_none() {
                    fasta_path = Some(&args[i]);
                } else if gfa_path.is_none() {
                    gfa_path = Some(&args[i]);
                } else {
                    eprintln!("error: unexpected argument: {}", args[i]);
                    usage();
                }
            }
        }
        i += 1;
    }

    let fasta_path = fasta_path.unwrap_or_else(|| {
        eprintln!("error: no input FASTA file specified");
        usage();
    });

    let seqs = parse_fasta(fasta_path).unwrap_or_else(|e| {
        eprintln!("error reading {}: {}", fasta_path, e);
        std::process::exit(1);
    });

    if seqs.is_empty() {
        eprintln!("error: no sequences found in {}", fasta_path);
        std::process::exit(1);
    }

    let max_len = seqs.iter().map(|s| s.len()).max().unwrap();
    if max_len > LIMBS * 32 {
        eprintln!(
            "error: sequence length {} exceeds capacity of LIMBS={} (max {}bp). Recompile with larger LIMBS.",
            max_len, LIMBS, LIMBS * 32
        );
        std::process::exit(1);
    }

    eprintln!(
        "loaded {} sequences (lengths {}-{}) with l_min={}",
        seqs.len(),
        seqs.iter().map(|s| s.len()).min().unwrap(),
        max_len,
        l_min
    );

    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let edges = build_overlap_graph::<LIMBS>(&seq_refs, l_min);
    eprintln!("found {} edges", edges.len());

    let result = match gfa_path {
        Some(path) => {
            let file = File::create(path).unwrap_or_else(|e| {
                eprintln!("error creating {}: {}", path, e);
                std::process::exit(1);
            });
            write_gfa(&seq_refs, &edges, BufWriter::new(file))
        }
        None => write_gfa(&seq_refs, &edges, BufWriter::new(std::io::stdout())),
    };

    if let Err(e) = result {
        eprintln!("error writing GFA: {}", e);
        std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fwd_fwd_overlap() {
        let s0: &[u8] = b"ACGTACGTACGT"; // tail "ACGTACGT"
        let s1: &[u8] = b"ACGTACGTGGGGGG"; // head "ACGTACGT"
        let edges = build_overlap_graph::<1>(&[s0, s1], 6);
        assert!(edges.iter().any(|e| e.from_id == 0
            && e.to_id == 1
            && e.from_strand == Strand::Fwd
            && e.to_strand == Strand::Fwd
            && e.overlap_len == 8));
    }

    #[test]
    fn no_internal_substring_match() {
        // S1 contains S0 as an internal substring -> NOT a valid suffix-prefix overlap
        let s0: &[u8] = b"ACGTACGT";
        let s1: &[u8] = b"GGGGACGTACGTGGGG";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6);
        // No edge between 0 and 1 in either direction
        assert!(!edges
            .iter()
            .any(|e| (e.from_id == 0 && e.to_id == 1) || (e.from_id == 1 && e.to_id == 0)));
    }

    #[test]
    fn fwd_rev_overlap() {
        let s0: &[u8] = b"AACCGGTTAACCGG";
        let s1: &[u8] = b"GGGGGGCCGGTTAA";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6);
        assert!(
            edges.iter().any(|e| e.from_id == 0
                && e.to_id == 1
                && e.from_strand == Strand::Fwd
                && e.to_strand == Strand::Rev
                && e.overlap_len == 8),
            "expected Fwd->Rev edge with overlap 8, got: {:?}",
            edges
        );
    }

    #[test]
    fn revcomp_u64_palindrome() {
        let acgt = 0b11_10_01_00u64;
        assert_eq!(revcomp_u64(acgt, 4), acgt);
    }

    #[test]
    fn revcomp_u64_non_palindrome() {
        let aaac = 0b01_00_00_00u64;
        let gttt = 0b11_11_11_10u64;
        assert_eq!(revcomp_u64(aaac, 4), gttt);
    }

    #[test]
    fn revcomp_u64_roundtrip() {
        let seq = 0b10_01_11_00u64;
        assert_eq!(revcomp_u64(revcomp_u64(seq, 4), 4), seq);
    }

    #[test]
    fn rev_fwd_overlap() {
        // Type 3: A- -> B+ means suffix of RC(A) matches prefix of B.
        // A = "AACCGGTTTT", RC(A) = "AAAACCGGTT"
        // RC(A)[4..10] = "CCGGTT", B[0..6] = "CCGGTT" -> overlap 6
        // In scan: p_fwd=0, rc_key = revcomp(A[0..6]) = "CCGGTT", matches B fwd prefix
        let a: &[u8] = b"AACCGGTTTT";
        let b: &[u8] = b"CCGGTTGGGG";
        let edges = build_overlap_graph::<1>(&[a, b], 6);
        assert!(
            edges.iter().any(|e| {
                (e.from_id == 0
                    && e.to_id == 1
                    && e.from_strand == Strand::Rev
                    && e.to_strand == Strand::Fwd
                    && e.overlap_len == 6)
                    || (e.from_id == 1
                        && e.to_id == 0
                        && e.from_strand == Strand::Rev
                        && e.to_strand == Strand::Fwd
                        && e.overlap_len == 6)
            }),
            "expected Rev->Fwd edge with overlap 6, got: {:?}",
            edges
        );
    }
}
