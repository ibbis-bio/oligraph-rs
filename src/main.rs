use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};

use clap::Parser;
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

#[derive(Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
pub enum AssemblyMethod {
    All,
    Pca,
}

const HELP_TEMPLATE: &str = "\
{name} {version}
{about}

{usage-heading} {usage}

Required:
{positionals}
Options:
{options}";

#[derive(Parser)]
#[command(
    name = "oligraph",
    version,
    about = "Overlap graph builder and contig assembler for oligonucleotide pools",
    help_template = HELP_TEMPLATE,
)]
struct Cli {
    /// Input FASTA file
    #[arg(short, long)]
    input: PathBuf,

    /// Output file prefix (writes .gfa, .fasta, .contigs.fasta)
    #[arg(short, long)]
    output: PathBuf,

    /// Minimum overlap length in bp
    #[arg(short = 'l', long, default_value_t = 15)]
    min_overlap: u32,

    /// Assembly method
    #[arg(short, long, value_enum, default_value_t = AssemblyMethod::All)]
    method: AssemblyMethod,

    /// Number of threads
    #[arg(short, long, default_value_t = 1)]
    threads: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Topology {
    Linear,
    Cyclic,
}

pub struct Contig {
    pub sequence: Vec<u8>,
    pub component: usize,
    pub path: Vec<(u32, Strand)>,
    pub topology: Topology,
    pub branches: u32,
}

fn rc_bytes(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b,
        })
        .collect()
}

struct UnionFind {
    parent: Vec<u32>,
    rank: Vec<u8>,
}

impl UnionFind {
    fn new(n: usize) -> Self {
        UnionFind {
            parent: (0..n as u32).collect(),
            rank: vec![0; n],
        }
    }

    fn find(&mut self, mut x: u32) -> u32 {
        while self.parent[x as usize] != x {
            self.parent[x as usize] = self.parent[self.parent[x as usize] as usize];
            x = self.parent[x as usize];
        }
        x
    }

    fn union(&mut self, a: u32, b: u32) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra == rb {
            return;
        }
        if self.rank[ra as usize] < self.rank[rb as usize] {
            self.parent[ra as usize] = rb;
        } else if self.rank[ra as usize] > self.rank[rb as usize] {
            self.parent[rb as usize] = ra;
        } else {
            self.parent[rb as usize] = ra;
            self.rank[ra as usize] += 1;
        }
    }
}

// ============================================================
// Build overlap graph
// ============================================================
//
// LIMBS: choose so that 64*LIMBS >= max sequence length.
// L_MIN: minimum overlap length (must be <= 32 for the seed-key fits-in-u64 invariant).

pub fn build_overlap_graph<const LIMBS: usize>(seqs: &[&[u8]], l_min: u32, assembly_method: AssemblyMethod) -> Vec<Edge> {
    assert!(
        (1..=32).contains(&l_min),
        "l_min must fit in a u64 seed (<=32)"
    );

    let n_seqs = seqs.len();

    // ---- 1. Pack forward + RC ---------------------------------------------
    // packed[0..n_seqs]: forward strands, packed[n_seqs..2*n_seqs]: RC (for verification)
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
                    if assembly_method == AssemblyMethod::Pca && strand_b == Strand::Fwd {
                        continue;
                    }
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
    assembly_method: AssemblyMethod,
) -> std::io::Result<()> {
    match assembly_method {
        AssemblyMethod::All => writeln!(w, "H\tVN:Z:1.0")?,
        AssemblyMethod::Pca => writeln!(w, "H\tVN:Z:1.0\tam:Z:pca")?,
    }
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

pub fn write_fasta<W: std::io::Write>(
    seqs: &[&[u8]],
    edges: &[Edge],
    mut w: W,
) -> std::io::Result<()> {
    let mut edges_by_from: Vec<Vec<&Edge>> = vec![Vec::new(); seqs.len()];
    for e in edges {
        edges_by_from[e.from_id as usize].push(e);
    }
    for (i, s) in seqs.iter().enumerate() {
        write!(w, ">{}", i)?;
        for e in &edges_by_from[i] {
            let f = match e.from_strand {
                Strand::Fwd => '+',
                Strand::Rev => '-',
            };
            let t = match e.to_strand {
                Strand::Fwd => '+',
                Strand::Rev => '-',
            };
            write!(w, " L:{}:{}:{}", f, e.to_id, t)?;
        }
        writeln!(w)?;
        writeln!(w, "{}", std::str::from_utf8(s).unwrap())?;
    }
    Ok(())
}

fn greedy_walk(
    start_id: u32,
    start_strand: Strand,
    adj: &[Vec<(u32, Strand, u32)>],
    visited: &mut [bool],
) -> (Vec<(u32, Strand)>, Vec<u32>, u32) {
    visited[start_id as usize] = true;
    let mut path = vec![(start_id, start_strand)];
    let mut overlaps: Vec<u32> = Vec::new();
    let mut branches: u32 = 0;

    let mut cur_id = start_id;
    let mut cur_strand = start_strand;
    loop {
        let idx = cur_id as usize * 2 + (cur_strand == Strand::Rev) as usize;
        let mut picked = None;
        let mut extra_unvisited = 0u32;
        for &(nid, nstrand, ov) in &adj[idx] {
            if visited[nid as usize] {
                continue;
            }
            if picked.is_none() {
                picked = Some((nid, nstrand, ov));
            } else {
                extra_unvisited += 1;
            }
        }
        if extra_unvisited > 0 {
            branches += 1;
        }
        match picked {
            Some((nid, nstrand, ov)) => {
                visited[nid as usize] = true;
                path.push((nid, nstrand));
                overlaps.push(ov);
                cur_id = nid;
                cur_strand = nstrand;
            }
            None => break,
        }
    }
    (path, overlaps, branches)
}

fn pick_start(
    comp: &[u32],
    adj: &[Vec<(u32, Strand, u32)>],
    visited: &[bool],
) -> (u32, Strand) {
    let mut best_tip: Option<(u32, Strand, u32)> = None;
    for &id in comp {
        if visited[id as usize] {
            continue;
        }
        for strand in [Strand::Fwd, Strand::Rev] {
            let idx = id as usize * 2 + (strand == Strand::Rev) as usize;
            let flip_idx = id as usize * 2 + (flip(strand) == Strand::Rev) as usize;
            let has_fwd = adj[idx].iter().any(|&(nid, _, _)| !visited[nid as usize]);
            let has_flip = adj[flip_idx].iter().any(|&(nid, _, _)| !visited[nid as usize]);
            if has_fwd && !has_flip {
                let max_ov = adj[idx]
                    .iter()
                    .filter(|&&(nid, _, _)| !visited[nid as usize])
                    .map(|&(_, _, ov)| ov)
                    .max()
                    .unwrap_or(0);
                if best_tip.is_none()
                    || max_ov > best_tip.unwrap().2
                    || (max_ov == best_tip.unwrap().2 && id < best_tip.unwrap().0)
                {
                    best_tip = Some((id, strand, max_ov));
                }
            }
        }
    }
    if let Some((id, strand, _)) = best_tip {
        return (id, strand);
    }
    let mut best: Option<(u32, Strand, u32)> = None;
    for &id in comp {
        if visited[id as usize] {
            continue;
        }
        for strand in [Strand::Fwd, Strand::Rev] {
            let idx = id as usize * 2 + (strand == Strand::Rev) as usize;
            for &(nid, _, ov) in &adj[idx] {
                if visited[nid as usize] {
                    continue;
                }
                if best.is_none() || ov > best.unwrap().2 || (ov == best.unwrap().2 && id < best.unwrap().0) {
                    best = Some((id, strand, ov));
                }
            }
        }
    }
    if let Some((id, strand, _)) = best {
        return (id, strand);
    }
    // All unvisited nodes have only visited neighbors — pick first unvisited
    for &id in comp {
        if !visited[id as usize] {
            return (id, Strand::Fwd);
        }
    }
    unreachable!("pick_start called with all nodes visited")
}

fn greedy_bidirectional_walk(
    comp: &[u32],
    adj: &[Vec<(u32, Strand, u32)>],
    visited: &mut [bool],
) -> (Vec<(u32, Strand)>, Vec<u32>, u32) {
    let (start_id, start_strand) = pick_start(comp, adj, visited);

    let (fwd_path, fwd_overlaps, fwd_branches) =
        greedy_walk(start_id, start_strand, adj, visited);

    let (bwd_path, bwd_overlaps, bwd_branches) =
        greedy_walk(start_id, flip(start_strand), adj, visited);

    let total_branches = fwd_branches + bwd_branches;

    if bwd_path.len() <= 1 {
        return (fwd_path, fwd_overlaps, total_branches);
    }

    let mut prefix_path: Vec<(u32, Strand)> = bwd_path[1..]
        .iter()
        .rev()
        .map(|&(id, s)| (id, flip(s)))
        .collect();
    let mut prefix_overlaps: Vec<u32> = bwd_overlaps.iter().rev().copied().collect();

    prefix_path.extend(fwd_path);
    prefix_overlaps.extend(fwd_overlaps);

    (prefix_path, prefix_overlaps, total_branches)
}

fn detect_topology(
    path: &[(u32, Strand)],
    adj: &[Vec<(u32, Strand, u32)>],
) -> Topology {
    if path.len() < 2 {
        return Topology::Linear;
    }
    let (last_id, last_strand) = path[path.len() - 1];
    let (first_id, first_strand) = path[0];
    let idx = last_id as usize * 2 + (last_strand == Strand::Rev) as usize;
    for &(nid, nstrand, _) in &adj[idx] {
        if nid == first_id && nstrand == first_strand {
            return Topology::Cyclic;
        }
    }
    Topology::Linear
}

fn stitch(
    path: &[(u32, Strand)],
    overlaps: &[u32],
    seqs: &[&[u8]],
) -> Vec<u8> {
    if path.is_empty() {
        return Vec::new();
    }

    fn oriented_seq(seqs: &[&[u8]], id: u32, strand: Strand) -> Vec<u8> {
        let s = seqs[id as usize];
        match strand {
            Strand::Fwd => s.to_vec(),
            Strand::Rev => rc_bytes(s),
        }
    }

    let total_estimate: usize = path.iter().map(|&(id, _)| seqs[id as usize].len()).sum::<usize>()
        - overlaps.iter().map(|&o| o as usize).sum::<usize>();
    let mut buf: Vec<u8> = Vec::with_capacity(total_estimate);

    let (first_id, first_strand) = path[0];
    buf.extend_from_slice(&oriented_seq(seqs, first_id, first_strand));

    for (i, &(id, strand)) in path[1..].iter().enumerate() {
        let ov = overlaps[i] as usize;
        let s = oriented_seq(seqs, id, strand);
        buf.extend_from_slice(&s[ov..]);
    }

    buf
}

pub fn assemble_contigs(seqs: &[&[u8]], edges: &[Edge]) -> Vec<Contig> {
    let n = seqs.len();
    if n == 0 || edges.is_empty() {
        return Vec::new();
    }

    let mut adj: Vec<Vec<(u32, Strand, u32)>> = vec![Vec::new(); n * 2];
    for e in edges {
        let from_idx = e.from_id as usize * 2 + (e.from_strand == Strand::Rev) as usize;
        adj[from_idx].push((e.to_id, e.to_strand, e.overlap_len));
        let mirror_from = e.to_id as usize * 2 + (flip(e.to_strand) == Strand::Rev) as usize;
        adj[mirror_from].push((e.from_id, flip(e.from_strand), e.overlap_len));
    }
    for slot in &mut adj {
        slot.sort_unstable_by(|a, b| {
            b.2.cmp(&a.2)
                .then(a.0.cmp(&b.0))
                .then((a.1 as u8).cmp(&(b.1 as u8)))
        });
    }

    let mut uf = UnionFind::new(n);
    let mut has_edge = vec![false; n];
    for e in edges {
        if e.from_id != e.to_id {
            uf.union(e.from_id, e.to_id);
            has_edge[e.from_id as usize] = true;
            has_edge[e.to_id as usize] = true;
        }
    }

    let mut comp_map: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
    for i in 0..n {
        if !has_edge[i] {
            continue;
        }
        let root = uf.find(i as u32);
        comp_map.entry(root).or_default().push(i as u32);
    }
    let mut components: Vec<Vec<u32>> = comp_map.into_values().collect();
    components.sort_unstable_by(|a, b| {
        b.len().cmp(&a.len()).then_with(|| {
            let min_a = a.iter().min().unwrap();
            let min_b = b.iter().min().unwrap();
            min_a.cmp(min_b)
        })
    });

    let mut visited = vec![false; n];
    let mut contigs: Vec<Contig> = Vec::new();

    for (comp_idx, comp) in components.iter().enumerate() {
        while comp.iter().any(|&id| !visited[id as usize]) {
            let (path, overlaps, branches) =
                greedy_bidirectional_walk(comp, &adj, &mut visited);

            let topology = detect_topology(&path, &adj);
            let sequence = stitch(&path, &overlaps, seqs);

            contigs.push(Contig {
                sequence,
                component: comp_idx,
                path,
                topology,
                branches,
            });
        }
    }

    contigs.sort_unstable_by(|a, b| {
        b.sequence.len().cmp(&a.sequence.len()).then(a.component.cmp(&b.component))
    });

    contigs
}

pub fn write_contigs_fasta<W: std::io::Write>(
    contigs: &[Contig],
    mut w: W,
) -> std::io::Result<()> {
    for (i, c) in contigs.iter().enumerate() {
        let topo = match c.topology {
            Topology::Linear => "linear",
            Topology::Cyclic => "cyclic",
        };
        let path_str: String = c
            .path
            .iter()
            .map(|&(id, s)| {
                let sc = match s {
                    Strand::Fwd => '+',
                    Strand::Rev => '-',
                };
                format!("{}{}", id, sc)
            })
            .collect::<Vec<_>>()
            .join(",");
        writeln!(
            w,
            ">contig_{} component={} oligos={} length={} topology={} branches={} path={}",
            i,
            c.component,
            c.path.len(),
            c.sequence.len(),
            topo,
            c.branches,
            path_str,
        )?;
        writeln!(w, "{}", std::str::from_utf8(&c.sequence).unwrap())?;
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
    eprintln!("Usage: oligraph-rs <input.fasta> [output.gfa] [-l <min_overlap>] [--assembly-method <all|pca>]");
    eprintln!("  input.fasta                Input FASTA file of sequences");
    eprintln!("  output.gfa                 Output GFA file (default: stdout, GFA only)");
    eprintln!("                             Also writes .fasta alongside the .gfa");
    eprintln!("  -l <min>                   Minimum overlap length (default: 20, max: 32)");
    eprintln!("  --assembly-method <method>  Assembly method filter (default: all)");
    std::process::exit(1);
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        usage();
    }

    let mut fasta_path: Option<&str> = None;
    let mut gfa_path: Option<&str> = None;
    let mut l_min: u32 = 20;
    let mut assembly_method = AssemblyMethod::All;

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
            "--assembly-method" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("error: --assembly-method requires a value");
                    usage();
                }
                assembly_method = match args[i].as_str() {
                    "all" => AssemblyMethod::All,
                    "pca" => AssemblyMethod::Pca,
                    _ => {
                        eprintln!("error: unknown assembly method '{}' (expected: all, pca)", args[i]);
                        usage();
                    }
                };
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
        "loaded {} sequences (lengths {}-{}) with l_min={}, assembly_method={}",
        seqs.len(),
        seqs.iter().map(|s| s.len()).min().unwrap(),
        max_len,
        l_min,
        match assembly_method {
            AssemblyMethod::All => "all",
            AssemblyMethod::Pca => "pca",
        }
    );

    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let edges = build_overlap_graph::<LIMBS>(&seq_refs, l_min, assembly_method);
    eprintln!("found {} edges", edges.len());

    match gfa_path {
        Some(path) => {
            let gfa_file = File::create(path).unwrap_or_else(|e| {
                eprintln!("error creating {}: {}", path, e);
                std::process::exit(1);
            });
            if let Err(e) = write_gfa(&seq_refs, &edges, BufWriter::new(gfa_file), assembly_method)
            {
                eprintln!("error writing GFA: {}", e);
                std::process::exit(1);
            }

            let fasta_path = Path::new(path).with_extension("fasta");
            let fasta_file = File::create(&fasta_path).unwrap_or_else(|e| {
                eprintln!("error creating {}: {}", fasta_path.display(), e);
                std::process::exit(1);
            });
            if let Err(e) = write_fasta(&seq_refs, &edges, BufWriter::new(fasta_file)) {
                eprintln!("error writing FASTA: {}", e);
                std::process::exit(1);
            }
            eprintln!("wrote {}", fasta_path.display());

            let contigs = assemble_contigs(&seq_refs, &edges);
            if !contigs.is_empty() {
                let contigs_path = Path::new(path).with_extension("contigs.fasta");
                let contigs_file = File::create(&contigs_path).unwrap_or_else(|e| {
                    eprintln!("error creating {}: {}", contigs_path.display(), e);
                    std::process::exit(1);
                });
                if let Err(e) = write_contigs_fasta(&contigs, BufWriter::new(contigs_file)) {
                    eprintln!("error writing contigs: {}", e);
                    std::process::exit(1);
                }
                eprintln!(
                    "wrote {} ({} contigs)",
                    contigs_path.display(),
                    contigs.len()
                );
            } else {
                eprintln!("no contigs assembled (no connected components with edges)");
            }
        }
        None => {
            if let Err(e) =
                write_gfa(&seq_refs, &edges, BufWriter::new(std::io::stdout()), assembly_method)
            {
                eprintln!("error writing GFA: {}", e);
                std::process::exit(1);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fwd_fwd_overlap() {
        let s0: &[u8] = b"ACGTACGTACGT"; // tail "ACGTACGT"
        let s1: &[u8] = b"ACGTACGTGGGGGG"; // head "ACGTACGT"
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::All);
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
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::All);
        // No edge between 0 and 1 in either direction
        assert!(!edges
            .iter()
            .any(|e| (e.from_id == 0 && e.to_id == 1) || (e.from_id == 1 && e.to_id == 0)));
    }

    #[test]
    fn fwd_rev_overlap() {
        let s0: &[u8] = b"AACCGGTTAACCGG";
        let s1: &[u8] = b"GGGGGGCCGGTTAA";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::All);
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
        let edges = build_overlap_graph::<1>(&[a, b], 6, AssemblyMethod::All);
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

    #[test]
    fn self_overlap_tandem_repeat() {
        let s: &[u8] = b"ACGTACGTACGT";
        let edges = build_overlap_graph::<1>(&[s], 4, AssemblyMethod::All);
        assert!(
            edges
                .iter()
                .any(|e| e.from_id == 0 && e.to_id == 0 && e.overlap_len == 8),
            "expected tandem self-overlap with overlap 8, got: {:?}",
            edges
        );
        assert!(
            !edges
                .iter()
                .any(|e| e.from_id == 0 && e.to_id == 0 && e.overlap_len == 12),
            "trivial full-length self-loop should be filtered, got: {:?}",
            edges
        );
    }

    #[test]
    fn self_rc_overlap() {
        // s = "AACCAACCGGTT" (12bp)
        // RC = comp: TTGGTTGGCCAA, rev: AACCGGTTGGTT (12bp)
        // s suffix at p=4: "AACCGGTT" (8bp). RC prefix[0..8]: "AACCGGTT". MATCH!
        // Edge: 0+ -> 0- with overlap 8.
        let s: &[u8] = b"AACCAACCGGTT";
        let edges = build_overlap_graph::<1>(&[s], 6, AssemblyMethod::All);
        assert!(
            edges.iter().any(|e| e.from_id == 0
                && e.to_id == 0
                && e.from_strand == Strand::Fwd
                && e.to_strand == Strand::Rev
                && e.overlap_len == 8),
            "expected self-RC edge 0+ -> 0- with overlap 8, got: {:?}",
            edges
        );
    }

    #[test]
    fn overlap_length_correctness() {
        let s0: &[u8] = b"ACGTACGTCCCC";
        let s1: &[u8] = b"ACGTCCCCGGGG";
        let s2: &[u8] = b"CCCCGGGGACGT";
        let edges = build_overlap_graph::<1>(&[s0, s1, s2], 4, AssemblyMethod::All);

        let has_edge = |from: u32, fs: Strand, to: u32, ts: Strand, ov: u32| -> bool {
            edges.iter().any(|e| {
                (e.from_id == from
                    && e.from_strand == fs
                    && e.to_id == to
                    && e.to_strand == ts
                    && e.overlap_len == ov)
                    || (e.from_id == to
                        && e.from_strand == flip(ts)
                        && e.to_id == from
                        && e.to_strand == flip(fs)
                        && e.overlap_len == ov)
            })
        };

        assert!(
            has_edge(0, Strand::Fwd, 1, Strand::Fwd, 8),
            "expected 0+ -> 1+ overlap 8, got: {:?}",
            edges
        );
        assert!(
            has_edge(1, Strand::Fwd, 2, Strand::Fwd, 8),
            "expected 1+ -> 2+ overlap 8, got: {:?}",
            edges
        );
        assert!(
            has_edge(0, Strand::Fwd, 2, Strand::Fwd, 4),
            "expected 0+ -> 2+ overlap 4, got: {:?}",
            edges
        );
        assert!(
            has_edge(1, Strand::Fwd, 0, Strand::Rev, 4),
            "expected 1+ -> 0- overlap 4, got: {:?}",
            edges
        );
    }

    #[test]
    fn pca_drops_fwd_fwd_overlap() {
        let s0: &[u8] = b"ACGTACGTACGT";
        let s1: &[u8] = b"ACGTACGTGGGGGG";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::Pca);
        assert!(
            !edges.iter().any(|e| e.from_strand == Strand::Fwd && e.to_strand == Strand::Fwd),
            "PCA mode should drop Type 1 (Fwd->Fwd) edges, got: {:?}",
            edges
        );
    }

    #[test]
    fn pca_keeps_fwd_rev_overlap() {
        let s0: &[u8] = b"AACCGGTTAACCGG";
        let s1: &[u8] = b"GGGGGGCCGGTTAA";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::Pca);
        assert!(
            edges.iter().any(|e| e.from_id == 0
                && e.to_id == 1
                && e.from_strand == Strand::Fwd
                && e.to_strand == Strand::Rev
                && e.overlap_len == 8),
            "PCA mode should keep Type 2 (Fwd->Rev) edges, got: {:?}",
            edges
        );
    }

    #[test]
    fn pca_keeps_rev_fwd_overlap() {
        let a: &[u8] = b"AACCGGTTTT";
        let b: &[u8] = b"CCGGTTGGGG";
        let edges = build_overlap_graph::<1>(&[a, b], 6, AssemblyMethod::Pca);
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
            "PCA mode should keep Type 3 (Rev->Fwd) edges, got: {:?}",
            edges
        );
    }

    #[test]
    fn gfa_header_includes_assembly_method_tag() {
        let s0: &[u8] = b"ACGTACGT";
        let seqs: &[&[u8]] = &[s0];
        let edges: &[Edge] = &[];

        let mut buf_all = Vec::new();
        write_gfa(seqs, edges, &mut buf_all, AssemblyMethod::All).unwrap();
        let header_all = std::str::from_utf8(&buf_all)
            .unwrap()
            .lines()
            .next()
            .unwrap();
        assert_eq!(header_all, "H\tVN:Z:1.0");

        let mut buf_pca = Vec::new();
        write_gfa(seqs, edges, &mut buf_pca, AssemblyMethod::Pca).unwrap();
        let header_pca = std::str::from_utf8(&buf_pca)
            .unwrap()
            .lines()
            .next()
            .unwrap();
        assert_eq!(header_pca, "H\tVN:Z:1.0\tam:Z:pca");
    }

    #[test]
    fn fasta_output_format() {
        let seqs: Vec<&[u8]> = vec![b"ACGTACGT", b"CGTACGTT", b"TTTTAAAA"];
        let edges = vec![
            Edge {
                from_id: 0,
                from_strand: Strand::Fwd,
                to_id: 1,
                to_strand: Strand::Fwd,
                overlap_len: 7,
            },
            Edge {
                from_id: 0,
                from_strand: Strand::Rev,
                to_id: 2,
                to_strand: Strand::Fwd,
                overlap_len: 4,
            },
        ];
        let mut buf = Vec::new();
        write_fasta(&seqs, &edges, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[0], ">0 L:+:1:+ L:-:2:+");
        assert_eq!(lines[1], "ACGTACGT");
        assert_eq!(lines[2], ">1");
        assert_eq!(lines[3], "CGTACGTT");
        assert_eq!(lines[4], ">2");
        assert_eq!(lines[5], "TTTTAAAA");
        assert_eq!(lines.len(), 6);
    }

    #[test]
    fn assemble_contigs_isolated_nodes_discarded() {
        let s0: &[u8] = b"ACGTACGTCCCCCC";
        let s1: &[u8] = b"CCCCCCGGGGGGGG";
        let s2: &[u8] = b"TTTTTTTTTTTTTT";
        let edges = build_overlap_graph::<1>(&[s0, s1, s2], 6, AssemblyMethod::All);
        let contigs = assemble_contigs(&[s0, s1, s2], &edges);
        for c in &contigs {
            assert!(!c.path.iter().any(|&(id, _)| id == 2),
                "isolated node 2 should not appear in any contig");
        }
        assert!(!contigs.is_empty(), "should produce at least one contig from the 0-1 edge");
    }

    #[test]
    fn union_find_basic() {
        let mut uf = UnionFind::new(5);
        uf.union(0, 1);
        uf.union(2, 3);
        uf.union(1, 3);
        assert_eq!(uf.find(0), uf.find(3));
        assert_ne!(uf.find(0), uf.find(4));
    }

    #[test]
    fn contig_single_edge_stitch() {
        let s0: &[u8] = b"ACGTACGTCCCCCC";
        let s1: &[u8] = b"CCCCCCGGGGGGGG";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::All);
        let contigs = assemble_contigs(&[s0, s1], &edges);
        assert_eq!(contigs.len(), 1);
        let c = &contigs[0];
        assert_eq!(c.path.len(), 2);
        // Assembler may walk in either direction; accept both orientations
        let fwd = b"ACGTACGTCCCCCCGGGGGGGG".to_vec();
        let rev = rc_bytes(&fwd);
        assert!(
            c.sequence == fwd || c.sequence == rev,
            "expected stitched sequence in either orientation, got: {}",
            String::from_utf8_lossy(&c.sequence)
        );
        assert_eq!(c.topology, Topology::Linear);
    }

    #[test]
    fn contig_linear_chain_three_nodes() {
        // A = "AACCGGTTAACCGG" (14bp), B = "TTTTTTCCGGTTAA" (14bp), C = "AAAAAACCCCCCCC"
        // The assembler produces a 3-node chain through A, B, C (possibly in RC orientation).
        // The chain stitches all three nodes regardless of walk direction.
        let a: &[u8] = b"AACCGGTTAACCGG";
        let b: &[u8] = b"TTTTTTCCGGTTAA";
        let c: &[u8] = b"AAAAAACCCCCCCC";
        let edges = build_overlap_graph::<1>(&[a, b, c], 6, AssemblyMethod::All);
        let contigs = assemble_contigs(&[a, b, c], &edges);
        // Expect a contig that includes all three nodes in its path
        let full_chain = contigs.iter().find(|ct| ct.path.len() >= 3);
        assert!(
            full_chain.is_some() || contigs.iter().any(|ct| {
                let ids: Vec<u32> = ct.path.iter().map(|&(id, _)| id).collect();
                ids.contains(&0) && ids.contains(&1) && ids.contains(&2)
            }),
            "expected a contig containing all three nodes, got: {:?}",
            contigs.iter().map(|ct| {
                let ids: Vec<u32> = ct.path.iter().map(|&(id, _)| id).collect();
                format!("ids={:?} seq={}", ids, String::from_utf8_lossy(&ct.sequence))
            }).collect::<Vec<_>>()
        );
    }

    #[test]
    fn contig_branch_picks_longest_overlap() {
        // A has overlaps to both B and C; greedy should pick the longer overlap first.
        // The first (longest) contig should contain A and its best neighbor.
        let a: &[u8] = b"ACGTACGTACCCCCCCCC";
        let b: &[u8] = b"CCCCCCCCCTTTTTTTT";
        let c: &[u8] = b"CCCCCCGGGGGGGGGGGG";
        let seqs: &[&[u8]] = &[a, b, c];
        let edges = build_overlap_graph::<1>(seqs, 6, AssemblyMethod::All);
        let contigs = assemble_contigs(seqs, &edges);
        // The first contig (longest sequence) should contain node 0 (A)
        let longest = &contigs[0];
        assert!(
            longest.path.iter().any(|&(id, _)| id == 0),
            "longest contig should include node 0, path: {:?}",
            longest.path
        );
        // At least 2 contigs expected since B and C can't both be on the same path with A
        // (one gets consumed first, the other becomes a separate contig)
        assert!(contigs.len() >= 2, "expected at least 2 contigs from branching graph, got {}", contigs.len());
    }

    #[test]
    fn contig_equal_overlap_tiebreak_by_id() {
        let a: &[u8] = b"ACGTACCCCCCC";
        let b: &[u8] = b"CCCCCCTTTTTT";
        let c: &[u8] = b"CCCCCCAAAAAA";
        let seqs: &[&[u8]] = &[a, b, c];
        let edges = build_overlap_graph::<1>(seqs, 6, AssemblyMethod::All);
        let contigs = assemble_contigs(seqs, &edges);
        let first = &contigs[0];
        let ids: Vec<u32> = first.path.iter().map(|&(id, _)| id).collect();
        assert!(
            ids.contains(&0) && ids.contains(&1),
            "equal overlap tiebreaker should pick lower seq_id (1 over 2), path ids: {:?}",
            ids
        );
    }

    #[test]
    fn contig_empty_graph() {
        let s0: &[u8] = b"ACGTACGT";
        let contigs = assemble_contigs(&[s0], &[]);
        assert!(contigs.is_empty());
    }

    #[test]
    fn contig_backward_walk_extends() {
        let a: &[u8] = b"AAAACCCCCCCC";
        let b: &[u8] = b"CCCCCCCCTTTT";
        let c: &[u8] = b"TTTTGGGGGGGG";
        let seqs: &[&[u8]] = &[a, b, c];
        let edges = build_overlap_graph::<1>(seqs, 4, AssemblyMethod::All);
        let contigs = assemble_contigs(seqs, &edges);
        assert_eq!(contigs.len(), 1, "expected 1 contig from linear chain, got {}", contigs.len());
        let ids: Vec<u32> = contigs[0].path.iter().map(|&(id, _)| id).collect();
        assert!(ids.contains(&0) && ids.contains(&1) && ids.contains(&2),
            "backward walk should capture full chain, path ids: {:?}", ids);
    }

    #[test]
    fn contig_multi_component_ordering() {
        // Two disconnected components: {s0,s1,s2} and {s3,s4}, plus isolated s5.
        // Use a sequence for s5 that won't overlap with anything at l_min=6.
        let s0: &[u8] = b"AAAAACCCCCCC";
        let s1: &[u8] = b"CCCCCCCGGGGG";
        let s2: &[u8] = b"GGGGGTTTTTTTT";
        let s3: &[u8] = b"ACACACACACAC";
        let s4: &[u8] = b"CACACACACGTGT";
        let s5: &[u8] = b"AGTCAGTCAGTC";
        let seqs: &[&[u8]] = &[s0, s1, s2, s3, s4, s5];
        let edges = build_overlap_graph::<1>(seqs, 6, AssemblyMethod::All);
        let contigs = assemble_contigs(seqs, &edges);
        assert!(contigs.len() >= 2, "expected at least 2 contigs from 2 components");
        assert!(!contigs.iter().any(|c| c.path.iter().any(|&(id, _)| id == 5)),
            "isolated node s5 should not appear");
        let comp0_contig = contigs.iter().find(|c| c.component == 0).unwrap();
        let comp1_contig = contigs.iter().find(|c| c.component == 1).unwrap();
        assert!(comp0_contig.path.len() >= comp1_contig.path.len(),
            "component 0 should be the larger component");
    }

    #[test]
    fn write_contigs_fasta_format() {
        let contigs = vec![
            Contig {
                sequence: b"ACGTACGTCCCCCCGGGGGG".to_vec(),
                component: 0,
                path: vec![(3, Strand::Fwd), (7, Strand::Rev), (1, Strand::Fwd)],
                topology: Topology::Linear,
                branches: 1,
            },
            Contig {
                sequence: b"TTTTAAAA".to_vec(),
                component: 1,
                path: vec![(5, Strand::Fwd)],
                topology: Topology::Cyclic,
                branches: 0,
            },
        ];
        let mut buf = Vec::new();
        write_contigs_fasta(&contigs, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = output.lines().collect();
        assert_eq!(lines[0], ">contig_0 component=0 oligos=3 length=20 topology=linear branches=1 path=3+,7-,1+");
        assert_eq!(lines[1], "ACGTACGTCCCCCCGGGGGG");
        assert_eq!(lines[2], ">contig_1 component=1 oligos=1 length=8 topology=cyclic branches=0 path=5+");
        assert_eq!(lines[3], "TTTTAAAA");
    }

    #[test]
    fn contig_self_overlap_no_infinite_loop() {
        let s: &[u8] = b"ACGTACGTACGT";
        let edges = build_overlap_graph::<1>(&[s], 4, AssemblyMethod::All);
        let contigs = assemble_contigs(&[s], &edges);
        assert!(contigs.len() <= 1);
        if !contigs.is_empty() {
            assert_eq!(contigs[0].path.len(), 1, "self-edge should not produce multi-node contig");
        }
    }
}
