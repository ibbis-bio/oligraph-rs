use rustc_hash::FxHashMap;

// ============================================================
// 2-bit encoding
// ============================================================
// A = 00, C = 01, G = 10, T = 11
// Stored little-endian in u64 limbs: position 0 in bits 0..2 of limb 0, etc.
// LIMBS: choose so that 32*LIMBS >= max sequence length (2 bits/base, 32 bases per u64 limb).
//
// We use a const generic over LIMBS so the same code handles different lengths.

pub const LIMBS: usize = 10;

const NUC_BAD: u8 = 0xFF;

#[inline]
fn nuc_to_2bit(b: u8) -> u8 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => NUC_BAD,
    }
}

#[inline]
fn revcomp_seed(key: u128, len: u32) -> u128 {
    let mut out = 0u128;
    for i in 0..len {
        let base = (key >> (2 * i)) & 0b11;
        out |= (base ^ 0b11) << (2 * (len - 1 - i));
    }
    out
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Packed<const LIMBS: usize> {
    limbs: [u64; LIMBS],
    len: u32,
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
        let mut out = [0u64; LIMBS];
        for i in 0..(self.len as usize) {
            let src = i / 32;
            let src_shift = (i % 32) * 2;
            let n = ((self.limbs[src] >> src_shift) & 0b11) ^ 0b11;
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

    #[inline]
    fn base_at(&self, i: usize) -> u64 {
        (self.limbs[i / 32] >> ((i % 32) * 2)) & 0b11
    }

    /// Extract up to 64 bases (128 bits) spanning the first two limbs as a u128.
    /// `mask` selects the relevant bits; e.g. for a 40-base seed use `(1u128 << 80) - 1`.
    #[inline]
    fn seed(&self, mask: u128) -> u128 {
        let lo = self.limbs[0] as u128;
        let hi = if LIMBS > 1 { self.limbs[1] as u128 } else { 0 };
        (lo | (hi << 64)) & mask
    }

    fn match_range<const M: usize>(
        &self,
        a_start: usize,
        other: &Packed<M>,
        b_start: usize,
        n: usize,
    ) -> bool {
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
pub fn flip(s: Strand) -> Strand {
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

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum AssemblyMethod {
    All,
    Pca,
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
// LIMBS: choose so that 32*LIMBS >= max sequence length (2 bits/base, 32 bases per u64 limb).
// L_MIN: minimum overlap length (must be <= 64 for the seed-key fits-in-u128 invariant).

pub fn build_overlap_graph<const LIMBS: usize>(
    seqs: &[&[u8]],
    l_min: u32,
    assembly_method: AssemblyMethod,
) -> Vec<Edge> {
    assert!(
        (1..=64).contains(&l_min),
        "l_min must fit in a u128 seed (<=64)"
    );

    let n_seqs = seqs.len();

    #[cfg(feature = "progress")]
    let t0 = std::time::Instant::now();

    // ---- 1. Pack forward + RC ---------------------------------------------
    // packed[0..n_seqs]: forward strands, packed[n_seqs..2*n_seqs]: RC (for verification)
    let mut packed: Vec<Packed<LIMBS>> = Vec::with_capacity(2 * n_seqs);
    for s in seqs {
        let p = Packed::<LIMBS>::from_bytes(s).expect("non-ACGT base");
        packed.push(p);
    }
    for i in 0..n_seqs {
        let rc = packed[i].revcomp();
        packed.push(rc);
    }

    // ---- 2. Index every prefix of length L_min -----------------------------
    //
    // Key: u128 holding the L_min 2-bit values (low 2*L_min bits).
    // Value: Vec<(seq_id, strand)> — tagged so we know which orientation matched.

    let mask: u128 = if l_min == 64 {
        u128::MAX
    } else {
        (1u128 << (2 * l_min)) - 1
    };

    let mut prefix_index: FxHashMap<u128, Vec<(u32, Strand)>> =
        FxHashMap::with_capacity_and_hasher(2 * n_seqs, Default::default());

    for i in 0..n_seqs {
        let fwd = &packed[i];
        if fwd.len < l_min {
            continue;
        }
        let fwd_key = fwd.seed(mask);
        prefix_index
            .entry(fwd_key)
            .or_default()
            .push((i as u32, Strand::Fwd));

        let rc = &packed[n_seqs + i];
        let rc_key = rc.seed(mask);
        prefix_index
            .entry(rc_key)
            .or_default()
            .push((i as u32, Strand::Rev));
    }

    #[cfg(feature = "progress")]
    {
        let t_index = std::time::Instant::now();
        let total_candidates: usize = prefix_index.values().map(|v| v.len()).sum();
        eprintln!(
            "[index]  {:.3}s  ({} unique seeds, {} indexed candidates)",
            t_index.duration_since(t0).as_secs_f64(),
            prefix_index.len(),
            total_candidates
        );
    }

    // ---- 3. Seed-and-extend (single forward scan, dual rolling seeds) ------
    //
    // For each forward sequence A, at each suffix position p:
    //   fwd_key = A_fwd[p..p+l_min]          -> types 1, 2
    //   rc_key  = revcomp(A_fwd[p..p+l_min])  -> type 3 (skip type 4)

    #[cfg(feature = "progress")]
    let pb = {
        let pb = indicatif::ProgressBar::new(n_seqs as u64);
        pb.set_style(
            indicatif::ProgressStyle::with_template(
                "{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} seqs ({eta})",
            )
            .unwrap()
            .progress_chars("#>-"),
        );
        pb
    };

    let mut raw: Vec<Edge> = Vec::new();

    for seq_a in 0..n_seqs {
        #[cfg(feature = "progress")]
        pb.inc(1);

        let pa = &packed[seq_a];
        if pa.len < l_min {
            continue;
        }

        let pa_rc = &packed[n_seqs + seq_a];
        let len_a = pa.len;
        let last_p = (len_a - l_min) as usize;

        let mut fwd_key: u128 = pa.seed(mask);
        let mut rc_key: u128 = revcomp_seed(fwd_key, l_min);

        for p in 0..=last_p {
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

            if let Some(candidates) = prefix_index.get(&rc_key) {
                let overlap = p as u32 + l_min;
                for &(seq_b, strand_b) in candidates {
                    if strand_b == Strand::Rev {
                        continue;
                    }
                    let pb_fwd = &packed[seq_b as usize];
                    if pb_fwd.len < overlap {
                        continue;
                    }
                    if seq_a as u32 == seq_b && overlap == len_a {
                        continue;
                    }
                    if pa_rc.match_range((len_a - p as u32) as usize, pb_fwd, l_min as usize, p) {
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

            if p < last_p {
                let new_base = pa.base_at(p + l_min as usize) as u128;
                fwd_key = (fwd_key >> 2) | (new_base << (2 * (l_min - 1)));
                rc_key = ((rc_key << 2) | (new_base ^ 0b11)) & mask;
            }
        }
    }

    #[cfg(feature = "progress")]
    pb.finish_with_message("done");

    #[cfg(feature = "progress")]
    {
        let t_scan = std::time::Instant::now();
        eprintln!(
            "[scan]   {:.3}s  ({} seqs, {} raw edges)",
            t_scan.duration_since(t0).as_secs_f64(),
            n_seqs,
            raw.len()
        );
    }

    // ---- 4. Dedup, keep max overlap per canonical edge ---------------------
    let raw_count = raw.len();
    let deduped = canonicalize_and_dedup(raw);

    #[cfg(feature = "progress")]
    eprintln!("[dedup]  ({} raw -> {} deduped)", raw_count, deduped.len());

    #[cfg(not(feature = "progress"))]
    let _ = raw_count;

    deduped
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
// Output writers
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

// ============================================================
// Contig assembly
// ============================================================

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

fn pick_start(comp: &[u32], adj: &[Vec<(u32, Strand, u32)>], visited: &[bool]) -> (u32, Strand) {
    let mut best_tip: Option<(u32, Strand, u32)> = None;
    for &id in comp {
        if visited[id as usize] {
            continue;
        }
        for strand in [Strand::Fwd, Strand::Rev] {
            let idx = id as usize * 2 + (strand == Strand::Rev) as usize;
            let flip_idx = id as usize * 2 + (flip(strand) == Strand::Rev) as usize;
            let has_fwd = adj[idx].iter().any(|&(nid, _, _)| !visited[nid as usize]);
            let has_flip = adj[flip_idx]
                .iter()
                .any(|&(nid, _, _)| !visited[nid as usize]);
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
                if best.is_none()
                    || ov > best.unwrap().2
                    || (ov == best.unwrap().2 && id < best.unwrap().0)
                {
                    best = Some((id, strand, ov));
                }
            }
        }
    }
    if let Some((id, strand, _)) = best {
        return (id, strand);
    }
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

    let (fwd_path, fwd_overlaps, fwd_branches) = greedy_walk(start_id, start_strand, adj, visited);

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

fn detect_topology(path: &[(u32, Strand)], adj: &[Vec<(u32, Strand, u32)>]) -> Topology {
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

fn stitch(path: &[(u32, Strand)], overlaps: &[u32], seqs: &[&[u8]]) -> Vec<u8> {
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

    let total_estimate: usize = path
        .iter()
        .map(|&(id, _)| seqs[id as usize].len())
        .sum::<usize>()
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
            let (path, overlaps, branches) = greedy_bidirectional_walk(comp, &adj, &mut visited);

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
        b.sequence
            .len()
            .cmp(&a.sequence.len())
            .then(a.component.cmp(&b.component))
    });

    contigs
}

pub fn write_contigs_fasta<W: std::io::Write>(contigs: &[Contig], mut w: W) -> std::io::Result<()> {
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

// ============================================================
// FASTA parsing (in-memory)
// ============================================================

/// Parse a FASTA-formatted string. Returns (sequences, count_skipped_non_acgt).
pub fn parse_fasta_str(content: &str) -> (Vec<Vec<u8>>, u32) {
    let mut seqs: Vec<Vec<u8>> = Vec::new();
    let mut current: Option<Vec<u8>> = None;
    let mut skipped = 0u32;

    let flush = |current: &mut Option<Vec<u8>>, seqs: &mut Vec<Vec<u8>>, skipped: &mut u32| {
        if let Some(seq) = current.take() {
            let has_non_acgt = seq.iter().any(|&b| !matches!(b, b'A' | b'C' | b'G' | b'T'));
            if has_non_acgt {
                *skipped += 1;
            } else if !seq.is_empty() {
                seqs.push(seq);
            }
        }
    };

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            flush(&mut current, &mut seqs, &mut skipped);
            current = Some(Vec::new());
        } else if let Some(ref mut seq) = current {
            seq.extend(line.as_bytes().iter().map(|b| b.to_ascii_uppercase()));
        }
    }
    flush(&mut current, &mut seqs, &mut skipped);
    (seqs, skipped)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fwd_fwd_overlap() {
        let s0: &[u8] = b"ACGTACGTACGT";
        let s1: &[u8] = b"ACGTACGTGGGGGG";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::All);
        assert!(edges.iter().any(|e| e.from_id == 0
            && e.to_id == 1
            && e.from_strand == Strand::Fwd
            && e.to_strand == Strand::Fwd
            && e.overlap_len == 8));
    }

    #[test]
    fn no_internal_substring_match() {
        let s0: &[u8] = b"ACGTACGT";
        let s1: &[u8] = b"GGGGACGTACGTGGGG";
        let edges = build_overlap_graph::<1>(&[s0, s1], 6, AssemblyMethod::All);
        assert!(
            !edges
                .iter()
                .any(|e| (e.from_id == 0 && e.to_id == 1) || (e.from_id == 1 && e.to_id == 0))
        );
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
    fn revcomp_seed_palindrome() {
        let acgt = 0b11_10_01_00u128;
        assert_eq!(revcomp_seed(acgt, 4), acgt);
    }

    #[test]
    fn revcomp_seed_non_palindrome() {
        let aaac = 0b01_00_00_00u128;
        let gttt = 0b11_11_11_10u128;
        assert_eq!(revcomp_seed(aaac, 4), gttt);
    }

    #[test]
    fn revcomp_seed_short_roundtrip() {
        let seq = 0b10_01_11_00u128;
        assert_eq!(revcomp_seed(revcomp_seed(seq, 4), 4), seq);
    }

    #[test]
    fn revcomp_seed_cross_limb() {
        // 34 A's (all zeros) -> 34 T's (all 0b11), needs 68 bits > u64
        let all_a: u128 = 0;
        let all_t: u128 = (1u128 << 68) - 1;
        assert_eq!(revcomp_seed(all_a, 34), all_t);
    }

    #[test]
    fn revcomp_seed_roundtrip_36() {
        // Non-palindromic 36-base key crossing u64 boundary
        let bases: [u128; 5] = [0, 0, 1, 2, 3];
        let mut key: u128 = 0;
        for i in 0..36u32 {
            key |= bases[i as usize % 5] << (2 * i);
        }
        assert_eq!(revcomp_seed(revcomp_seed(key, 36), 36), key);
    }

    #[test]
    fn packed_seed_single_limb() {
        let seq = b"ACGTACGTACGTACGT";
        let p = Packed::<2>::from_bytes(seq).unwrap();
        let mask: u128 = (1u128 << 32) - 1;
        let seed = p.seed(mask);
        assert_eq!(seed, p.limbs[0] as u128 & mask);
    }

    #[test]
    fn packed_seed_cross_limb() {
        // 40 bases: 32 in limb[0], 8 in limb[1]
        let seq = b"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
        let p = Packed::<2>::from_bytes(seq).unwrap();
        let mask: u128 = (1u128 << 80) - 1;
        let seed = p.seed(mask);
        assert_eq!(seed & 0b11, 0); // position 0: A=0
        assert_eq!((seed >> 2) & 0b11, 1); // position 1: C=1
        assert_eq!((seed >> 64) & 0b11, 0); // position 32: A=0 (crosses into limb[1])
        assert_eq!((seed >> 66) & 0b11, 1); // position 33: C=1
    }

    #[test]
    fn cross_limb_overlap_l_min_35() {
        // 50bp sequences with 40bp overlap, l_min=35 requires u128 seeds (70 bits)
        let s0: &[u8] = b"TTTTTTTTTTACGTCGATCGATCGATCGATCGATCGATCGATCGATCGAT";
        let s1: &[u8] = b"ACGTCGATCGATCGATCGATCGATCGATCGATCGATCGATGGGGGGGGGG";
        let edges = build_overlap_graph::<2>(&[s0, s1], 35, AssemblyMethod::All);
        assert!(
            edges.iter().any(|e| e.from_id == 0
                && e.to_id == 1
                && e.from_strand == Strand::Fwd
                && e.to_strand == Strand::Fwd
                && e.overlap_len == 40),
            "expected Fwd->Fwd edge with overlap 40 at l_min=35, got: {:?}",
            edges
        );
    }

    #[test]
    fn rev_fwd_overlap() {
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
            !edges
                .iter()
                .any(|e| e.from_strand == Strand::Fwd && e.to_strand == Strand::Fwd),
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
            assert!(
                !c.path.iter().any(|&(id, _)| id == 2),
                "isolated node 2 should not appear in any contig"
            );
        }
        assert!(
            !contigs.is_empty(),
            "should produce at least one contig from the 0-1 edge"
        );
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
        let a: &[u8] = b"AACCGGTTAACCGG";
        let b: &[u8] = b"TTTTTTCCGGTTAA";
        let c: &[u8] = b"AAAAAACCCCCCCC";
        let edges = build_overlap_graph::<1>(&[a, b, c], 6, AssemblyMethod::All);
        let contigs = assemble_contigs(&[a, b, c], &edges);
        let full_chain = contigs.iter().find(|ct| ct.path.len() >= 3);
        assert!(
            full_chain.is_some()
                || contigs.iter().any(|ct| {
                    let ids: Vec<u32> = ct.path.iter().map(|&(id, _)| id).collect();
                    ids.contains(&0) && ids.contains(&1) && ids.contains(&2)
                }),
            "expected a contig containing all three nodes, got: {:?}",
            contigs
                .iter()
                .map(|ct| {
                    let ids: Vec<u32> = ct.path.iter().map(|&(id, _)| id).collect();
                    format!(
                        "ids={:?} seq={}",
                        ids,
                        String::from_utf8_lossy(&ct.sequence)
                    )
                })
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn contig_branch_picks_longest_overlap() {
        let a: &[u8] = b"ACGTACGTACCCCCCCCC";
        let b: &[u8] = b"CCCCCCCCCTTTTTTTT";
        let c: &[u8] = b"CCCCCCGGGGGGGGGGGG";
        let seqs: &[&[u8]] = &[a, b, c];
        let edges = build_overlap_graph::<1>(seqs, 6, AssemblyMethod::All);
        let contigs = assemble_contigs(seqs, &edges);
        let longest = &contigs[0];
        assert!(
            longest.path.iter().any(|&(id, _)| id == 0),
            "longest contig should include node 0, path: {:?}",
            longest.path
        );
        assert!(
            contigs.len() >= 2,
            "expected at least 2 contigs from branching graph, got {}",
            contigs.len()
        );
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
        assert_eq!(
            contigs.len(),
            1,
            "expected 1 contig from linear chain, got {}",
            contigs.len()
        );
        let ids: Vec<u32> = contigs[0].path.iter().map(|&(id, _)| id).collect();
        assert!(
            ids.contains(&0) && ids.contains(&1) && ids.contains(&2),
            "backward walk should capture full chain, path ids: {:?}",
            ids
        );
    }

    #[test]
    fn contig_multi_component_ordering() {
        let s0: &[u8] = b"AAAAACCCCCCC";
        let s1: &[u8] = b"CCCCCCCGGGGG";
        let s2: &[u8] = b"GGGGGTTTTTTTT";
        let s3: &[u8] = b"ACACACACACAC";
        let s4: &[u8] = b"CACACACACGTGT";
        let s5: &[u8] = b"AGTCAGTCAGTC";
        let seqs: &[&[u8]] = &[s0, s1, s2, s3, s4, s5];
        let edges = build_overlap_graph::<1>(seqs, 6, AssemblyMethod::All);
        let contigs = assemble_contigs(seqs, &edges);
        assert!(
            contigs.len() >= 2,
            "expected at least 2 contigs from 2 components"
        );
        assert!(
            !contigs
                .iter()
                .any(|c| c.path.iter().any(|&(id, _)| id == 5)),
            "isolated node s5 should not appear"
        );
        let comp0_contig = contigs.iter().find(|c| c.component == 0).unwrap();
        let comp1_contig = contigs.iter().find(|c| c.component == 1).unwrap();
        assert!(
            comp0_contig.path.len() >= comp1_contig.path.len(),
            "component 0 should be the larger component"
        );
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
        assert_eq!(
            lines[0],
            ">contig_0 component=0 oligos=3 length=20 topology=linear branches=1 path=3+,7-,1+"
        );
        assert_eq!(lines[1], "ACGTACGTCCCCCCGGGGGG");
        assert_eq!(
            lines[2],
            ">contig_1 component=1 oligos=1 length=8 topology=cyclic branches=0 path=5+"
        );
        assert_eq!(lines[3], "TTTTAAAA");
    }

    #[test]
    fn contig_self_overlap_no_infinite_loop() {
        let s: &[u8] = b"ACGTACGTACGT";
        let edges = build_overlap_graph::<1>(&[s], 4, AssemblyMethod::All);
        let contigs = assemble_contigs(&[s], &edges);
        assert!(contigs.len() <= 1);
        if !contigs.is_empty() {
            assert_eq!(
                contigs[0].path.len(),
                1,
                "self-edge should not produce multi-node contig"
            );
        }
    }

    #[test]
    fn parse_fasta_str_basic() {
        let content = ">a\nACGT\n>b\nTTTT\nAAAA\n";
        let (seqs, skipped) = parse_fasta_str(content);
        assert_eq!(seqs.len(), 2);
        assert_eq!(seqs[0], b"ACGT");
        assert_eq!(seqs[1], b"TTTTAAAA");
        assert_eq!(skipped, 0);
    }

    #[test]
    fn parse_fasta_str_skips_non_acgt() {
        let content = ">a\nACGT\n>b\nACNT\n>c\nCCCC\n";
        let (seqs, skipped) = parse_fasta_str(content);
        assert_eq!(seqs.len(), 2);
        assert_eq!(skipped, 1);
    }

    #[test]
    fn parse_fasta_str_lowercase_uppercased() {
        let content = ">a\nacgt\n";
        let (seqs, _) = parse_fasta_str(content);
        assert_eq!(seqs[0], b"ACGT");
    }
}
