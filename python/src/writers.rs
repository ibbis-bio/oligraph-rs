//! Label-aware variants of the core output writers.
//!
//! The core writers (`oligraph_rs::write_gfa`, `write_fasta`, `write_contigs_fasta`)
//! identify sequences by their positional index, because `parse_fasta_str` throws
//! FASTA headers away. The Python bindings keep the record names, so these variants
//! take a `labels` slice and substitute it wherever the core writes an index.
//!
//! Passing `default_labels(n)` reproduces the core writers byte for byte — asserted
//! by the tests at the bottom of this file, so the two cannot silently drift.

use std::io::Write;

use oligraph_rs::{AssemblyMethod, Contig, Edge, Strand, Topology};

#[inline]
fn strand_char(s: Strand) -> char {
    match s {
        Strand::Fwd => '+',
        Strand::Rev => '-',
    }
}

/// `["0", "1", ..., "n-1"]` — the labels that reproduce the core writers' output.
pub fn default_labels(n: usize) -> Vec<String> {
    (0..n).map(|i| i.to_string()).collect()
}

/// GFA 1.0. Mirrors `oligraph_rs::write_gfa` (lib.rs:445).
pub fn write_gfa_labelled<W: Write>(
    seqs: &[&[u8]],
    labels: &[String],
    edges: &[Edge],
    mut w: W,
    assembly_method: AssemblyMethod,
) -> std::io::Result<()> {
    match assembly_method {
        AssemblyMethod::All => writeln!(w, "H\tVN:Z:1.0")?,
        AssemblyMethod::Pca => writeln!(w, "H\tVN:Z:1.0\tam:Z:pca")?,
    }
    for (i, s) in seqs.iter().enumerate() {
        writeln!(w, "S\t{}\t{}", labels[i], String::from_utf8_lossy(s))?;
    }
    for e in edges {
        writeln!(
            w,
            "L\t{}\t{}\t{}\t{}\t{}M",
            labels[e.from_id as usize],
            strand_char(e.from_strand),
            labels[e.to_id as usize],
            strand_char(e.to_strand),
            e.overlap_len
        )?;
    }
    Ok(())
}

/// Per-sequence FASTA with edge annotations in the header.
/// Mirrors `oligraph_rs::write_fasta` (lib.rs:476).
pub fn write_fasta_labelled<W: Write>(
    seqs: &[&[u8]],
    labels: &[String],
    edges: &[Edge],
    mut w: W,
) -> std::io::Result<()> {
    let mut edges_by_from: Vec<Vec<&Edge>> = vec![Vec::new(); seqs.len()];
    for e in edges {
        edges_by_from[e.from_id as usize].push(e);
    }
    for (i, s) in seqs.iter().enumerate() {
        write!(w, ">{}", labels[i])?;
        for e in &edges_by_from[i] {
            write!(
                w,
                " L:{}:{}:{}",
                strand_char(e.from_strand),
                labels[e.to_id as usize],
                strand_char(e.to_strand)
            )?;
        }
        writeln!(w)?;
        writeln!(w, "{}", String::from_utf8_lossy(s))?;
    }
    Ok(())
}

/// Assembled contigs as FASTA. Mirrors `oligraph_rs::write_contigs_fasta` (lib.rs:771).
pub fn write_contigs_fasta_labelled<W: Write>(
    contigs: &[Contig],
    labels: &[String],
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
            .map(|&(id, s)| format!("{}{}", labels[id as usize], strand_char(s)))
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
        writeln!(w, "{}", String::from_utf8_lossy(&c.sequence))?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use oligraph_rs::{
        LIMBS, assemble_contigs, build_overlap_graph, write_contigs_fasta, write_fasta, write_gfa,
    };

    /// Three oligos chained by 8 bp overlaps, so both edges and contigs exist.
    fn fixture() -> Vec<Vec<u8>> {
        vec![
            b"ACGTACGTACGTAAAA".to_vec(),
            b"ACGTAAAACCCCGGGG".to_vec(),
            b"CCCCGGGGTTTTATAT".to_vec(),
        ]
    }

    fn build(method: AssemblyMethod) -> (Vec<Vec<u8>>, Vec<Edge>) {
        let seqs = fixture();
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let edges = build_overlap_graph::<LIMBS>(&refs, 6, method).unwrap();
        (seqs, edges)
    }

    #[test]
    fn gfa_matches_core_writer() {
        for method in [AssemblyMethod::All, AssemblyMethod::Pca] {
            let (seqs, edges) = build(method);
            let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

            let mut core = Vec::new();
            write_gfa(&refs, &edges, &mut core, method).unwrap();

            let mut ours = Vec::new();
            write_gfa_labelled(&refs, &default_labels(refs.len()), &edges, &mut ours, method)
                .unwrap();

            assert_eq!(
                String::from_utf8(ours).unwrap(),
                String::from_utf8(core).unwrap()
            );
        }
    }

    #[test]
    fn fasta_matches_core_writer() {
        let (seqs, edges) = build(AssemblyMethod::All);
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();

        let mut core = Vec::new();
        write_fasta(&refs, &edges, &mut core).unwrap();

        let mut ours = Vec::new();
        write_fasta_labelled(&refs, &default_labels(refs.len()), &edges, &mut ours).unwrap();

        assert_eq!(
            String::from_utf8(ours).unwrap(),
            String::from_utf8(core).unwrap()
        );
    }

    #[test]
    fn contigs_fasta_matches_core_writer() {
        let (seqs, edges) = build(AssemblyMethod::All);
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let contigs = assemble_contigs(&refs, &edges).unwrap();
        assert!(!contigs.is_empty(), "fixture should assemble at least one contig");

        let mut core = Vec::new();
        write_contigs_fasta(&contigs, &mut core).unwrap();

        let mut ours = Vec::new();
        write_contigs_fasta_labelled(&contigs, &default_labels(refs.len()), &mut ours).unwrap();

        assert_eq!(
            String::from_utf8(ours).unwrap(),
            String::from_utf8(core).unwrap()
        );
    }

    #[test]
    fn labels_are_substituted() {
        let (seqs, edges) = build(AssemblyMethod::All);
        let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
        let labels: Vec<String> = vec!["oligo_a".into(), "oligo_b".into(), "oligo_c".into()];

        let mut out = Vec::new();
        write_gfa_labelled(
            &refs,
            &labels,
            &edges,
            &mut out,
            AssemblyMethod::All,
        )
        .unwrap();
        let text = String::from_utf8(out).unwrap();

        assert!(text.contains("S\toligo_a\t"));
        assert!(text.contains("S\toligo_c\t"));
        assert!(!edges.is_empty());
        for line in text.lines().filter(|l| l.starts_with('L')) {
            let fields: Vec<&str> = line.split('\t').collect();
            assert!(labels.contains(&fields[1].to_string()), "unlabelled: {line}");
            assert!(labels.contains(&fields[3].to_string()), "unlabelled: {line}");
        }
    }
}
