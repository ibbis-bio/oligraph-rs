# OliGraph

A graph-based screening tool for large oligonucleotide orders. OliGraph discovers overlapping relationships between DNA sequences, builds a bi-directed overlap graph, and assembles overlapping oligos into contigs. Designed for validating PCA (Polymerase Cycling Assembly) designs and detecting unintended cross-annealing in oligo pools.

## Installation

Requires Rust (edition 2024).

```sh
cargo build --release
```

The binary is at `target/release/oligraph-rs`.

## Usage

```
oligraph-rs <input.fasta> [output.gfa] [-l <min_overlap>] [--assembly-method <all|pca>]
```

| Argument | Description |
|---|---|
| `input.fasta` | Input FASTA file of oligonucleotide sequences |
| `output.gfa` | Output file path (omit to write GFA to stdout) |
| `-l <min>` | Minimum overlap length in bp (default: 20, max: 32) |
| `--assembly-method` | `all` (default) keeps all edge types; `pca` keeps only 3'-3' annealing overlaps |

### Example

```sh
# Screen an oligo pool, write overlap graph and assembled contigs
oligraph-rs oligos.fasta results.gfa -l 20 --assembly-method pca
```

This produces three files:

- `results.gfa` — overlap graph in GFA 1.0 format
- `results.fasta` — per-sequence FASTA with edge annotations in headers
- `results.contigs.fasta` — assembled contigs from connected components

### GFA output

```
H	VN:Z:1.0	am:Z:pca
S	0	TCACGGGGGTGGTTCCAATCTTAGTCGAG...
S	1	GGACACGGTTTGACTTACCTTTCGACACC...
L	0	+	2	-	60M
L	1	+	2	+	20M
```

Segments (`S`) are the input sequences. Links (`L`) are overlaps between sequence suffixes and prefixes, with strand orientation (`+`/`-`) reflecting forward or reverse-complement.

### Contig output

```
>contig_0 component=0 oligos=5 length=240 topology=linear branches=1 path=3+,7-,1+,0+,2-
ACGTACGT...
```

| Field | Description |
|---|---|
| `component` | Connected component ID (0 = largest) |
| `oligos` | Number of sequences in the contig |
| `length` | Assembled sequence length |
| `topology` | `linear` or `cyclic` |
| `branches` | Branch points where the greedy walk chose between multiple neighbors |
| `path` | Ordered node IDs with strand (`+`/`-`) |

## How it works

1. **2-bit packing** — sequences are encoded as 2 bits per base in `u64` limbs for fast comparison (up to 320 bp per sequence).

2. **Seed-and-extend overlap detection** — a rolling seed of length `l_min` indexes all sequence prefixes. Each suffix position is scanned against the index and verified base-by-base to find exact overlaps.

3. **Bi-directed graph model** — each sequence is a node that can be traversed in forward or reverse-complement orientation. Edges connect suffix-to-prefix overlaps across four orientations (Types 1-3), following the BCALM2 bi-directed graph convention. Mirror-symmetric edges are canonicalized and deduplicated, keeping the longest overlap per pair.

4. **Greedy contig assembly** — connected components are identified via union-find. Within each component, a bidirectional greedy walk extends from a start node, always choosing the neighbor with the longest overlap. The walk detects cyclic topology and counts branch points.

### Edge types and PCA filtering

| Type | From | To | Description |
|---|---|---|---|
| 1 | A+ | B+ | suffix(A) = prefix(B) |
| 2 | A+ | B- | suffix(A) = prefix(revcomp(B)) |
| 3 | A- | B+ | suffix(revcomp(A)) = prefix(B) |

In `--assembly-method pca` mode, only Type 2 edges are retained. These correspond to 3'-end annealing — the physical mechanism of PCA — filtering out overlaps that would not participate in assembly.

## Dependencies

- [rustc-hash](https://crates.io/crates/rustc-hash) — fast non-cryptographic hashing
- [rayon](https://crates.io/crates/rayon) — data parallelism (reserved for future use)
- [indicatif](https://crates.io/crates/indicatif) — progress bars

## License

TBD
