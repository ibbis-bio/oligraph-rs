# OliGraph

A graph-based screening tool for large oligonucleotide orders. OliGraph finds overlaps between DNA sequences, builds a bi-directed overlap graph, and assembles contigs. Built for validating PCA (Polymerase Cycling Assembly) designs and detecting unintended cross-annealing in oligo pools.

Available as a CLI tool and a browser-based web app (all computation runs locally, no server required).

## CLI usage

```
oligraph-rs -i <input.fasta> -o <output_prefix> [-l <min_overlap>] [-m <all|pca>]
```

| Flag | Description |
|---|---|
| `-i, --input` | Input FASTA file of oligonucleotide sequences |
| `-o, --output` | Output file prefix (writes `.gfa`, `.fasta`, `.contigs.fasta`) |
| `-l, --min-overlap` | Minimum overlap length in bp (default: 20, range: 1–64) |
| `-m, --method` | `all` (default) keeps all edge types; `pca` keeps only 3'-end annealing overlaps |

### Example

```sh
# Screen an oligo pool, write overlap graph and assembled contigs
oligraph-rs -i oligos.fasta -o results -l 20 -m pca
```

This produces three files:

- `results.gfa`: overlap graph in GFA 1.0 format
- `results.fasta`: per-sequence FASTA with edge annotations in headers
- `results.contigs.fasta`: assembled contigs from connected components

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
| `branches` | Branch points where the greedy walk chose between multiple neighbours |
| `path` | Ordered node IDs with strand (`+`/`-`) |

## How it works

1. **2-bit packing**: sequences are encoded as 2 bits per base in `u64` limbs for fast comparison (up to 320 bp per sequence).

2. **Seed-and-extend overlap detection**: a rolling seed of length `l_min` indexes all sequence prefixes. Each suffix position is scanned against the index and verified base-by-base to find exact overlaps.

3. **Bi-directed graph model**: each sequence is a node that can be traversed in forward or reverse-complement orientation. Edges connect suffix-to-prefix overlaps across three effective orientations (Fwd→Fwd, Fwd→Rev, Rev→Fwd), following the BCALM2 bi-directed graph convention. Mirror-symmetric edges are canonicalised and deduplicated, keeping the longest overlap per pair.

4. **Greedy contig assembly**: connected components are identified via union-find. Within each component, a bidirectional greedy walk extends from a start node, always choosing the neighbour with the longest overlap. The walk detects cyclic topology and counts branch points.

### Edge types and PCA filtering

Four edge kinds correspond to strand orientation pairs:

| Kind | From | To | Description |
|---|---|---|---|
| Fwd→Fwd | A+ | B+ | suffix(A) = prefix(B) |
| Fwd→Rev | A+ | B− | suffix(A) = prefix(revcomp(B)) |
| Rev→Fwd | A− | B+ | suffix(revcomp(A)) = prefix(B) |
| Rev→Rev | A− | B− | suffix(revcomp(A)) = prefix(revcomp(B)) |

In practice only the first three kinds are produced; Rev→Rev edges are excluded during overlap detection because they are mirror-symmetric with Fwd→Fwd.

In `-m pca` mode, Fwd→Fwd edges are dropped. Fwd→Rev and Rev→Fwd edges represent 3'-end annealing (the physical mechanism of PCA), so only overlaps that participate in assembly remain.

## Web app

The Leptos/WASM frontend runs entirely in the browser. No data leaves the client.

- Upload FASTA, adjust minimum overlap (1–64 bp) and assembly method in real time
- Interactive SVG graph with pan, zoom, and node dragging
- Edges colour-coded by kind with stroke width scaled by overlap length
- Component-based node colouring with bidirectional highlighting on hover
- Isolated nodes hidden by default (toggle to show)
- Contig results table with per-contig FASTA download

## Dependencies

### CLI (`oligraph-rs`)

- [rustc-hash](https://crates.io/crates/rustc-hash): fast non-cryptographic hashing
- [clap](https://crates.io/crates/clap): command-line argument parsing
- [indicatif](https://crates.io/crates/indicatif): progress bars (optional, enabled by default)

### Web (`oligraph-web`)

- [leptos](https://crates.io/crates/leptos): reactive WASM UI framework
- [web-sys](https://crates.io/crates/web-sys) / [gloo](https://crates.io/crates/gloo): browser API bindings

## Licence

TBD
