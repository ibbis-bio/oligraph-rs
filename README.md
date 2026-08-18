# OliGraph

A graph-based screening tool for large oligonucleotide orders. OliGraph finds overlaps between DNA sequences, builds a bi-directed overlap graph, and assembles contigs. Built for validating PCA (Polymerase Cycling Assembly) designs and detecting unintended cross-annealing in oligo pools.

Available as a CLI tool, a Python package, and a browser-based web app (all computation runs locally, no server required).

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

## Python usage

### Install

```sh
pip install oligraph
```

Wheels cover Python 3.9+ on Linux, macOS and Windows and need no Rust toolchain. The same wheels are attached to each `py-v*` [release](https://github.com/AgentK9/oligraph-rs/releases), so a specific build can also be installed by URL.

To build from source instead (needs Rust 1.85 or later):

```sh
pip install git+https://github.com/AgentK9/oligraph-rs
```

### Quick start

Pass a dict of sequences, and get back a graph whose edges and contigs are labelled with your own names:

```python
import oligraph

g = oligraph.build({"oligo_0": "ACGTACGT...", "oligo_1": "..."}, min_overlap=20, method="pca")

g.edges[0]          # Edge(oligo_0+ -> oligo_1-, overlap_len=20)
g.contigs[0].path   # [('oligo_0', Strand.Fwd), ('oligo_1', Strand.Rev), ...]
g.stats             # Stats(n_sequences=5, n_edges=8, n_isolated=0, n_contigs=1, ...)
g["oligo_0"]        # the sequence, by name
```

`build()` also accepts an iterable of `(name, sequence)` pairs, or of bare sequences (named `"0"`, `"1"`, ...).

### FASTA and FASTQ files

Files are read with [needletail](https://crates.io/crates/needletail), so FASTA and FASTQ, single- and multi-line, plain and gzipped all work. The record name is the header up to the first whitespace.

```python
g = oligraph.build_from_fasta("pool.fasta.gz", min_overlap=20)
```

### Output

```python
g.to_gfa()                       # GFA 1.0, segments named after your records
g.to_gfa(names=False)            # positional IDs — byte-identical to the CLI's .gfa
g.to_fasta(), g.contigs_fasta()  # the CLI's other two outputs
g.write_gfa("results.gfa")       # stream to a file instead

import pandas as pd
pd.DataFrame(g.edge_table())     # edges as a dataframe; pandas is not a dependency
```

### Handling messy input

The core returns [typed errors](#errors) rather than panicking; the bindings translate them into ordinary Python exceptions:

| Condition | Behaviour |
|---|---|
| Record is empty or holds a non-ACGT base (including `N`) | Dropped, and its name is listed in `g.skipped`. Pass `on_invalid="raise"` for a `ValueError` |
| Sequence longer than 320 bp (`oligraph.MAX_SEQUENCE_LENGTH`) | `ValueError`; the core has to be recompiled with a larger `LIMBS` |
| `min_overlap` outside 1–64 (`oligraph.MAX_OVERLAP`) | `ValueError` |
| `min_overlap` longer than the shortest sequence | `UserWarning`; those sequences simply form no edges |

Because unusable records are dropped, `g.ids` is what maps a graph index back to your data — not the original file order.

### Development

Building the extension locally needs a Python with a shared `libpython` (the one Xcode ships on macOS does not have one):

```sh
uv venv --python 3.12 && uv pip install -e '.[test]' && .venv/bin/python -m pytest tests
```

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

### Errors

The library is total: no input makes it panic. `build_overlap_graph` and `assemble_contigs` return `Result<_, Error>`, and the writers return `io::Result`.

```rust
pub enum Error {
    MinOverlapOutOfRange { min_overlap: u32 },
    Sequence { index: usize, source: SequenceError },
    EdgeOutOfRange { id: u32, n_seqs: usize },
}

pub enum SequenceError {
    InvalidBase { position: usize, base: u8 },
    TooLong { len: usize, capacity: usize },
}
```

`Error::Sequence` reports *which* input failed and where, so callers can name the offending record rather than just the limit that was hit.

The same rules are available as standalone checks, for callers that want to filter records instead of failing the whole batch:

```rust
oligraph_rs::validate_sequence::<LIMBS>(seq)?;   // ACGT only, <= 32 * LIMBS bases
oligraph_rs::validate_min_overlap(l_min)?;       // MIN_OVERLAP ..= MAX_OVERLAP
```

The Python bindings use exactly these, which is why the 320 bp and 1–64 limits are stated in only one place.

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

### Python (`oligraph`)

- [pyo3](https://crates.io/crates/pyo3): Rust bindings for the CPython API (abi3, Python 3.9+)
- [needletail](https://crates.io/crates/needletail): FASTA/FASTQ parsing (gzip only, so no C toolchain is needed to build from source)
- [maturin](https://www.maturin.rs): build backend

### Web (`oligraph-web`)

- [leptos](https://crates.io/crates/leptos): reactive WASM UI framework
- [web-sys](https://crates.io/crates/web-sys) / [gloo](https://crates.io/crates/gloo): browser API bindings

## Licence

[MIT](LICENSE).
