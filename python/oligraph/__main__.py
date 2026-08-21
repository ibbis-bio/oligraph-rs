"""Command-line interface for OliGraph.

Installed as the ``oligraph`` command, and runnable as ``python -m oligraph``.

Flags mirror the Rust CLI (``oligraph-rs``) so the two are interchangeable, with
one deliberate difference: this one labels its output with the FASTA record names
rather than positional indices. Pass ``--positional-ids`` for byte-identical
output to the Rust CLI.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from . import __version__, build_from_fasta


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="oligraph",
        description=(
            "Overlap graph builder and contig assembler for oligonucleotide pools."
        ),
    )
    p.add_argument(
        "-i",
        "--input",
        required=True,
        type=Path,
        metavar="FASTA",
        help="input FASTA/FASTQ file; plain or gzipped",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        metavar="PREFIX",
        help=(
            "output prefix; writes PREFIX.gfa, PREFIX.fasta and "
            "PREFIX.contigs.fasta. Omit to write the GFA to stdout"
        ),
    )
    # No 1-64 bound is declared here on purpose: the range lives in the Rust core
    # (validate_min_overlap) and is reported through the ValueError below, so it
    # is never restated in two places.
    p.add_argument(
        "-l",
        "--min-overlap",
        type=int,
        default=20,
        metavar="BP",
        help="minimum overlap length in bases (default: 20, range: 1-64)",
    )
    p.add_argument(
        "-m",
        "--method",
        choices=("all", "pca"),
        default="all",
        help="'all' keeps every edge type; 'pca' keeps only 3'-end annealing",
    )
    p.add_argument(
        "--on-invalid",
        choices=("skip", "raise"),
        default="skip",
        help=(
            "what to do with an empty or non-ACGT record (default: skip, which "
            "drops it and reports the name on stderr)"
        ),
    )
    p.add_argument(
        "--positional-ids",
        action="store_true",
        help="label output with positional indices instead of record names",
    )
    p.add_argument(
        "-q", "--quiet", action="store_true", help="suppress progress on stderr"
    )
    p.add_argument("--version", action="version", version=f"oligraph {__version__}")
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = _build_parser().parse_args(argv)

    def log(msg: str) -> None:
        if not args.quiet:
            print(msg, file=sys.stderr)

    # The core reports bad input as typed errors; surface them as a one-line
    # message rather than a traceback.
    try:
        g = build_from_fasta(
            args.input,
            min_overlap=args.min_overlap,
            method=args.method,
            on_invalid=args.on_invalid,
        )
    except (OSError, ValueError) as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    names = not args.positional_ids
    stats = g.stats
    log(
        f"loaded {stats.n_sequences} sequences "
        f"(lengths {stats.len_min}-{stats.len_max}) with min_overlap={g.min_overlap}"
    )
    if g.skipped:
        log(f"skipped {len(g.skipped)} unusable records: {', '.join(g.skipped)}")
    log(f"found {stats.n_edges} edges, {stats.n_isolated} isolated sequences")

    if args.output is None:
        sys.stdout.write(g.to_gfa(names=names))
        return 0

    prefix = args.output
    if prefix.parent != Path(""):
        prefix.parent.mkdir(parents=True, exist_ok=True)

    try:
        gfa = prefix.with_suffix(".gfa")
        g.write_gfa(gfa, names=names)
        log(f"wrote {gfa}")

        fasta = prefix.with_suffix(".fasta")
        g.write_fasta(fasta, names=names)
        log(f"wrote {fasta}")

        if g.contigs:
            contigs = prefix.with_suffix(".contigs.fasta")
            g.write_contigs_fasta(contigs, names=names)
            log(f"wrote {contigs} ({stats.n_contigs} contigs)")
        else:
            log("no contigs assembled (no connected components with edges)")
    except OSError as e:
        print(f"error: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
