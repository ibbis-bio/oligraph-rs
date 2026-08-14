"""OliGraph — overlap graph builder and contig assembler for oligonucleotide pools.

Find overlaps between DNA sequences, build a bi-directed overlap graph, and
assemble contigs. Built for validating PCA (Polymerase Cycling Assembly) designs
and detecting unintended cross-annealing in oligo pools.

    >>> import oligraph
    >>> g = oligraph.build({"a": "ACGTACGTACGT", "b": "ACGTACGTGGGGGG"}, min_overlap=6)
    >>> g.edges[0]
    Edge(a+ -> b+, overlap_len=8)

Sequences can come from a dict, from any iterable of ``(name, sequence)`` pairs or
bare sequences, or straight from a FASTA/FASTQ file::

    >>> g = oligraph.build_from_fasta("pool.fasta.gz")   # doctest: +SKIP

Two limits are inherited from the Rust core: sequences may be at most
``MAX_SEQUENCE_LENGTH`` (320) bases, and ``min_overlap`` must be between 1 and
``MAX_OVERLAP`` (64). Both raise ``ValueError`` rather than crashing.
"""

from __future__ import annotations

import os
from typing import Any, Iterable, List, Mapping, Tuple, Union

from . import _oligraph as _oligraph
from ._oligraph import (
    MAX_OVERLAP,
    MAX_SEQUENCE_LENGTH,
    Contig,
    Edge,
    OverlapGraph,
    Stats,
    Strand,
    Topology,
)
from ._oligraph import __version__ as __version__

__all__ = [
    "MAX_OVERLAP",
    "MAX_SEQUENCE_LENGTH",
    "Contig",
    "Edge",
    "OverlapGraph",
    "Stats",
    "Strand",
    "Topology",
    "__version__",
    "build",
    "build_from_fasta",
]

#: Anything :func:`build` will take: a name -> sequence mapping, an iterable of
#: ``(name, sequence)`` pairs, or an iterable of bare sequences.
Sequences = Union[
    Mapping[str, Any],
    Iterable[Union[str, bytes, Tuple[str, Any]]],
]

_BYTES_LIKE = (bytes, bytearray, memoryview)


def _as_text(value: Any, where: str) -> str:
    """Coerce one sequence to ``str``, rejecting anything ambiguous."""
    if isinstance(value, str):
        return value
    if isinstance(value, _BYTES_LIKE):
        return bytes(value).decode("ascii", errors="replace")
    raise TypeError(
        f"{where}: expected a str or bytes sequence, got {type(value).__name__}. "
        f"Biopython records need converting explicitly, e.g. "
        f"{{r.id: str(r.seq) for r in SeqIO.parse(path, 'fasta')}}"
    )


def _normalise(sequences: Sequences) -> Tuple[List[str], List[str]]:
    """Flatten any accepted input into parallel name and sequence lists.

    Dict order (guaranteed since Python 3.7) fixes the indices used by
    :attr:`Edge.from_index` and :attr:`Contig.path_indices`.
    """
    if isinstance(sequences, Mapping):
        names = [str(name) for name in sequences]
        seqs = [
            _as_text(seq, f"sequence {name!r}") for name, seq in sequences.items()
        ]
        return names, seqs

    if isinstance(sequences, (str, *_BYTES_LIKE)):
        raise TypeError(
            "build() takes a collection of sequences, not a single sequence. "
            "Pass a dict like {'name': 'ACGT...'}, a list of sequences, or use "
            "build_from_fasta() to read a file."
        )

    names: List[str] = []
    seqs: List[str] = []
    for i, item in enumerate(sequences):
        if isinstance(item, (str, *_BYTES_LIKE)):
            names.append(str(i))
            seqs.append(_as_text(item, f"sequence at position {i}"))
            continue
        try:
            name, seq = item
        except (TypeError, ValueError):
            raise TypeError(
                f"item at position {i} is neither a sequence nor a "
                f"(name, sequence) pair: {item!r}"
            ) from None
        names.append(str(name))
        seqs.append(_as_text(seq, f"sequence {name!r}"))

    return names, seqs


def build(
    sequences: Sequences,
    *,
    min_overlap: int = 20,
    method: str = "all",
    on_invalid: str = "skip",
) -> OverlapGraph:
    """Build an overlap graph from in-memory sequences.

    Args:
        sequences: A ``{name: sequence}`` mapping, an iterable of
            ``(name, sequence)`` pairs, or an iterable of bare sequences (which
            get named ``"0"``, ``"1"``, ...). Sequences may be ``str`` or
            ``bytes`` and are uppercased for you.
        min_overlap: Minimum overlap length in bases, 1 to ``MAX_OVERLAP`` (64).
        method: ``"all"`` keeps every edge type; ``"pca"`` keeps only 3'-end
            annealing overlaps, matching the CLI's ``-m`` flag.
        on_invalid: What to do with a record that is empty or holds a non-ACGT
            base (including ``N``). ``"skip"`` drops it and records its name in
            :attr:`OverlapGraph.skipped`; ``"raise"`` raises ``ValueError``.

    Returns:
        An :class:`OverlapGraph` whose edges and contigs carry the input names.

    Raises:
        ValueError: If ``min_overlap`` is out of range, ``method`` or
            ``on_invalid`` is not a recognised value, a sequence is longer than
            ``MAX_SEQUENCE_LENGTH`` (320) bases, no usable sequences remain, or
            ``on_invalid="raise"`` and a record is unusable.
        TypeError: If a sequence is not str- or bytes-like.

    Warns:
        UserWarning: If ``min_overlap`` is longer than the shortest sequence,
            since those sequences cannot form any edge.

    Example:
        >>> g = build({"a": "AACCGGTTAACCGG", "b": "GGGGGGCCGGTTAA"}, min_overlap=6)
        >>> g.edges[0].to_strand is Strand.Rev
        True
    """
    names, seqs = _normalise(sequences)
    return _oligraph.build_from_pairs(
        names,
        seqs,
        min_overlap=min_overlap,
        method=method,
        on_invalid=on_invalid,
    )


def build_from_fasta(
    path: Union[str, "os.PathLike[str]"],
    *,
    min_overlap: int = 20,
    method: str = "all",
    on_invalid: str = "skip",
) -> OverlapGraph:
    """Build an overlap graph from a FASTA or FASTQ file.

    Records are read with `needletail <https://crates.io/crates/needletail>`_, so
    FASTA and FASTQ, single- and multi-line, plain and gzipped all work. The
    record name is the header up to the first whitespace; any description after
    it is ignored. Duplicate names are kept as-is, and lookup by name
    (``g["dup"]``) returns the first occurrence.

    Args:
        path: Path to the file. ``.gz`` is detected from the file's contents, not
            its extension.
        min_overlap: Minimum overlap length in bases, 1 to ``MAX_OVERLAP`` (64).
        method: ``"all"`` or ``"pca"``.
        on_invalid: ``"skip"`` or ``"raise"``; see :func:`build`.

    Returns:
        An :class:`OverlapGraph` whose edges and contigs carry the record names.

    Raises:
        FileNotFoundError: If ``path`` does not exist.
        OSError: If the file cannot be read.
        ValueError: For a malformed file, a file with no records, or any of the
            conditions listed in :func:`build`.

    Warns:
        UserWarning: If ``min_overlap`` is longer than the shortest record.
    """
    return _oligraph.build_from_fasta(
        os.fspath(path),
        min_overlap=min_overlap,
        method=method,
        on_invalid=on_invalid,
    )
