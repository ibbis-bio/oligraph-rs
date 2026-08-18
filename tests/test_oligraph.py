"""Tests for the oligraph Python bindings.

Several fixtures are ports of the Rust unit tests in `oligraph-rs/src/lib.rs`, so
the expected edge counts and overlap lengths are already known-good.
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import oligraph
from oligraph import Strand, Topology

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA = Path(__file__).resolve().parent / "data"
POOL_FASTA = DATA / "pool.fasta"


# ============================================================
# Overlap detection (ports of the Rust unit tests)
# ============================================================


def test_fwd_fwd_overlap():
    """Port of `fwd_fwd_overlap` (lib.rs:846)."""
    g = oligraph.build({"a": "ACGTACGTACGT", "b": "ACGTACGTGGGGGG"}, min_overlap=6)
    assert any(
        e.from_name == "a"
        and e.to_name == "b"
        and e.from_strand is Strand.Fwd
        and e.to_strand is Strand.Fwd
        and e.overlap_len == 8
        for e in g.edges
    ), g.edges


def test_fwd_rev_overlap():
    """Port of `fwd_rev_overlap` (lib.rs:870)."""
    g = oligraph.build({"a": "AACCGGTTAACCGG", "b": "GGGGGGCCGGTTAA"}, min_overlap=6)
    assert any(
        e.from_name == "a"
        and e.to_name == "b"
        and e.from_strand is Strand.Fwd
        and e.to_strand is Strand.Rev
        and e.overlap_len == 8
        for e in g.edges
    ), g.edges


def test_no_internal_substring_match():
    """Port of `no_internal_substring_match` (lib.rs:858).

    An internal substring is not an overlap: only suffix-to-prefix counts.
    """
    g = oligraph.build({"a": "ACGTACGT", "b": "GGGGACGTACGTGGGG"}, min_overlap=6)
    assert not [e for e in g.edges if {e.from_name, e.to_name} == {"a", "b"}]


def test_lowercase_is_uppercased():
    lower = oligraph.build({"a": "acgtacgtacgt", "b": "acgtacgtgggggg"}, min_overlap=6)
    upper = oligraph.build({"a": "ACGTACGTACGT", "b": "ACGTACGTGGGGGG"}, min_overlap=6)
    assert lower.sequences == upper.sequences
    assert len(lower.edges) == len(upper.edges)


# ============================================================
# Names
# ============================================================


def test_edge_names_match_indices():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    assert g.ids == [f"oligo_{i}" for i in range(5)]
    assert g.edges, "the pool should produce edges"
    for e in g.edges:
        assert g.ids[e.from_index] == e.from_name
        assert g.ids[e.to_index] == e.to_name


def test_contig_paths_carry_names():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    assert g.contigs
    for c in g.contigs:
        assert [name for name, _ in c.path] == [g.ids[i] for i, _ in c.path_indices]
        assert [s for _, s in c.path] == [s for _, s in c.path_indices]


def test_fasta_description_is_stripped_from_name():
    """`>oligo_0 fragment 0` names the record `oligo_0`."""
    g = oligraph.build_from_fasta(POOL_FASTA)
    assert all(" " not in name for name in g.ids)


def test_lookup_by_name_and_index():
    g = oligraph.build({"a": "ACGTACGTACGT", "b": "ACGTACGTGGGGGG"}, min_overlap=6)
    assert g["a"] == "ACGTACGTACGT"
    assert g[0] == "ACGTACGTACGT"
    assert g[1] == g["b"]
    assert "a" in g
    assert "nope" not in g
    assert len(g) == 2
    with pytest.raises(KeyError):
        g["nope"]
    with pytest.raises(IndexError):
        g[99]
    with pytest.raises(TypeError):
        g[None]


def test_duplicate_names_keep_first_for_lookup():
    g = oligraph.build(
        [("dup", "ACGTACGTACGT"), ("dup", "TTTTTTTTGGGG")], min_overlap=6
    )
    assert g.ids == ["dup", "dup"]
    assert g["dup"] == "ACGTACGTACGT"


# ============================================================
# Input coercion
# ============================================================


def test_accepts_pairs_and_bare_sequences():
    seqs = ["ACGTACGTACGT", "ACGTACGTGGGGGG"]
    from_dict = oligraph.build({"0": seqs[0], "1": seqs[1]}, min_overlap=6)
    from_pairs = oligraph.build([("0", seqs[0]), ("1", seqs[1])], min_overlap=6)
    from_bare = oligraph.build(seqs, min_overlap=6)

    assert from_bare.ids == ["0", "1"]
    for other in (from_pairs, from_bare):
        assert other.ids == from_dict.ids
        assert other.sequences == from_dict.sequences
        assert other.to_gfa() == from_dict.to_gfa()


def test_accepts_bytes():
    g = oligraph.build({"a": b"ACGTACGTACGT", "b": bytearray(b"ACGTACGTGGGGGG")}, min_overlap=6)
    assert g.sequences == ["ACGTACGTACGT", "ACGTACGTGGGGGG"]


def test_generator_input():
    g = oligraph.build(
        (f"ACGTACGT{'A' * i}" for i in range(3)),
        min_overlap=6,
    )
    assert g.ids == ["0", "1", "2"]


def test_single_string_is_rejected():
    with pytest.raises(TypeError, match="collection of sequences"):
        oligraph.build("ACGTACGTACGT")


def test_non_sequence_value_is_rejected():
    with pytest.raises(TypeError, match="str or bytes"):
        oligraph.build({"a": 42})


# ============================================================
# File input
# ============================================================


def _write_variants(tmp_path: Path) -> dict:
    records = []
    name = None
    seq: list = []
    for line in POOL_FASTA.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name, seq = line[1:], []
        else:
            seq.append(line)
    records.append((name, "".join(seq)))

    gz = tmp_path / "pool.fasta.gz"
    with gzip.open(gz, "wt") as fh:
        fh.write(POOL_FASTA.read_text())

    fq = tmp_path / "pool.fastq"
    fq.write_text(
        "".join(f"@{n}\n{s}\n+\n{'I' * len(s)}\n" for n, s in records)
    )

    wrapped = tmp_path / "pool.multiline.fasta"
    wrapped.write_text(
        "".join(
            f">{n}\n" + "".join(f"{s[i:i + 30]}\n" for i in range(0, len(s), 30))
            for n, s in records
        )
    )
    return {"gz": gz, "fastq": fq, "multiline": wrapped}


def test_gzip_fastq_and_multiline_agree_with_plain_fasta(tmp_path):
    variants = _write_variants(tmp_path)
    baseline = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)

    for label, path in variants.items():
        g = oligraph.build_from_fasta(path, min_overlap=20)
        assert g.sequences == baseline.sequences, label
        assert g.to_gfa() == baseline.to_gfa(), label


def test_crlf_line_endings(tmp_path):
    """A FASTA saved on Windows has \r\n; the \r must not read as a base."""
    crlf = tmp_path / "crlf.fasta"
    crlf.write_bytes(POOL_FASTA.read_text().replace("\n", "\r\n").encode())

    g = oligraph.build_from_fasta(crlf, min_overlap=20)
    baseline = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)

    assert g.skipped == []
    assert g.ids == baseline.ids
    assert g.sequences == baseline.sequences
    assert g.to_gfa() == baseline.to_gfa()


def test_missing_file_raises_filenotfound(tmp_path):
    with pytest.raises(FileNotFoundError):
        oligraph.build_from_fasta(tmp_path / "nope.fasta")


def test_file_with_no_records(tmp_path):
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    with pytest.raises(ValueError, match="no sequence records found.*file is empty"):
        oligraph.build_from_fasta(empty)


def test_accepts_pathlike_and_str(tmp_path):
    by_path = oligraph.build_from_fasta(POOL_FASTA)
    by_str = oligraph.build_from_fasta(str(POOL_FASTA))
    assert by_path.to_gfa() == by_str.to_gfa()


# ============================================================
# Assembly method
# ============================================================


def _edge_key(e):
    return (e.from_index, e.from_strand, e.to_index, e.to_strand, e.overlap_len)


def _is_fwd_fwd(e):
    return e.from_strand is Strand.Fwd and e.to_strand is Strand.Fwd


def test_pca_drops_fwd_fwd_edges():
    """lib.rs:325 skips Fwd targets in PCA mode, keeping only 3'-end annealing.

    Checking ``to_strand`` alone would be wrong: ``canonicalize_and_dedup``
    (lib.rs:416) may emit an edge as its mirror, flipping both strands. Fwd->Fwd
    is the edge *kind* that PCA mode excludes, which is also what the web app
    asserts (web/src/main.rs:139).
    """
    seqs = {
        "a": "ACGTACGTACGT",
        "b": "ACGTACGTGGGGGG",  # a -> b anneals Fwd -> Fwd
        "c": "AACCGGTTAACCGG",
        "d": "GGGGGGCCGGTTAA",  # c -> d anneals Fwd -> Rev
    }
    all_g = oligraph.build(seqs, min_overlap=6, method="all")
    pca = oligraph.build(seqs, min_overlap=6, method="pca")

    assert [e for e in all_g.edges if _is_fwd_fwd(e)]
    assert not [e for e in pca.edges if _is_fwd_fwd(e)]
    assert {_edge_key(e) for e in pca.edges} < {_edge_key(e) for e in all_g.edges}
    assert pca.method == "pca"


def test_pca_edges_are_a_subset_on_a_real_pool():
    all_edges = {
        _edge_key(e) for e in oligraph.build_from_fasta(POOL_FASTA, method="all").edges
    }
    pca_edges = {
        _edge_key(e) for e in oligraph.build_from_fasta(POOL_FASTA, method="pca").edges
    }
    assert pca_edges
    assert pca_edges <= all_edges
    assert not [
        e
        for e in oligraph.build_from_fasta(POOL_FASTA, method="pca").edges
        if _is_fwd_fwd(e)
    ]


def test_method_is_case_insensitive():
    assert oligraph.build_from_fasta(POOL_FASTA, method="PCA").method == "pca"


def test_invalid_method():
    with pytest.raises(ValueError, match='expected "all" or "pca"'):
        oligraph.build({"a": "ACGT" * 8}, method="greedy")


# ============================================================
# Validation — these surface the core's typed errors as ValueError
# ============================================================


@pytest.mark.parametrize("min_overlap", [0, 65, 1000])
def test_min_overlap_out_of_range(min_overlap):
    with pytest.raises(ValueError, match="must be between 1 and 64"):
        oligraph.build({"a": "ACGT" * 8}, min_overlap=min_overlap)


def test_sequence_too_long():
    too_long = "ACGT" * 100  # 400 bp, over the 320 bp cap
    with pytest.raises(
        ValueError, match="400 bases, but the packing holds at most 320"
    ):
        oligraph.build({"big": too_long})


def test_limits_come_from_the_core():
    """The binding must not restate the core's limits, or they can drift apart."""
    assert oligraph.MAX_SEQUENCE_LENGTH == 320
    assert oligraph.MAX_OVERLAP == 64

    at_limit = "A" * oligraph.MAX_SEQUENCE_LENGTH
    over = at_limit + "A"
    assert oligraph.build({"a": at_limit, "b": at_limit}).stats.len_max == 320
    with pytest.raises(ValueError, match=str(oligraph.MAX_SEQUENCE_LENGTH)):
        oligraph.build({"a": over})

    oligraph.build({"a": at_limit}, min_overlap=oligraph.MAX_OVERLAP)
    with pytest.raises(ValueError):
        oligraph.build({"a": at_limit}, min_overlap=oligraph.MAX_OVERLAP + 1)


def test_max_length_is_accepted():
    at_limit = "ACGT" * (oligraph.MAX_SEQUENCE_LENGTH // 4)
    g = oligraph.build({"a": at_limit, "b": at_limit})
    assert g.stats.len_max == oligraph.MAX_SEQUENCE_LENGTH


def test_empty_input():
    with pytest.raises(ValueError, match="no sequences given"):
        oligraph.build({})


def test_non_acgt_is_skipped_by_default():
    g = oligraph.build(
        {"good": "ACGTACGTACGT", "ambiguous": "ACGTNNNNACGT", "also_good": "ACGTACGTGGGGGG"},
        min_overlap=6,
    )
    assert g.skipped == ["ambiguous"]
    assert g.ids == ["good", "also_good"]
    # Indices stay contiguous over the survivors.
    assert [e.from_index for e in g.edges] == [
        g.ids.index(e.from_name) for e in g.edges
    ]
    assert g.stats.n_skipped == 1
    assert g.stats.n_sequences == 2


def test_non_acgt_can_raise():
    with pytest.raises(ValueError, match="non-ACGT base 'N' at position 4"):
        oligraph.build(
            {"good": "ACGTACGTACGT", "ambiguous": "ACGTNNNNACGT"},
            min_overlap=6,
            on_invalid="raise",
        )


def test_empty_record_is_skipped():
    g = oligraph.build(
        {"a": "ACGTACGTACGT", "blank": "", "b": "ACGTACGTGGGGGG"}, min_overlap=6
    )
    assert g.skipped == ["blank"]


def test_all_records_unusable():
    with pytest.raises(ValueError, match="no usable sequences"):
        oligraph.build({"a": "NNNN", "b": "NNNN"})


def test_invalid_on_invalid():
    with pytest.raises(ValueError, match='expected "skip" or "raise"'):
        oligraph.build({"a": "ACGT" * 8}, on_invalid="explode")


def test_min_overlap_longer_than_shortest_sequence_warns():
    long_seq = "CGTACCATCGACAGTCAAGCTCTGTGTCGTCTAGGACGGGACGGCAGGACTACTAAGTTA"
    with pytest.warns(UserWarning, match="cannot form any edge"):
        g = oligraph.build({"short": "ACGTACGTACGT", "long": long_seq}, min_overlap=20)
    assert not [e for e in g.edges if "short" in (e.from_name, e.to_name)]


def test_no_warning_when_all_sequences_are_long_enough(recwarn):
    oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    assert [w for w in recwarn if issubclass(w.category, UserWarning)] == []


# ============================================================
# Results and output
# ============================================================


def test_stats():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    s = g.stats
    assert s.n_sequences == 5
    assert s.n_skipped == 0
    assert s.len_min == s.len_max == 60
    assert s.n_edges == len(g.edges)
    assert s.n_contigs == len(g.contigs)
    assert s.to_dict()["n_sequences"] == 5
    assert "n_isolated" in repr(s)


def test_contigs_are_sorted_longest_first():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    lengths = [len(c) for c in g.contigs]
    assert lengths == sorted(lengths, reverse=True)
    for c in g.contigs:
        assert len(c) == len(c.sequence)
        assert c.n_oligos == len(c.path)
        assert c.topology in (Topology.Linear, Topology.Cyclic)
        assert set(c.sequence) <= set("ACGT")


def test_contigs_are_cached():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    first = [c.sequence for c in g.contigs]
    second = [c.sequence for c in g.contigs]
    assert first == second


def test_edge_table_columns():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    table = g.edge_table()
    assert set(table) == {
        "from_index",
        "from_name",
        "from_strand",
        "to_index",
        "to_name",
        "to_strand",
        "overlap_len",
    }
    assert all(len(col) == len(g.edges) for col in table.values())
    assert set(table["from_strand"]) <= {"+", "-"}
    assert table["from_name"][0] == g.edges[0].from_name


def test_strand_and_topology_str():
    assert str(Strand.Fwd) == "+"
    assert str(Strand.Rev) == "-"
    assert Strand.Fwd.symbol == "+"
    assert str(Topology.Linear) == "linear"
    assert str(Topology.Cyclic) == "cyclic"


def test_strand_and_topology_are_hashable():
    """Defining `__eq__` nulls `__hash__` unless it is provided explicitly."""
    assert {Strand.Fwd, Strand.Rev, Strand.Fwd} == {Strand.Fwd, Strand.Rev}
    assert {Topology.Linear: 1}[Topology.Linear] == 1


def test_repr_is_informative():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    assert "OverlapGraph" in repr(g)
    assert "min_overlap=20" in repr(g)
    assert "->" in repr(g.edges[0])
    assert "topology=" in repr(g.contigs[0])


def test_gfa_uses_names_by_default():
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    named = g.to_gfa()
    indexed = g.to_gfa(names=False)

    assert named.startswith("H\tVN:Z:1.0\n")
    assert "S\toligo_0\t" in named
    assert "S\t0\t" in indexed
    assert "oligo_0" not in indexed
    assert len(named.splitlines()) == len(indexed.splitlines())


def test_pca_gfa_header():
    g = oligraph.build_from_fasta(POOL_FASTA, method="pca")
    assert g.to_gfa().startswith("H\tVN:Z:1.0\tam:Z:pca\n")


@pytest.mark.parametrize(
    "render, write, stem",
    [
        ("to_gfa", "write_gfa", "out.gfa"),
        ("to_fasta", "write_fasta", "out.fasta"),
        ("contigs_fasta", "write_contigs_fasta", "out.contigs.fasta"),
    ],
)
def test_write_methods_match_string_methods(tmp_path, render, write, stem):
    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20)
    path = tmp_path / stem
    getattr(g, write)(path)
    assert path.read_text() == getattr(g, render)()


# ============================================================
# Parity with the CLI
# ============================================================


def _cli_binary() -> Path | None:
    exe = ".exe" if sys.platform == "win32" else ""
    for profile in ("release", "debug"):
        candidate = REPO_ROOT / "target" / profile / f"oligraph-rs{exe}"
        if candidate.exists():
            return candidate
    found = shutil.which("oligraph-rs")
    return Path(found) if found else None


@pytest.mark.skipif(_cli_binary() is None, reason="oligraph-rs CLI not built")
def test_output_matches_cli_byte_for_byte(tmp_path):
    """`names=False` must reproduce the CLI's three output files exactly."""
    cli = _cli_binary()
    prefix = tmp_path / "cli"
    subprocess.run(
        [str(cli), "-i", str(POOL_FASTA), "-o", str(prefix), "-l", "20", "-m", "all"],
        check=True,
        capture_output=True,
    )

    g = oligraph.build_from_fasta(POOL_FASTA, min_overlap=20, method="all")
    assert g.to_gfa(names=False) == prefix.with_suffix(".gfa").read_text()
    assert g.to_fasta(names=False) == prefix.with_suffix(".fasta").read_text()
    assert (
        g.contigs_fasta(names=False)
        == prefix.with_suffix(".contigs.fasta").read_text()
    )


def test_version_is_exposed():
    assert oligraph.__version__
    assert sys.version_info >= (3, 9)
