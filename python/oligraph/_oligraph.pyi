"""Type stubs for the compiled Rust extension.

The public API is re-exported from ``oligraph/__init__.py``; import from
``oligraph`` rather than from here.
"""

from enum import Enum
from typing import Any, Dict, List, Sequence, Tuple, Union
import os

__version__: str

#: Longest sequence the core's bit packing holds: 10 limbs x 32 bases.
MAX_SEQUENCE_LENGTH: int
#: Largest permitted ``min_overlap``: the seed key must fit in a u128.
MAX_OVERLAP: int

class Strand(Enum):
    Fwd = ...
    Rev = ...

    @property
    def symbol(self) -> str:
        """``"+"`` for :attr:`Fwd`, ``"-"`` for :attr:`Rev`."""

class Topology(Enum):
    Linear = ...
    Cyclic = ...

class Edge:
    """One overlap between two oligos."""

    @property
    def from_index(self) -> int:
        """Position of the source sequence in :attr:`OverlapGraph.ids`."""

    @property
    def from_name(self) -> str: ...
    @property
    def from_strand(self) -> Strand: ...
    @property
    def to_index(self) -> int: ...
    @property
    def to_name(self) -> str: ...
    @property
    def to_strand(self) -> Strand: ...
    @property
    def overlap_len(self) -> int:
        """Length of the overlap in bases."""

class Contig:
    """A contig walked out of one connected component of the overlap graph."""

    @property
    def sequence(self) -> str: ...
    @property
    def component(self) -> int:
        """Connected-component id; 0 is the largest component."""

    @property
    def path(self) -> List[Tuple[str, Strand]]:
        """The oligos walked, in order, as ``(name, strand)``."""

    @property
    def path_indices(self) -> List[Tuple[int, Strand]]:
        """The same walk as ``(index, strand)``, indexing into :attr:`OverlapGraph.ids`."""

    @property
    def topology(self) -> Topology: ...
    @property
    def branches(self) -> int:
        """Branch points where the greedy walk chose between several neighbours."""

    @property
    def n_oligos(self) -> int: ...
    def __len__(self) -> int:
        """Length of the assembled sequence in bases."""

class Stats:
    """Summary counts for a graph."""

    @property
    def n_sequences(self) -> int: ...
    @property
    def n_skipped(self) -> int: ...
    @property
    def len_min(self) -> int: ...
    @property
    def len_max(self) -> int: ...
    @property
    def n_edges(self) -> int: ...
    @property
    def n_isolated(self) -> int:
        """Sequences with no overlap to any other sequence."""

    @property
    def n_contigs(self) -> int: ...
    def to_dict(self) -> Dict[str, int]: ...

class OverlapGraph:
    """The overlap graph for a set of oligos, with names preserved."""

    @property
    def ids(self) -> List[str]:
        """Sequence names, in the order the graph indexes them."""

    @property
    def sequences(self) -> List[str]:
        """Sequences, uppercased, parallel to :attr:`ids`."""

    @property
    def skipped(self) -> List[str]:
        """Names of records dropped as empty or non-ACGT."""

    @property
    def method(self) -> str: ...
    @property
    def min_overlap(self) -> int: ...
    @property
    def edges(self) -> List[Edge]: ...
    @property
    def contigs(self) -> List[Contig]:
        """Assembled contigs, longest first. Computed on first access, then cached."""

    @property
    def stats(self) -> Stats: ...
    def edge_table(self) -> Dict[str, List[Any]]:
        """Edges as parallel columns, ready for ``pandas.DataFrame(g.edge_table())``."""

    def to_gfa(self, *, names: bool = True) -> str:
        """The graph in GFA 1.0. ``names=False`` reproduces the CLI's ``.gfa`` exactly."""

    def to_fasta(self, *, names: bool = True) -> str:
        """Per-sequence FASTA with edge annotations in the headers."""

    def contigs_fasta(self, *, names: bool = True) -> str:
        """Assembled contigs as FASTA."""

    def write_gfa(
        self, path: Union[str, "os.PathLike[str]"], *, names: bool = True
    ) -> None: ...
    def write_fasta(
        self, path: Union[str, "os.PathLike[str]"], *, names: bool = True
    ) -> None: ...
    def write_contigs_fasta(
        self, path: Union[str, "os.PathLike[str]"], *, names: bool = True
    ) -> None: ...
    def __len__(self) -> int: ...
    def __contains__(self, name: str) -> bool: ...
    def __getitem__(self, key: Union[str, int]) -> str: ...

def build_from_pairs(
    names: Sequence[str],
    sequences: Sequence[str],
    *,
    min_overlap: int = 20,
    method: str = "all",
    on_invalid: str = "skip",
) -> OverlapGraph: ...
def build_from_fasta(
    path: Union[str, "os.PathLike[str]"],
    *,
    min_overlap: int = 20,
    method: str = "all",
    on_invalid: str = "skip",
) -> OverlapGraph: ...
