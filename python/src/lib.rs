//! Python bindings for OliGraph.
//!
//! The core returns typed errors (`oligraph_rs::Error`, `SequenceError`) rather
//! than panicking, so this layer translates rather than pre-validates: an
//! `Error` becomes a `ValueError` carrying the core's own message. The limits
//! themselves — the 320 bp capacity, the 1..=64 overlap range — are never
//! restated here; they come from `validate_sequence` and `validate_min_overlap`.
//!
//! What this layer does add is policy the core has no opinion on: `on_invalid`
//! decides whether an unusable record is dropped or fails the call, and record
//! names are kept so every edge, contig path and output file can be labelled
//! with them (the core identifies sequences by position, because
//! `parse_fasta_str` throws FASTA headers away).

use std::collections::HashMap;
use std::ffi::CString;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;

use oligraph_rs::{
    AssemblyMethod, Contig, Edge, LIMBS, Packed, SequenceError, Strand, Topology, assemble_contigs,
    build_overlap_graph, validate_min_overlap, validate_sequence,
};
use pyo3::exceptions::{PyIndexError, PyKeyError, PyTypeError, PyUserWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

mod writers;

use writers::{
    default_labels, write_contigs_fasta_labelled, write_fasta_labelled, write_gfa_labelled,
};

/// Longest sequence the core's fixed-width bit packing can hold: 10 limbs x 32
/// bases = 320 bp. Taken from the core so the two cannot disagree.
pub const MAX_SEQUENCE_LENGTH: usize = Packed::<LIMBS>::CAPACITY;

/// Largest permitted `min_overlap`: the seed key has to fit in a u128 at 2 bits
/// per base.
pub const MAX_OVERLAP: u32 = oligraph_rs::MAX_OVERLAP;

// ============================================================
// Errors
// ============================================================

/// Raised inside `Python::detach` blocks, where no `Python` token is available to
/// build a `PyErr` with.
enum BuildError {
    Value(String),
    Io(std::io::Error),
}

impl From<std::io::Error> for BuildError {
    fn from(e: std::io::Error) -> Self {
        BuildError::Io(e)
    }
}

impl From<BuildError> for PyErr {
    fn from(e: BuildError) -> PyErr {
        match e {
            BuildError::Value(msg) => PyValueError::new_err(msg),
            // PyO3 maps io::Error onto the matching OSError subclass
            // (FileNotFoundError, PermissionError, ...).
            BuildError::Io(err) => err.into(),
        }
    }
}

impl From<oligraph_rs::Error> for BuildError {
    fn from(e: oligraph_rs::Error) -> Self {
        BuildError::Value(e.to_string())
    }
}

// ============================================================
// Mirror enums
// ============================================================
//
// The upstream enums can't be used as pyclasses directly: `Strand` (lib.rs:111)
// has no `Ord` and `AssemblyMethod` (lib.rs:134) derives neither `Debug` nor
// `Hash`, so `#[pyclass(eq, eq_int)]` won't apply to them.

// `hash` is not optional here: setting `eq` makes Python null out `__hash__`, and
// strands are natural set members and dict keys.
#[pyclass(eq, eq_int, hash, frozen, from_py_object, module = "oligraph", name = "Strand")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PyStrand {
    Fwd,
    Rev,
}

#[pymethods]
impl PyStrand {
    /// `"+"` or `"-"`, matching the GFA and FASTA output.
    fn __str__(&self) -> &'static str {
        self.symbol()
    }

    /// `"+"` or `"-"`.
    #[getter]
    fn symbol(&self) -> &'static str {
        match self {
            PyStrand::Fwd => "+",
            PyStrand::Rev => "-",
        }
    }
}

impl From<Strand> for PyStrand {
    fn from(s: Strand) -> Self {
        match s {
            Strand::Fwd => PyStrand::Fwd,
            Strand::Rev => PyStrand::Rev,
        }
    }
}

#[pyclass(eq, eq_int, hash, frozen, from_py_object, module = "oligraph", name = "Topology")]
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum PyTopology {
    Linear,
    Cyclic,
}

#[pymethods]
impl PyTopology {
    fn __str__(&self) -> &'static str {
        match self {
            PyTopology::Linear => "linear",
            PyTopology::Cyclic => "cyclic",
        }
    }
}

impl From<Topology> for PyTopology {
    fn from(t: Topology) -> Self {
        match t {
            Topology::Linear => PyTopology::Linear,
            Topology::Cyclic => PyTopology::Cyclic,
        }
    }
}

fn parse_method(method: &str) -> PyResult<AssemblyMethod> {
    match method.to_ascii_lowercase().as_str() {
        "all" => Ok(AssemblyMethod::All),
        "pca" => Ok(AssemblyMethod::Pca),
        other => Err(PyValueError::new_err(format!(
            "invalid method {other:?}: expected \"all\" or \"pca\""
        ))),
    }
}

fn method_name(method: AssemblyMethod) -> &'static str {
    match method {
        AssemblyMethod::All => "all",
        AssemblyMethod::Pca => "pca",
    }
}

#[derive(Clone, Copy)]
enum OnInvalid {
    Skip,
    Raise,
}

fn parse_on_invalid(value: &str) -> PyResult<OnInvalid> {
    match value.to_ascii_lowercase().as_str() {
        "skip" => Ok(OnInvalid::Skip),
        "raise" => Ok(OnInvalid::Raise),
        other => Err(PyValueError::new_err(format!(
            "invalid on_invalid {other:?}: expected \"skip\" or \"raise\""
        ))),
    }
}

// ============================================================
// Result types
// ============================================================

/// One overlap between two oligos.
#[pyclass(frozen, skip_from_py_object, module = "oligraph", name = "Edge")]
#[derive(Clone)]
pub struct PyEdge {
    /// Position of the source sequence in `OverlapGraph.ids`.
    #[pyo3(get)]
    pub from_index: usize,
    /// Name of the source sequence.
    #[pyo3(get)]
    pub from_name: String,
    #[pyo3(get)]
    pub from_strand: PyStrand,
    #[pyo3(get)]
    pub to_index: usize,
    #[pyo3(get)]
    pub to_name: String,
    #[pyo3(get)]
    pub to_strand: PyStrand,
    /// Length of the overlap in bases.
    #[pyo3(get)]
    pub overlap_len: u32,
}

#[pymethods]
impl PyEdge {
    fn __repr__(&self) -> String {
        format!(
            "Edge({}{} -> {}{}, overlap_len={})",
            self.from_name,
            self.from_strand.symbol(),
            self.to_name,
            self.to_strand.symbol(),
            self.overlap_len
        )
    }
}

impl PyEdge {
    fn from_core(e: &Edge, names: &[String]) -> Self {
        PyEdge {
            from_index: e.from_id as usize,
            from_name: names[e.from_id as usize].clone(),
            from_strand: e.from_strand.into(),
            to_index: e.to_id as usize,
            to_name: names[e.to_id as usize].clone(),
            to_strand: e.to_strand.into(),
            overlap_len: e.overlap_len,
        }
    }
}

/// A contig walked out of one connected component of the overlap graph.
#[pyclass(frozen, skip_from_py_object, module = "oligraph", name = "Contig")]
#[derive(Clone)]
pub struct PyContig {
    /// The stitched sequence.
    #[pyo3(get)]
    pub sequence: String,
    /// Connected-component id; 0 is the largest component.
    #[pyo3(get)]
    pub component: usize,
    /// Branch points where the greedy walk had to choose between neighbours.
    #[pyo3(get)]
    pub branches: u32,
    topology: PyTopology,
    path: Vec<(String, PyStrand)>,
    path_indices: Vec<(usize, PyStrand)>,
}

#[pymethods]
impl PyContig {
    /// `"linear"` or `"cyclic"`.
    #[getter]
    fn topology(&self) -> PyTopology {
        self.topology
    }

    /// The oligos walked, in order, as `(name, strand)`.
    #[getter]
    fn path(&self) -> Vec<(String, PyStrand)> {
        self.path.clone()
    }

    /// The same walk as `(index, strand)`, indexing into `OverlapGraph.ids`.
    #[getter]
    fn path_indices(&self) -> Vec<(usize, PyStrand)> {
        self.path_indices.clone()
    }

    /// Number of oligos in the walk.
    #[getter]
    fn n_oligos(&self) -> usize {
        self.path.len()
    }

    /// Length of the assembled sequence in bases.
    fn __len__(&self) -> usize {
        self.sequence.len()
    }

    fn __repr__(&self) -> String {
        format!(
            "Contig(component={}, oligos={}, length={}, topology={}, branches={})",
            self.component,
            self.path.len(),
            self.sequence.len(),
            self.topology.__str__(),
            self.branches
        )
    }
}

impl PyContig {
    fn from_core(c: &Contig, names: &[String]) -> Self {
        PyContig {
            // Bases were validated as ACGT on the way in, so this is always UTF-8.
            sequence: String::from_utf8_lossy(&c.sequence).into_owned(),
            component: c.component,
            branches: c.branches,
            topology: c.topology.into(),
            path: c
                .path
                .iter()
                .map(|&(id, s)| (names[id as usize].clone(), s.into()))
                .collect(),
            path_indices: c
                .path
                .iter()
                .map(|&(id, s)| (id as usize, s.into()))
                .collect(),
        }
    }
}

/// Summary counts for a graph.
#[pyclass(frozen, skip_from_py_object, module = "oligraph", name = "Stats")]
#[derive(Clone)]
pub struct PyStats {
    /// Sequences that made it into the graph.
    #[pyo3(get)]
    pub n_sequences: usize,
    /// Records dropped for being empty or containing non-ACGT bases.
    #[pyo3(get)]
    pub n_skipped: usize,
    #[pyo3(get)]
    pub len_min: usize,
    #[pyo3(get)]
    pub len_max: usize,
    #[pyo3(get)]
    pub n_edges: usize,
    /// Sequences with no overlap to any other sequence.
    #[pyo3(get)]
    pub n_isolated: usize,
    #[pyo3(get)]
    pub n_contigs: usize,
}

#[pymethods]
impl PyStats {
    fn to_dict<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let d = PyDict::new(py);
        d.set_item("n_sequences", self.n_sequences)?;
        d.set_item("n_skipped", self.n_skipped)?;
        d.set_item("len_min", self.len_min)?;
        d.set_item("len_max", self.len_max)?;
        d.set_item("n_edges", self.n_edges)?;
        d.set_item("n_isolated", self.n_isolated)?;
        d.set_item("n_contigs", self.n_contigs)?;
        Ok(d)
    }

    fn __repr__(&self) -> String {
        format!(
            "Stats(n_sequences={}, n_skipped={}, len_min={}, len_max={}, \
             n_edges={}, n_isolated={}, n_contigs={})",
            self.n_sequences,
            self.n_skipped,
            self.len_min,
            self.len_max,
            self.n_edges,
            self.n_isolated,
            self.n_contigs
        )
    }
}

// ============================================================
// OverlapGraph
// ============================================================

/// The overlap graph for a set of oligos, with names preserved.
#[pyclass(frozen, module = "oligraph", name = "OverlapGraph")]
pub struct PyOverlapGraph {
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
    edges: Vec<Edge>,
    skipped: Vec<String>,
    method: AssemblyMethod,
    min_overlap: u32,
    index: HashMap<String, usize>,
    contigs: OnceLock<Vec<Contig>>,
}

impl PyOverlapGraph {
    fn seq_refs(&self) -> Vec<&[u8]> {
        self.seqs.iter().map(|s| s.as_slice()).collect()
    }

    /// Runs `assemble_contigs` on first call, then reuses the result.
    ///
    /// `OnceLock::get_or_init` cannot carry an error out, so compute first and
    /// store second; a benign race just recomputes and discards.
    fn core_contigs(&self, py: Python<'_>) -> PyResult<&[Contig]> {
        if let Some(contigs) = self.contigs.get() {
            return Ok(contigs);
        }
        let contigs = py
            .detach(|| assemble_contigs(&self.seq_refs(), &self.edges))
            .map_err(BuildError::from)?;
        Ok(self.contigs.get_or_init(|| contigs))
    }

    fn labels(&self, names: bool) -> Vec<String> {
        if names {
            self.names.clone()
        } else {
            default_labels(self.names.len())
        }
    }

    /// Render one of the writers into a `String`. Kept on the calling thread:
    /// formatting is trivial next to `build_overlap_graph`, and the file-writing
    /// variants below release the GIL for the actual I/O.
    fn render<F>(&self, f: F) -> PyResult<String>
    where
        F: FnOnce(&mut Vec<u8>) -> std::io::Result<()>,
    {
        let mut buf = Vec::new();
        f(&mut buf)?;
        String::from_utf8(buf).map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

#[pymethods]
impl PyOverlapGraph {
    /// Sequence names, in the order the graph indexes them.
    #[getter]
    fn ids(&self) -> Vec<String> {
        self.names.clone()
    }

    /// Sequences, uppercased, parallel to `ids`.
    #[getter]
    fn sequences(&self) -> Vec<String> {
        self.seqs
            .iter()
            .map(|s| String::from_utf8_lossy(s).into_owned())
            .collect()
    }

    /// Names of records dropped for being empty or containing non-ACGT bases.
    #[getter]
    fn skipped(&self) -> Vec<String> {
        self.skipped.clone()
    }

    /// `"all"` or `"pca"`.
    #[getter]
    fn method(&self) -> &'static str {
        method_name(self.method)
    }

    /// The `min_overlap` this graph was built with.
    #[getter]
    fn min_overlap(&self) -> u32 {
        self.min_overlap
    }

    #[getter]
    fn edges(&self) -> Vec<PyEdge> {
        self.edges
            .iter()
            .map(|e| PyEdge::from_core(e, &self.names))
            .collect()
    }

    /// Assembled contigs, longest first. Computed on first access, then cached.
    #[getter]
    fn contigs(&self, py: Python<'_>) -> PyResult<Vec<PyContig>> {
        Ok(self
            .core_contigs(py)?
            .iter()
            .map(|c| PyContig::from_core(c, &self.names))
            .collect())
    }

    #[getter]
    fn stats(&self, py: Python<'_>) -> PyResult<PyStats> {
        // Mirrors web/src/main.rs:78-85.
        let mut has_edge = vec![false; self.seqs.len()];
        for e in &self.edges {
            if e.from_id != e.to_id {
                has_edge[e.from_id as usize] = true;
                has_edge[e.to_id as usize] = true;
            }
        }
        Ok(PyStats {
            n_sequences: self.seqs.len(),
            n_skipped: self.skipped.len(),
            len_min: self.seqs.iter().map(|s| s.len()).min().unwrap_or(0),
            len_max: self.seqs.iter().map(|s| s.len()).max().unwrap_or(0),
            n_edges: self.edges.len(),
            n_isolated: has_edge.iter().filter(|&&v| !v).count(),
            n_contigs: self.core_contigs(py)?.len(),
        })
    }

    /// Edges as parallel columns, ready for `pandas.DataFrame(g.edge_table())`.
    fn edge_table<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        let n = self.edges.len();
        let mut from_index = Vec::with_capacity(n);
        let mut from_name = Vec::with_capacity(n);
        let mut from_strand = Vec::with_capacity(n);
        let mut to_index = Vec::with_capacity(n);
        let mut to_name = Vec::with_capacity(n);
        let mut to_strand = Vec::with_capacity(n);
        let mut overlap_len = Vec::with_capacity(n);

        for e in &self.edges {
            from_index.push(e.from_id as usize);
            from_name.push(self.names[e.from_id as usize].clone());
            from_strand.push(PyStrand::from(e.from_strand).symbol());
            to_index.push(e.to_id as usize);
            to_name.push(self.names[e.to_id as usize].clone());
            to_strand.push(PyStrand::from(e.to_strand).symbol());
            overlap_len.push(e.overlap_len);
        }

        let d = PyDict::new(py);
        d.set_item("from_index", from_index)?;
        d.set_item("from_name", from_name)?;
        d.set_item("from_strand", from_strand)?;
        d.set_item("to_index", to_index)?;
        d.set_item("to_name", to_name)?;
        d.set_item("to_strand", to_strand)?;
        d.set_item("overlap_len", overlap_len)?;
        Ok(d)
    }

    /// The overlap graph in GFA 1.0. With `names=False`, segment names are
    /// positional indices, matching the CLI's `.gfa` output exactly.
    #[pyo3(signature = (*, names = true))]
    fn to_gfa(&self, names: bool) -> PyResult<String> {
        let labels = self.labels(names);
        let refs = self.seq_refs();
        let method = self.method;
        let edges = &self.edges;
        self.render(|buf| write_gfa_labelled(&refs, &labels, edges, buf, method))
    }

    /// Per-sequence FASTA with edge annotations in the headers, matching the CLI's
    /// `.fasta` output when `names=False`.
    #[pyo3(signature = (*, names = true))]
    fn to_fasta(&self, names: bool) -> PyResult<String> {
        let labels = self.labels(names);
        let refs = self.seq_refs();
        let edges = &self.edges;
        self.render(|buf| write_fasta_labelled(&refs, &labels, edges, buf))
    }

    /// Assembled contigs as FASTA, matching the CLI's `.contigs.fasta` output when
    /// `names=False`.
    #[pyo3(signature = (*, names = true))]
    fn contigs_fasta(&self, py: Python<'_>, names: bool) -> PyResult<String> {
        let labels = self.labels(names);
        let contigs = self.core_contigs(py)?;
        self.render(|buf| write_contigs_fasta_labelled(contigs, &labels, buf))
    }

    /// Write the GFA straight to `path` without building the whole string first.
    #[pyo3(signature = (path, *, names = true))]
    fn write_gfa(&self, py: Python<'_>, path: PathBuf, names: bool) -> PyResult<()> {
        let labels = self.labels(names);
        let refs = self.seq_refs();
        let method = self.method;
        py.detach(|| -> std::io::Result<()> {
            let f = BufWriter::new(File::create(&path)?);
            write_gfa_labelled(&refs, &labels, &self.edges, f, method)
        })?;
        Ok(())
    }

    #[pyo3(signature = (path, *, names = true))]
    fn write_fasta(&self, py: Python<'_>, path: PathBuf, names: bool) -> PyResult<()> {
        let labels = self.labels(names);
        let refs = self.seq_refs();
        py.detach(|| -> std::io::Result<()> {
            let f = BufWriter::new(File::create(&path)?);
            write_fasta_labelled(&refs, &labels, &self.edges, f)
        })?;
        Ok(())
    }

    #[pyo3(signature = (path, *, names = true))]
    fn write_contigs_fasta(&self, py: Python<'_>, path: PathBuf, names: bool) -> PyResult<()> {
        let labels = self.labels(names);
        let contigs = self.core_contigs(py)?;
        py.detach(|| -> std::io::Result<()> {
            let f = BufWriter::new(File::create(&path)?);
            write_contigs_fasta_labelled(contigs, &labels, f)
        })?;
        Ok(())
    }

    /// Number of sequences in the graph.
    fn __len__(&self) -> usize {
        self.seqs.len()
    }

    fn __contains__(&self, name: &str) -> bool {
        self.index.contains_key(name)
    }

    /// `g["oligo_1"]` or `g[0]` -> the sequence.
    fn __getitem__(&self, key: &Bound<'_, PyAny>) -> PyResult<String> {
        if let Ok(name) = key.extract::<String>() {
            return match self.index.get(&name) {
                Some(&i) => Ok(String::from_utf8_lossy(&self.seqs[i]).into_owned()),
                None => Err(PyKeyError::new_err(name)),
            };
        }
        if let Ok(i) = key.extract::<usize>() {
            return match self.seqs.get(i) {
                Some(s) => Ok(String::from_utf8_lossy(s).into_owned()),
                None => Err(PyIndexError::new_err(format!(
                    "index {i} out of range for {} sequences",
                    self.seqs.len()
                ))),
            };
        }
        Err(PyTypeError::new_err(
            "OverlapGraph indices must be a sequence name (str) or a position (int)",
        ))
    }

    fn __repr__(&self) -> String {
        format!(
            "<OverlapGraph {} sequences, {} edges, min_overlap={}, method='{}'>",
            self.seqs.len(),
            self.edges.len(),
            self.min_overlap,
            method_name(self.method)
        )
    }
}

// ============================================================
// Validation
// ============================================================

struct Prepared {
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
    skipped: Vec<String>,
}

/// Uppercase and apply the `on_invalid` policy.
///
/// The rules themselves come from `oligraph_rs::validate_sequence`, so the limits
/// live in one place. What is added here is a *policy* the core has no opinion on:
/// which records to drop and carry on with, versus which to fail the whole call.
fn prepare(
    names: Vec<String>,
    seqs: Vec<Vec<u8>>,
    on_invalid: OnInvalid,
) -> Result<Prepared, BuildError> {
    if names.len() != seqs.len() {
        return Err(BuildError::Value(format!(
            "got {} names for {} sequences",
            names.len(),
            seqs.len()
        )));
    }
    if names.is_empty() {
        return Err(BuildError::Value("no sequences given".to_string()));
    }

    let total = names.len();
    let mut kept_names: Vec<String> = Vec::with_capacity(total);
    let mut kept_seqs: Vec<Vec<u8>> = Vec::with_capacity(total);
    let mut skipped: Vec<String> = Vec::new();

    for (name, seq) in names.into_iter().zip(seqs) {
        let seq: Vec<u8> = seq.iter().map(|b| b.to_ascii_uppercase()).collect();

        let problem = match validate_sequence::<LIMBS>(&seq) {
            // A capacity limit rather than dirty data, so this always fails the
            // call — `on_invalid` does not apply. Adding the record's name and a
            // remedy is the only thing done to the core's message here.
            Err(e @ SequenceError::TooLong { .. }) => {
                return Err(BuildError::Value(format!(
                    "record {name:?}: {e}. Raising that limit means recompiling \
                     the Rust core with a larger LIMBS"
                )));
            }
            Err(e) => Some(e.to_string()),
            // The core tolerates empty sequences (they simply form no edges), but
            // silently keeping a blank record is not useful to a caller.
            Ok(()) if seq.is_empty() => Some("sequence is empty".to_string()),
            Ok(()) => None,
        };

        match problem {
            None => {
                kept_names.push(name);
                kept_seqs.push(seq);
            }
            Some(why) => match on_invalid {
                OnInvalid::Skip => skipped.push(name),
                OnInvalid::Raise => {
                    return Err(BuildError::Value(format!("record {name:?}: {why}")));
                }
            },
        }
    }

    if kept_seqs.is_empty() {
        return Err(BuildError::Value(format!(
            "no usable sequences: all {total} records were empty or contained \
             non-ACGT bases"
        )));
    }

    Ok(Prepared {
        names: kept_names,
        seqs: kept_seqs,
        skipped,
    })
}

/// Reject a bad `min_overlap` before reading a file or copying sequences.
/// `build_overlap_graph` checks it too; this only moves the failure earlier.
fn check_min_overlap(min_overlap: u32) -> PyResult<()> {
    validate_min_overlap(min_overlap).map_err(|e| PyErr::from(BuildError::from(e)))
}

/// Sequences shorter than `min_overlap` are skipped by the core (lib.rs:255) and
/// simply produce no edges. Silently returning an empty graph is confusing, so
/// warn instead — the web app treats this as a hard error, which is too strict for
/// a pool where one oligo happens to be short.
fn warn_short_sequences(py: Python<'_>, seqs: &[Vec<u8>], min_overlap: u32) -> PyResult<()> {
    let too_short = seqs
        .iter()
        .filter(|s| (s.len() as u32) < min_overlap)
        .count();
    if too_short == 0 {
        return Ok(());
    }
    let len_min = seqs.iter().map(|s| s.len()).min().unwrap_or(0);
    let msg = CString::new(format!(
        "min_overlap ({min_overlap}) is longer than {too_short} of {} sequences \
         (shortest is {len_min} bp); those sequences cannot form any edge",
        seqs.len()
    ))
    .map_err(|e| PyValueError::new_err(e.to_string()))?;
    let category = py.get_type::<PyUserWarning>();
    // stacklevel 2 skips the wrapper in oligraph/__init__.py and points at the
    // caller's own code.
    PyErr::warn(py, &category, msg.as_c_str(), 2)
}

/// Validate, then run the core. Follows main.rs:62-95.
fn build_graph(
    py: Python<'_>,
    prepared: Prepared,
    min_overlap: u32,
    method: AssemblyMethod,
) -> PyResult<PyOverlapGraph> {
    warn_short_sequences(py, &prepared.seqs, min_overlap)?;

    let Prepared {
        names,
        seqs,
        skipped,
    } = prepared;

    // Pure CPU work touching no Python objects, so let other threads run.
    let edges = py
        .detach(|| {
            let refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
            build_overlap_graph::<LIMBS>(&refs, min_overlap, method)
        })
        .map_err(BuildError::from)?;

    // First occurrence wins, so duplicate FASTA record names stay usable.
    let mut index = HashMap::with_capacity(names.len());
    for (i, name) in names.iter().enumerate() {
        index.entry(name.clone()).or_insert(i);
    }

    Ok(PyOverlapGraph {
        names,
        seqs,
        edges,
        skipped,
        method,
        min_overlap,
        index,
        contigs: OnceLock::new(),
    })
}

// ============================================================
// Entry points
// ============================================================

/// Build a graph from parallel name and sequence lists.
///
/// The Python wrapper in `oligraph/__init__.py` normalises dicts, pair iterables
/// and bare sequence iterables down to this.
#[pyfunction]
#[pyo3(signature = (names, sequences, *, min_overlap = 20, method = "all", on_invalid = "skip"))]
fn build_from_pairs(
    py: Python<'_>,
    names: Vec<String>,
    sequences: Vec<String>,
    min_overlap: u32,
    method: &str,
    on_invalid: &str,
) -> PyResult<PyOverlapGraph> {
    check_min_overlap(min_overlap)?;
    let method = parse_method(method)?;
    let on_invalid = parse_on_invalid(on_invalid)?;

    let seqs: Vec<Vec<u8>> = sequences.into_iter().map(String::into_bytes).collect();
    let prepared = prepare(names, seqs, on_invalid)?;
    build_graph(py, prepared, min_overlap, method)
}

/// Read records with needletail, then build a graph.
///
/// Handles FASTA and FASTQ, single- or multi-line, plain or gzipped.
#[pyfunction]
#[pyo3(signature = (path, *, min_overlap = 20, method = "all", on_invalid = "skip"))]
fn build_from_fasta(
    py: Python<'_>,
    path: PathBuf,
    min_overlap: u32,
    method: &str,
    on_invalid: &str,
) -> PyResult<PyOverlapGraph> {
    check_min_overlap(min_overlap)?;
    let method = parse_method(method)?;
    let on_invalid = parse_on_invalid(on_invalid)?;

    let prepared = py.detach(|| -> Result<Prepared, BuildError> {
        let (names, seqs) = read_records(&path)?;
        prepare(names, seqs, on_invalid)
    })?;

    build_graph(py, prepared, min_overlap, method)
}

fn read_records(path: &Path) -> Result<(Vec<String>, Vec<Vec<u8>>), BuildError> {
    // Probe first so a missing or unreadable file yields a real io::Error, which
    // PyO3 turns into FileNotFoundError/PermissionError rather than a generic
    // parse failure.
    let probe = File::open(path)?;

    // needletail reports an empty file as "Failed to read the first two bytes",
    // which is not a useful thing to show a caller.
    if probe.metadata()?.len() == 0 {
        return Err(BuildError::Value(format!(
            "no sequence records found in {}: the file is empty",
            path.display()
        )));
    }

    let mut reader = needletail::parse_fastx_file(path).map_err(|e| {
        BuildError::Value(format!("could not parse {}: {e}", path.display()))
    })?;

    let mut names = Vec::new();
    let mut seqs = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.map_err(|e| {
            BuildError::Value(format!("malformed record in {}: {e}", path.display()))
        })?;
        // `id()` is the whole header line; the record name is everything up to the
        // first whitespace, the rest is a free-text description.
        let id = record.id();
        let name = id
            .split(|&b| b == b' ' || b == b'\t')
            .next()
            .unwrap_or(id);
        names.push(String::from_utf8_lossy(name).into_owned());
        // `seq()` already strips line breaks; casing and validation happen in
        // `prepare` so file and in-memory input behave identically.
        seqs.push(record.seq().into_owned());
    }

    if names.is_empty() {
        return Err(BuildError::Value(format!(
            "no sequence records found in {}",
            path.display()
        )));
    }
    Ok((names, seqs))
}

#[pymodule]
fn _oligraph(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("MAX_SEQUENCE_LENGTH", MAX_SEQUENCE_LENGTH)?;
    m.add("MAX_OVERLAP", MAX_OVERLAP)?;
    m.add_class::<PyStrand>()?;
    m.add_class::<PyTopology>()?;
    m.add_class::<PyEdge>()?;
    m.add_class::<PyContig>()?;
    m.add_class::<PyStats>()?;
    m.add_class::<PyOverlapGraph>()?;
    m.add_function(wrap_pyfunction!(build_from_pairs, m)?)?;
    m.add_function(wrap_pyfunction!(build_from_fasta, m)?)?;
    Ok(())
}
