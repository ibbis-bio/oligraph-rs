use std::fs::{self, File};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use clap::Parser;
use oligraph_rs::{
    AssemblyMethod, LIMBS, assemble_contigs, build_overlap_graph, parse_fasta_str,
    write_contigs_fasta, write_fasta, write_gfa,
};

#[derive(clap::ValueEnum, Clone, Copy)]
enum CliMethod {
    All,
    Pca,
}

impl From<CliMethod> for AssemblyMethod {
    fn from(m: CliMethod) -> Self {
        match m {
            CliMethod::All => AssemblyMethod::All,
            CliMethod::Pca => AssemblyMethod::Pca,
        }
    }
}

#[derive(Parser)]
#[command(
    name = "oligraph",
    version,
    about = "Overlap graph builder and contig assembler for oligonucleotide pools",
)]
struct Cli {
    /// Input FASTA file
    #[arg(short, long)]
    input: PathBuf,

    /// Output file prefix (writes .gfa, .fasta, .contigs.fasta)
    #[arg(short, long)]
    output: PathBuf,

    /// Minimum overlap length in bp [1-64]
    #[arg(short = 'l', long, default_value_t = 20, value_parser = clap::value_parser!(u32).range(1..=64))]
    min_overlap: u32,

    /// Assembly method
    #[arg(short, long, value_enum, default_value_t = CliMethod::All)]
    method: CliMethod,
}

fn parse_fasta(path: &Path) -> std::io::Result<Vec<Vec<u8>>> {
    let content = fs::read_to_string(path)?;
    let (seqs, skipped) = parse_fasta_str(&content);
    if skipped > 0 {
        eprintln!("skipped {} sequences with non-ACGT bases", skipped);
    }
    Ok(seqs)
}

fn main() {
    let cli = Cli::parse();

    let seqs = parse_fasta(&cli.input).unwrap_or_else(|e| {
        eprintln!("error reading {}: {}", cli.input.display(), e);
        std::process::exit(1);
    });

    if seqs.is_empty() {
        eprintln!("error: no sequences found in {}", cli.input.display());
        std::process::exit(1);
    }

    let max_len = seqs.iter().map(|s| s.len()).max().unwrap();
    if max_len > LIMBS * 32 {
        eprintln!(
            "error: sequence length {} exceeds capacity of LIMBS={} (max {}bp). Recompile with larger LIMBS.",
            max_len, LIMBS, LIMBS * 32
        );
        std::process::exit(1);
    }

    let l_min = cli.min_overlap;
    let assembly_method = AssemblyMethod::from(cli.method);

    eprintln!(
        "loaded {} sequences (lengths {}-{}) with l_min={}",
        seqs.len(),
        seqs.iter().map(|s| s.len()).min().unwrap(),
        max_len,
        l_min,
    );

    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let edges = build_overlap_graph::<LIMBS>(&seq_refs, l_min, assembly_method);
    eprintln!("found {} edges", edges.len());

    let prefix = &cli.output;
    if let Some(parent) = prefix.parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).unwrap_or_else(|e| {
                eprintln!("error creating directory {}: {}", parent.display(), e);
                std::process::exit(1);
            });
        }
    }

    let gfa_path = prefix.with_extension("gfa");
    let gfa_file = File::create(&gfa_path).unwrap_or_else(|e| {
        eprintln!("error creating {}: {}", gfa_path.display(), e);
        std::process::exit(1);
    });
    if let Err(e) = write_gfa(&seq_refs, &edges, BufWriter::new(gfa_file), assembly_method) {
        eprintln!("error writing GFA: {}", e);
        std::process::exit(1);
    }
    eprintln!("wrote {}", gfa_path.display());

    let fasta_out_path = prefix.with_extension("fasta");
    let fasta_file = File::create(&fasta_out_path).unwrap_or_else(|e| {
        eprintln!("error creating {}: {}", fasta_out_path.display(), e);
        std::process::exit(1);
    });
    if let Err(e) = write_fasta(&seq_refs, &edges, BufWriter::new(fasta_file)) {
        eprintln!("error writing FASTA: {}", e);
        std::process::exit(1);
    }
    eprintln!("wrote {}", fasta_out_path.display());

    let contigs = assemble_contigs(&seq_refs, &edges);
    if !contigs.is_empty() {
        let contigs_path = prefix.with_extension("contigs.fasta");
        let contigs_file = File::create(&contigs_path).unwrap_or_else(|e| {
            eprintln!("error creating {}: {}", contigs_path.display(), e);
            std::process::exit(1);
        });
        if let Err(e) = write_contigs_fasta(&contigs, BufWriter::new(contigs_file)) {
            eprintln!("error writing contigs: {}", e);
            std::process::exit(1);
        }
        eprintln!("wrote {} ({} contigs)", contigs_path.display(), contigs.len());
    } else {
        eprintln!("no contigs assembled (no connected components with edges)");
    }
}
