use std::env;
use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

use oligraph_rs::{
    AssemblyMethod, LIMBS, assemble_contigs, build_overlap_graph, parse_fasta_str, write_contigs_fasta,
    write_fasta, write_gfa,
};

fn parse_fasta(path: &str) -> std::io::Result<Vec<Vec<u8>>> {
    let content = fs::read_to_string(path)?;
    let (seqs, skipped) = parse_fasta_str(&content);
    if skipped > 0 {
        eprintln!("skipped {} sequences with non-ACGT bases", skipped);
    }
    Ok(seqs)
}

fn usage() -> ! {
    eprintln!("Usage: oligraph-rs <input.fasta> [output.gfa] [-l <min_overlap>] [--assembly-method <all|pca>]");
    eprintln!("  input.fasta                Input FASTA file of sequences");
    eprintln!("  output.gfa                 Output GFA file (default: stdout, GFA only)");
    eprintln!("                             Also writes .fasta alongside the .gfa");
    eprintln!("  -l <min>                   Minimum overlap length (default: 20, max: 32)");
    eprintln!("  --assembly-method <method>  Assembly method filter (default: all)");
    std::process::exit(1);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 2 {
        usage();
    }

    let mut fasta_path: Option<&str> = None;
    let mut gfa_path: Option<&str> = None;
    let mut l_min: u32 = 20;
    let mut assembly_method = AssemblyMethod::All;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-l" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("error: -l requires a value");
                    usage();
                }
                l_min = args[i].parse().unwrap_or_else(|_| {
                    eprintln!("error: invalid value for -l: {}", args[i]);
                    usage();
                });
            }
            "--assembly-method" => {
                i += 1;
                if i >= args.len() {
                    eprintln!("error: --assembly-method requires a value");
                    usage();
                }
                assembly_method = match args[i].as_str() {
                    "all" => AssemblyMethod::All,
                    "pca" => AssemblyMethod::Pca,
                    _ => {
                        eprintln!("error: unknown assembly method '{}' (expected: all, pca)", args[i]);
                        usage();
                    }
                };
            }
            "-" | "--" => {
                eprintln!("error: unknown flag: {}", args[i]);
                usage();
            }
            arg if arg.starts_with('-') => {
                eprintln!("error: unknown flag: {}", arg);
                usage();
            }
            _ => {
                if fasta_path.is_none() {
                    fasta_path = Some(&args[i]);
                } else if gfa_path.is_none() {
                    gfa_path = Some(&args[i]);
                } else {
                    eprintln!("error: unexpected argument: {}", args[i]);
                    usage();
                }
            }
        }
        i += 1;
    }

    let fasta_path = fasta_path.unwrap_or_else(|| {
        eprintln!("error: no input FASTA file specified");
        usage();
    });

    let seqs = parse_fasta(fasta_path).unwrap_or_else(|e| {
        eprintln!("error reading {}: {}", fasta_path, e);
        std::process::exit(1);
    });

    if seqs.is_empty() {
        eprintln!("error: no sequences found in {}", fasta_path);
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

    eprintln!(
        "loaded {} sequences (lengths {}-{}) with l_min={}, assembly_method={}",
        seqs.len(),
        seqs.iter().map(|s| s.len()).min().unwrap(),
        max_len,
        l_min,
        match assembly_method {
            AssemblyMethod::All => "all",
            AssemblyMethod::Pca => "pca",
        }
    );

    let seq_refs: Vec<&[u8]> = seqs.iter().map(|s| s.as_slice()).collect();
    let edges = build_overlap_graph::<LIMBS>(&seq_refs, l_min, assembly_method);
    eprintln!("found {} edges", edges.len());

    match gfa_path {
        Some(path) => {
            let gfa_file = File::create(path).unwrap_or_else(|e| {
                eprintln!("error creating {}: {}", path, e);
                std::process::exit(1);
            });
            if let Err(e) = write_gfa(&seq_refs, &edges, BufWriter::new(gfa_file), assembly_method)
            {
                eprintln!("error writing GFA: {}", e);
                std::process::exit(1);
            }

            let fasta_path = Path::new(path).with_extension("fasta");
            let fasta_file = File::create(&fasta_path).unwrap_or_else(|e| {
                eprintln!("error creating {}: {}", fasta_path.display(), e);
                std::process::exit(1);
            });
            if let Err(e) = write_fasta(&seq_refs, &edges, BufWriter::new(fasta_file)) {
                eprintln!("error writing FASTA: {}", e);
                std::process::exit(1);
            }
            eprintln!("wrote {}", fasta_path.display());

            let contigs = assemble_contigs(&seq_refs, &edges);
            if !contigs.is_empty() {
                let contigs_path = Path::new(path).with_extension("contigs.fasta");
                let contigs_file = File::create(&contigs_path).unwrap_or_else(|e| {
                    eprintln!("error creating {}: {}", contigs_path.display(), e);
                    std::process::exit(1);
                });
                if let Err(e) = write_contigs_fasta(&contigs, BufWriter::new(contigs_file)) {
                    eprintln!("error writing contigs: {}", e);
                    std::process::exit(1);
                }
                eprintln!(
                    "wrote {} ({} contigs)",
                    contigs_path.display(),
                    contigs.len()
                );
            } else {
                eprintln!("no contigs assembled (no connected components with edges)");
            }
        }
        None => {
            if let Err(e) =
                write_gfa(&seq_refs, &edges, BufWriter::new(std::io::stdout()), assembly_method)
            {
                eprintln!("error writing GFA: {}", e);
                std::process::exit(1);
            }
        }
    }
}
