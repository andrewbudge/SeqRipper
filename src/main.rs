mod genbank;
mod genes;

use std::collections::HashMap;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use clap::Parser;
use rayon::prelude::*;

const VERSION: &str = "2.0.0";

/// SeqRipper - Fast, concurrent tool for extracting genes from GenBank files
#[derive(Parser)]
#[command(name = "seqripper", version = VERSION, about)]
struct Cli {
    /// Comma-separated genes to extract (e.g., COI,12S,CYTB)
    #[arg(short = 'g', value_delimiter = ',')]
    genes: Vec<String>,

    /// Extract ALL annotated genes
    #[arg(long = "all")]
    extract_all: bool,

    /// Skip tRNA genes (only with --all)
    #[arg(long = "no-trna")]
    skip_trna: bool,

    /// Show detailed per-file output
    #[arg(short = 'v', long = "verbose")]
    verbose: bool,

    /// Suppress all output
    #[arg(short = 'q', long = "quiet")]
    quiet: bool,

    /// Input GenBank files followed by output directory
    #[arg(required = true, num_args = 2..)]
    paths: Vec<PathBuf>,
}

#[derive(Debug)]
struct ExtractionResult {
    accession: String,
    gene: String,
    length_bp: usize,
    status: String,
    filename: String,
}

/// Thread-safe unique filename tracker.
struct FilenameTracker {
    used: Mutex<HashMap<String, usize>>,
}

impl FilenameTracker {
    fn new() -> Self {
        Self {
            used: Mutex::new(HashMap::new()),
        }
    }

    fn get_unique(&self, base: &str) -> (String, bool) {
        let mut map = self.used.lock().unwrap();
        let count = map.entry(base.to_string()).or_insert(0);
        *count += 1;
        if *count == 1 {
            (base.to_string(), false)
        } else {
            let ext = Path::new(base)
                .extension()
                .map(|e| format!(".{}", e.to_string_lossy()))
                .unwrap_or_default();
            let stem = base.strip_suffix(&ext).unwrap_or(base);
            (format!("{}_{}{}", stem, *count - 1, ext), true)
        }
    }
}

/// Write a FASTA file with 60-character line wrapping.
fn write_fasta(path: &Path, header: &str, sequence: &str) -> std::io::Result<()> {
    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, ">{}", header)?;
    for chunk in sequence.as_bytes().chunks(60) {
        w.write_all(chunk)?;
        w.write_all(b"\n")?;
    }
    w.flush()
}

/// Write TSV log file.
fn write_log(path: &Path, results: &mut [ExtractionResult]) -> std::io::Result<()> {
    results.sort_by(|a, b| {
        a.accession
            .cmp(&b.accession)
            .then(a.gene.cmp(&b.gene))
    });

    let file = fs::File::create(path)?;
    let mut w = BufWriter::new(file);
    writeln!(w, "Accession\tGene\tLength_bp\tStatus\tFilename")?;
    for r in results.iter() {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}",
            r.accession, r.gene, r.length_bp, r.status, r.filename
        )?;
    }
    w.flush()
}

/// Process a single GenBank file and return extraction results.
fn process_file(
    path: &Path,
    extract_all: bool,
    skip_trna: bool,
    target_genes: &HashMap<String, Vec<String>>,
    output_dir: &Path,
    tracker: &FilenameTracker,
    verbose: bool,
) -> Vec<ExtractionResult> {
    let mut results = Vec::new();

    let records = match genbank::parse_genbank(path) {
        Ok(r) => r,
        Err(e) => {
            if verbose {
                eprintln!("  Error parsing {}: {}", path.display(), e);
            }
            results.push(ExtractionResult {
                accession: path
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .into_owned(),
                gene: "-".to_string(),
                length_bp: 0,
                status: "error".to_string(),
                filename: "-".to_string(),
            });
            return results;
        }
    };

    for record in &records {
        if verbose {
            eprintln!(
                "  Processing: {} ({})",
                record.accession, record.organism
            );
        }

        let mut found_genes: HashMap<String, bool> = HashMap::new();

        for feature in &record.features {
            let gene_name = match genes::get_gene_name(feature) {
                Some(n) => n,
                None => continue,
            };

            if skip_trna && genes::is_trna(feature, gene_name) {
                continue;
            }

            let target_name = if extract_all {
                genes::clean_name(gene_name)
            } else {
                match genes::match_gene(gene_name, target_genes) {
                    Some(name) => name,
                    None => continue,
                }
            };

            if found_genes.contains_key(&target_name) {
                continue;
            }

            let sequence = genes::extract_sequence(record, feature);
            if sequence.is_empty() {
                continue;
            }

            let base_filename = format!("{}_{}.fasta", record.accession, target_name);
            let (unique_filename, is_dup) = tracker.get_unique(&base_filename);

            if is_dup && verbose {
                eprintln!(
                    "    Duplicate: {} -> {}",
                    base_filename, unique_filename
                );
            }

            let output_path = output_dir.join(&unique_filename);
            let loc_str = genes::get_location_string(feature);
            let header = format!(
                "{}:{} {} [gene={}]",
                record.accession, loc_str, record.organism, target_name
            );

            if let Err(e) = write_fasta(&output_path, &header, &sequence) {
                if verbose {
                    eprintln!("    Failed to write {}: {}", unique_filename, e);
                }
                results.push(ExtractionResult {
                    accession: record.accession.clone(),
                    gene: target_name,
                    length_bp: 0,
                    status: "error".to_string(),
                    filename: "-".to_string(),
                });
                continue;
            }

            found_genes.insert(target_name.clone(), true);
            let seq_len = sequence.len();

            results.push(ExtractionResult {
                accession: record.accession.clone(),
                gene: target_name.clone(),
                length_bp: seq_len,
                status: "extracted".to_string(),
                filename: unique_filename.clone(),
            });

            if verbose {
                eprintln!(
                    "    {} -> {} ({} bp)",
                    target_name, unique_filename, seq_len
                );
            }
        }

        // Report missing genes in specific gene mode
        if !extract_all {
            for gene in target_genes.keys() {
                if !found_genes.contains_key(gene) {
                    results.push(ExtractionResult {
                        accession: record.accession.clone(),
                        gene: gene.clone(),
                        length_bp: 0,
                        status: "missing".to_string(),
                        filename: "-".to_string(),
                    });
                    if verbose {
                        eprintln!("    {}: not found", gene);
                    }
                }
            }
        }
    }

    results
}

fn main() {
    let cli = Cli::parse();

    // Validate gene selection
    if !cli.extract_all && cli.genes.is_empty() {
        eprintln!("Error: either -g (genes) or --all is required");
        std::process::exit(1);
    }

    if cli.skip_trna && !cli.extract_all {
        eprintln!("Warning: --no-trna only works with --all mode");
    }

    // Split paths: last is output dir, rest are input files
    let (input_paths, output_dir) = match cli.paths.split_last() {
        Some((last, rest)) if !rest.is_empty() => (rest, last.clone()),
        _ => {
            eprintln!("Error: need at least one input file and an output directory");
            std::process::exit(1);
        }
    };

    // Verify input files
    let valid_files: Vec<&PathBuf> = input_paths
        .iter()
        .filter(|p| {
            if p.exists() {
                true
            } else {
                if !cli.quiet {
                    eprintln!("Warning: file not found: {}", p.display());
                }
                false
            }
        })
        .collect();

    if valid_files.is_empty() {
        eprintln!("Error: no valid input files found");
        std::process::exit(1);
    }

    // Create output directory
    if let Err(e) = fs::create_dir_all(&output_dir) {
        eprintln!("Error creating output directory: {}", e);
        std::process::exit(1);
    }

    let target_genes = if cli.extract_all {
        HashMap::new()
    } else {
        genes::parse_target_genes(&cli.genes)
    };

    let tracker = FilenameTracker::new();
    let log_path = output_dir.join("seqripper_log.tsv");

    // Process files in parallel with rayon
    let mut all_results: Vec<ExtractionResult> = valid_files
        .par_iter()
        .flat_map(|path| {
            process_file(
                path,
                cli.extract_all,
                cli.skip_trna,
                &target_genes,
                &output_dir,
                &tracker,
                cli.verbose,
            )
        })
        .collect();

    // Write log
    if let Err(e) = write_log(&log_path, &mut all_results) {
        eprintln!("Error writing log file: {}", e);
        std::process::exit(1);
    }

    // Summary
    if !cli.quiet {
        let extracted = all_results.iter().filter(|r| r.status == "extracted").count();
        let errors = all_results.iter().filter(|r| r.status == "error").count();

        eprintln!("SUMMARY:");
        eprintln!("Sequences extracted: {}", extracted);
        if errors > 0 {
            eprintln!("Errors: {}", errors);
        }
        eprintln!("Output written to: {}", output_dir.display());
    }
}
