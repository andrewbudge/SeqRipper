# SeqRipper

Fast, concurrent tool for extracting genes from GenBank files. Written in Rust.

## Installation

```bash
cargo install --git https://github.com/andrewbudge/SeqRipper.git
```

Or build from source:
```bash
git clone https://github.com/andrewbudge/SeqRipper.git
cd SeqRipper
cargo build --release
```

The binary will be at `target/release/seqripper`.

## Usage

Extract specific genes:
```bash
seqripper -g COI,CYTB,12S,16S mitogenomes/*.gb output/
```

Extract all genes:
```bash
seqripper --all mitogenomes/*.gb output/
```

Extract all genes, skip tRNAs:
```bash
seqripper --all --no-trna mitogenomes/*.gb output/
```

## Options

- `-g GENES` - Comma-separated genes to extract
- `--all` - Extract all annotated genes
- `--no-trna` - Skip tRNA genes (with --all)
- `-v, --verbose` - Detailed output
- `-q, --quiet` - Suppress output

## Supported Genes

12S, 16S, COI, COII, COIII, CYTB, ND1-ND6, ND4L, ATP6, ATP8

Handles multiple naming conventions (COI/COXI/COX1/cytochrome oxidase subunit 1).

## Output

- Individual FASTA files: `<accession>_<gene>.fasta`
- Log file: `seqripper_log.tsv`

## Performance

Processes 34 mitochondrial genomes in ~13ms using parallel processing via rayon.
