# SeqRipper

SeqRipper is a program to extract genes from mitocondiral genomes in genbank format (.gb). SeqRipper is written in golang, and can be complied to a single binary with no depedencies. 

## Installation

```bash
go install github.com/andrewbudge/seqripper@latest
```

Or build from source:
```bash
go build -o seqripper seqripper.go
```

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
