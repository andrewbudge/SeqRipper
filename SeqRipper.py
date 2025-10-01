#!/usr/bin/env python3

import argparse
import sys
import os
from pathlib import Path
from Bio import Entrez, SeqIO
from typing import List, Dict, Optional
import time

__version__ = "1.0.0"
__author__ = "Andrew Budge"

# Config data
GENE_SYNONYMS = {
    "12S": ["12S", "12S rRNA", "s-rRNA", "12S ribosomal RNA", "rrnS"],
    "16S": ["16S", "16S rRNA", "l-rRNA", "16S ribosomal RNA", "rrnL"],
    "COI": ["COI", "COXI", "COX1", "CO1", "cytochrome oxidase subunit 1", 
            "cytochrome c oxidase subunit I", "cytochrome oxidase subunit I"],
    "COII": ["COII", "COXII", "COX2", "CO2", "cytochrome oxidase subunit 2",
             "cytochrome c oxidase subunit II"],
    "COIII": ["COIII", "COXIII", "COX3", "CO3", "cytochrome oxidase subunit 3"],
    "CYTB": ["CYTB", "Cyt b", "cytochrome b", "cob", "CYB"],
    "ND1": ["ND1", "NAD1", "NADH1", "NADH dehydrogenase subunit 1"],
    "ND2": ["ND2", "NAD2", "NADH2", "NADH dehydrogenase subunit 2"],
    "ND3": ["ND3", "NAD3", "NADH3", "NADH dehydrogenase subunit 3"],
    "ND4": ["ND4", "NAD4", "NADH4", "NADH dehydrogenase subunit 4"],
    "ND4L": ["ND4L", "NAD4L", "NADH4L", "NADH dehydrogenase subunit 4L"],
    "ND5": ["ND5", "NAD5", "NADH5", "NADH dehydrogenase subunit 5"],
    "ND6": ["ND6", "NAD6", "NADH6", "NADH dehydrogenase subunit 6"],
    "ATP6": ["ATP6", "ATPase6", "ATPase 6", "ATP synthase F0 subunit 6"],
    "ATP8": ["ATP8", "ATPase8", "ATPase 8", "ATP synthase F0 subunit 8"],
}

FEATURE_TYPES = ["rRNA", "CDS", "gene", "tRNA", "misc_RNA"]
RATE_LIMIT = 0.34  # this is to keep NCBI happy

# Functions

def clean_name(name: str) -> str:
    """Clean a name for use in filenames"""
    return name.replace(" ", "_").replace("/", "_")


def get_gene_name(feature) -> Optional[str]:
    """Extract gene name from a GenBank feature"""
    if "gene" in feature.qualifiers:
        return feature.qualifiers["gene"][0]
    elif "product" in feature.qualifiers:
        return feature.qualifiers["product"][0]
    return None


def match_gene(annotation: str, target_genes: Dict[str, List[str]]) -> Optional[str]:
    """
    Match an annotation name to target genes using synonyms
    Returns the matched gene name or None
    """
    annotation_lower = annotation.lower().strip()
    
    for gene, synonyms in target_genes.items():
        # Try exact match first
        for syn in synonyms:
            if annotation_lower == syn.lower():
                return gene
        # Then try partial match
        for syn in synonyms:
            if syn.lower() in annotation_lower:
                return gene
    
    return None


def write_fasta(output_file: Path, header: str, sequence: str) -> None:
    """Write sequence data to the FASTA"""
    with open(output_file, 'w') as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")


def load_accessions(input_path: str) -> List[str]:
    """
    Read in the Acc. numbers. Either a single or a list
    """
    if os.path.isfile(input_path):
        print(f"Reading accessions from file: {input_path}")
        with open(input_path, 'r') as f:
            accessions = [line.strip() for line in f 
                         if line.strip() and not line.startswith('#')]
        print(f"Found {len(accessions)} accession(s)")
        return accessions
    
    # Single accession
    print(f"Processing single accession: {input_path}")
    return [input_path]


def parse_genes(gene_list: List[str]) -> Dict[str, List[str]]:
    """Build gene dictionary with synonyms from user input"""
    gene_dict = {}
    
    for gene in gene_list:
        gene_upper = gene.upper()
        
        if gene_upper in GENE_SYNONYMS:
            gene_dict[gene_upper] = GENE_SYNONYMS[gene_upper]
            print(f"  {gene_upper}: using {len(GENE_SYNONYMS[gene_upper])} synonyms")
        else:
            # Custom gene - create basic synonyms
            gene_dict[gene_upper] = [gene, gene.upper(), gene.lower()]
            print(f"  {gene_upper}: custom gene")
    
    return gene_dict


# Main extraction function

def extract_genes(accession: str, target_genes: Dict[str, List[str]], 
                 output_dir: Path, email: str, extract_all: bool = False,
                 skip_trna: bool = False) -> Dict[str, bool]:
    """Extract genes from a GenBank accession"""
    Entrez.email = email
    results = {} if extract_all else {gene: False for gene in target_genes.keys()}
    
    print(f"\n{'='*60}")
    print(f"Processing: {accession}")
    print(f"{'='*60}")
    
   
    try:
        print("  Fetching from NCBI...", end=' ')
        handle = Entrez.efetch(db="nuccore", id=accession, 
                              rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()
        print("✓")
        
    except ValueError:
        print("✗")
        print(f"  ✗ ERROR: Accession '{accession}' not found in NCBI")
        print(f"     Check the accession number and try again")
        return results
        
    except IOError:
        print("✗")
        print(f"  ✗ ERROR: Network error while fetching {accession}")
        print(f"     Check your internet connection")
        return results
        
    except Exception:
        print("✗")
        print(f"  ✗ ERROR: Failed to fetch {accession}")
        return results
    
    try:
        species = clean_name(record.annotations.get("organism", "Unknown"))
        print(f"  Species: {species}")
        
        mode = "all genes" if extract_all else "specified genes"
        suffix = " (excluding tRNAs)" if skip_trna else ""
        print(f"  Extracting {mode}{suffix}...")
        
        found_genes = set()
        
        for feature in record.features:
            if feature.type not in FEATURE_TYPES:
                continue
            
            if skip_trna and feature.type == "tRNA":
                continue
            
            gene_name = get_gene_name(feature)
            if not gene_name:
                continue
            
            if extract_all:
                target_name = clean_name(gene_name)
            else:
                target_name = match_gene(gene_name, target_genes)
            
            if not target_name or target_name in found_genes:
                continue
            
            sequence = str(feature.extract(record.seq))
            filename = f"{accession}_{target_name}.fasta"
            output_file = output_dir / filename
            
            # Write fasta out with new header
            # Get start and end bp position
            start = int(feature.location.start) + 1
            end = int(feature.location.end)

            # Get organism name
            organism = record.annotations.get("organism", "Unknown organism")

            # Create new header
            header = f"{accession}:{start}-{end} {organism} [gene={target_name}]"
            # Write fasta with extracted sequnce and new header
            write_fasta(output_file, header, sequence)

            print(f"    ✓ {target_name:20s} → {filename} ({len(sequence)} bp)")
            found_genes.add(target_name)
            results[target_name] = True
        
        if extract_all:
            print(f"  Total: {len(found_genes)} genes extracted")
        else:
            missing = set(target_genes.keys()) - found_genes
            if missing:
                print(f"    ✗ Missing: {', '.join(sorted(missing))}")
        
        time.sleep(RATE_LIMIT)
        return results
        
    except Exception as e:
        # CHANGE 6: Better error for parsing failures
        print(f"  ✗ ERROR: Failed to extract genes from {accession}: {e}")
        return results


# CLI details and handeling

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        prog='seqripper',
        description='SeqRipper - Extract genes from GenBank sequences',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract specific genes
  seqripper -i PQ283265.2 -g 12S 16S COI -o results/
  
  # Extract from multiple accessions
  seqripper -i accessions.txt -g COI CYTB ND2 -o results/
  
  # Extract ALL genes
  seqripper -i PQ283265.2 --all -o results/
  
  # Extract ALL except tRNAs
  seqripper -i accessions.txt --all --no-trna -o results/

Output:
  Creates files: ACCESSION_GENE.fasta
  Example: PQ283265_COI.fasta, PQ283265_12S.fasta
        """)
    
    parser.add_argument('-i', '--input', required=True,
                       help='Single accession OR file with accessions (one per line)')
    
    parser.add_argument('-o', '--output', required=True,
                       help='Output directory')
    
    gene_group = parser.add_mutually_exclusive_group(required=True)
    gene_group.add_argument('-g', '--genes', nargs='+',
                           help='Genes to extract (space-separated)')
    gene_group.add_argument('--all', action='store_true',
                           help='Extract ALL annotated genes')
    
    parser.add_argument('--no-trna', action='store_true',
                       help='Skip tRNA genes (only with --all)')
    
    parser.add_argument('-e', '--email', 
                       default=os.environ.get('NCBI_EMAIL', 'anonymous@example.com'),
                       help='Email for NCBI (or set NCBI_EMAIL env variable)')
    
    parser.add_argument('--version', action='version',
                       version=f'SeqRipper {__version__}')
    
    return parser.parse_args()


def main() -> None:
    """Main entry point"""
    args = parse_args()
    
    # Setup
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Print header
    print("="*60)
    print(f"SeqRipper v{__version__}")
    print("="*60)
    
    # Load inputs
    try:
        accessions = load_accessions(args.input)
    except Exception as e:
        print(f"ERROR loading accessions: {e}")
        sys.exit(1)
    
    if args.all:
        print(f"\nMode: Extract ALL genes" + 
              (" (excluding tRNAs)" if args.no_trna else ""))
        target_genes = {}
    else:
        if args.no_trna:
            print("\nWarning: --no-trna only works with --all mode")
        print("\nGenes to extract:")
        target_genes = parse_genes(args.genes)
    
    print(f"\nOutput: {output_dir.absolute()}")
    
    # Process accessions
    print("\n" + "="*60)
    print("PROCESSING")
    print("="*60)
    
    summary = {}
    for i, accession in enumerate(accessions, 1):
        print(f"\n[{i}/{len(accessions)}]")
        results = extract_genes(accession, target_genes, output_dir, 
                               args.email, args.all, args.no_trna)
        summary[accession] = results
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    total = 0
    for accession, results in summary.items():
        count = sum(results.values())
        total += count
        print(f"{accession}: {count} genes")
    
    print(f"\nTotal sequences: {total}")
    print(f"Output directory: {output_dir.absolute()}")
    print("\n✓ Done!")


if __name__ == "__main__":
    main()
