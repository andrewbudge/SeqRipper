# SeqRipper: Extract genes from GenBank mitochondrial genomes

SeqRipper is a command-line tool that extracts specific genes from mitochondrial genomes using GenBank accession numbers. It reads annotation data directly from NCBI and outputs individual FASTA files for each gene-accession combination. SeqRipper directly interfaces with NCBI's API, so there is no need to wait for a separate database to update.

## Installation

**Requirements:** Python 3.7+, BioPython 1.80+

```bash
git clone https://github.com/andrewbudge/SeqRipper.git
cd SeqRipper
pip install -r requirements.txt
```

Then run SeqRipper with:

```bash
python SeqRipper.py -i INPUT -g GENES -o OUTPUT
```

## Usage

### Basic Syntax

```bash
python SeqRipper.py -i INPUT -g GENES [GENES ...] -o OUTPUT
python SeqRipper.py -i INPUT --all -o OUTPUT
```

### Required Arguments

1. **Input accession(s)**: Either a single GenBank accession number OR a text file with one accession per line
2. **Genes to extract**: Space-separated gene abbreviations (see Supported Genes below)
3. **Output directory**: A directory to save the extracted FASTA files

### Options

```
-i, --input       Single accession number OR file with accessions (one per line)
-o, --output      Output directory for FASTA files
-g, --genes       Genes to extract (space-separated: COI, 16S, ND2, etc.)
--all             Extract all annotated genes
--no-trna         Skip tRNA genes (only with --all)
-e, --email       Email for NCBI (or set NCBI_EMAIL env variable)
```

### Examples

### Extract specific genes from a single accession

```bash
python SeqRipper.py -i NC_012920.1 -g 12S 16S COI CYTB -o results/
```

### Extract genes from multiple accessions

Create a file `accessions.txt`:

```
NC_012920.1
NC_001807.4
NC_002008.4
```

Then run:

```bash
python SeqRipper.py -i accessions.txt -g COI ND2 -o results/
```
## Output Format

SeqRipper creates **one FASTA file per accession-gene combination**.

**Example output file** (`NC_012920.1_COI.fasta`):

```
>Homo_sapiens_NC_012920.1 [gene=COI]
ATGTTCGCCGACCGTTGACTATTCTCTACAAACCACAAAGACATTGGAACACTATACTAT
ATTCTTCGGTGCCTGAGCCGGAATAGTAGGCACGGCCCTAAGCCTCCTTATTCGAGCCGA
GCTAAGCCAACCCGGATCACTTCTAGGTAACGACCACATCTACAACGTTATCGTCACAGC
...
```

**File naming convention:** `{ACCESSION}_{GENE}.fasta`

## Supported Genes

SeqRipper can extract 15 common mitochondrial genes. It automatically handles multiple naming conventions:

|Gene|Abbreviations Recognized|
|---|---|
|**12S**|12S, 12S rRNA, s-rRNA, 12S ribosomal RNA, rrnS|
|**16S**|16S, 16S rRNA, l-rRNA, 16S ribosomal RNA, rrnL|
|**COI**|COI, COXI, COX1, CO1, cytochrome oxidase subunit 1, cytochrome c oxidase subunit I|
|**COII**|COII, COXII, COX2, CO2, cytochrome oxidase subunit 2, cytochrome c oxidase subunit II|
|**COIII**|COIII, COXIII, COX3, CO3, cytochrome oxidase subunit 3|
|**CYTB**|CYTB, Cyt b, cytochrome b, cob, CYB|
|**ND1**|ND1, NAD1, NADH1, NADH dehydrogenase subunit 1|
|**ND2**|ND2, NAD2, NADH2, NADH dehydrogenase subunit 2|
|**ND3**|ND3, NAD3, NADH3, NADH dehydrogenase subunit 3|
|**ND4**|ND4, NAD4, NADH4, NADH dehydrogenase subunit 4|
|**ND4L**|ND4L, NAD4L, NADH4L, NADH dehydrogenase subunit 4L|
|**ND5**|ND5, NAD5, NADH5, NADH dehydrogenase subunit 5|
|**ND6**|ND6, NAD6, NADH6, NADH dehydrogenase subunit 6|
|**ATP6**|ATP6, ATPase6, ATPase 6, ATP synthase F0 subunit 6|
|**ATP8**|ATP8, ATPase8, ATPase 8, ATP synthase F0 subunit 8|

**Note:** Gene names are case-insensitive. `COI`, `coi`, and `CoI` all work.

## Citation

If SeqRipper helps your research, please cite:

```
Budge, A. (2025). SeqRipper: Extract genes from GenBank mitochondrial genomes.
GitHub: https://github.com/andrewbudge/SeqRipper
```

## Contributing

Found a bug? Have a feature request? Open an issue at:  
https://github.com/andrewbudge/SeqRipper/issues

## License

MIT License - see LICENSE file for details.
