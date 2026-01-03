# HLA-Splitter

**Description:** A computational framework leveraging HLA genotypes to demultiplex scRNA-seq data.

## External Dependencies

The following external tools must be installed and accessible in your system PATH:

* `samtools` (>= v1.21)
* `cellranger` (>= v7.0.0) - *10x Genomics*
* `kallisto` (>= v0.46.1)
* `bustools` (>= v0.43.1)

> **Note:** We recommend installing these tools via Conda/Mamba.

## Python Dependencies

Python dependencies are listed in `pyproject.toml`. Installing this package via `pip` will handle these dependencies automatically.

## HLA Database

**IPD-IMGT/HLA database** is a specialist sequence database for sequences of the human histocompatibility complex.

We use the file `hla_nuc.fasta` which contains the nucleotide coding sequences (CDS) to build the HLA reference.

* **Download:** [hla_nuc.fasta](https://github.com/ANHIG/IMGTHLA/blob/Latest/hla_nuc.fasta)
* **Version:** The logic is compatible with version 3.60.0+.
* **Tip:** It is best to check and update this file from the source before use.

## Installation

```bash
# Clone the repository
git clone [https://github.com/Chenlab-ZJU/HLA-Splitter.git](https://github.com/Chenlab-ZJU/HLA-Splitter.git)

# Enter the directory
cd HLA-Splitter

# Install the package
pip install .

# Or for development mode:
# pip install -e .
```

## Usage

```Bash
HLA-Splitter \
    --bam /path/to/your/possorted_genome_bam.bam \
    --barcodes /path/to/your/filtered_feature_bc_matrix/barcodes.tsv.gz \
    --hladb /path/to/IMGTHLA_latest/hla_nuc.fasta \
    --hlalist /path/to/your/samples_hla_alleles.csv \
    --outdir /path/to/output_directory \
    --threads 16
```

## Input Data Requirements

- **`-B, --bam`**: Required. The indexed BAM file (e.g., `possorted_genome_bam.bam` from Cell Ranger output).
- **`-b, --barcodes`**: Required. The cell barcodes file (e.g., `filtered_feature_bc_matrix/barcodes.tsv.gz` from Cell Ranger).
- **`-H, --hladb`**: Required. The IMGT/HLA nucleotide FASTA file (e.g., `IMGTHLA_latest/hla_nuc.fasta`).
- **`-L, --hlalist`**: Required. A CSV file listing HLA alleles for each expected sample (Donor).
  - **Header row:** Must contain sample names.
  - **Precision:** The HLA genotype needs to be of **2-field precision** (e.g., `A*11:01`). Low-precision HLA genes may lead to poor resolution results.

**Example format for `--hlalist`:**

```csv
        Donor1      Donor2      Donor3      Donor4
A1      A*11:01     A*11:01     A*02:06     A*11:01
A2      A*24:02     A*24:02     A*33:03     A*24:02
B1      B*15:02     B*40:01     B*48:01     B*13:01
B2      B*51:01     B*55:02     B*58:01     B*15:01
C1      C*08:01     C*07:02     C*03:02     C*03:04
C2      C*14:02     C*12:03     C*08:03     C*07:02
DRB11   DRB1*04:05  DRB1*11:01  DRB1*12:01  DRB1*09:01
DRB12   DRB1*12:02  DRB1*15:01  DRB1*13:02  DRB1*15:01
DQB11   DQB1*03:01  DQB1*06:01  DQB1*06:09  DQB1*06:01
DQB12   DQB1*04:01  DQB1*03:01  DQB1*03:01  DQB1*03:03
```

- **`-o, --outdir`**: Path to the output directory.
- **`-t, --threads`**: Number of threads to use for compatible tools (default: 8).

## Output

The output directory contains the following structure and main files:

```Plaintext
output_directory/
├── fastq/                             # Extracted chr6 reads converted to FASTQ
├── sc_output/                         # Alignment results using Kallisto
│   └── filter_counts/                 # Filtered results
│       ├── HLA_type.mtx               # Cell count matrix of HLA alleles
│       └── output.barcodes.txt        # Cell barcodes
├── sorted_genotypes.txt               # Sorted HLA alleles
├── HLA_demultiplex.csv                # Final demultiplexing result CSV
├── HLA_TSNE_Adjusted_Prediction.pdf   # t-SNE plot of demultiplexed cells
└── Countplot_Adjusted_Prediction.pdf  # Count plot of demultiplexed cells
```

**Example of `HLA_demultiplex.csv`:**

| **Cell barcode** | **Prediction** | **Adjusted_prediction** | **Certainty** |
| ---------------- | -------------- | ----------------------- | ------------- |
| AAACCTGAGAAGGCCT | Donor30        | Donor30                 | 0.91          |
| AAACCTGAGACCACGA | Donor1         | Donor1                  | 0.98          |
| AAACCTGAGACTGTAA | Donor18        | Donor18                 | 0.95          |
| AAACCTGAGATCGGGT | Donor3         | Donor3                  | 0.99          |
| AAACCTGAGATCTGCT | Donor23        | Donor23                 | 0.99          |
| AAACCTGAGCCAGTAG | Donor28        | Donor28                 | 0.93          |
| AAACCTGCAAGCCGCT | Donor15        | Undefined               | 0.2           |

**Example Plot (`HLA_TSNE_Adjusted_Prediction.pdf`):**

![image](https://github.com/user-attachments/assets/b1476531-e4ae-4379-bb3b-7972e3a8e36f)

**Example Plot (`Countplot_Adjusted_Prediction.pdf`):**

- ![image](https://github.com/user-attachments/assets/a0998765-7c0e-40f5-9f3b-b2c311fd3802)

## How it Works

1. **Installation:** The user installs the package using `pip`.
2. **Command Availability:** This makes the `HLA-Splitter` command available (defined in `pyproject.toml` under `[project.scripts]`).
3. **Execution:** Running `HLA-Splitter` executes the `main()` function in `hla_demultiplex/cli.py`.
4. **Workflow:**
   - `cli.py` parses arguments and calls functions in `steps.py` and `analysis.py`.
   - `steps.py` uses `utils.run_command` to execute the necessary external bioinformatics tools (`samtools`, `kallisto`, etc.).
   - `analysis.py` handles the logic for assigning cells to donors based on the alignment results.
5. **Logging:** The tool provides logging to display progress information and error details.
