# HLA-Splitter

**HLA-Splitter** is a computational framework leveraging HLA genotypes to demultiplex single-cell RNA-sequencing (scRNA-seq) data.

------

## Installation

We recommend using **Conda** or **Mamba** to manage environments and external bioinformatics dependencies.

### 1. Create Environment and Install Dependencies

```Bash
# Create a new environment
conda create -n hla_splitter_env python=3.9 -y
conda activate hla_splitter_env

# Install external bioinformatics tools (kallisto version =0.46.1)
conda install -c bioconda -c conda-forge samtools>=1.21 kallisto=0.46.1 bustools>=0.43.1 -y
    
# Mannual install Cell Ranger (version>=7.0.0)
# visit “https://www.10xgenomics.com/support/software/cell-ranger/downloads”
```

### 2. Install HLA-Splitter

```Bash
# Clone the repository
git clone https://github.com/Chenlab-ZJU/HLA-Splitter.git
cd HLA-Splitter

# Install in editable mode (recommended for development) or standard mode
pip install .
```

### 3. Run test data (optional)

```bash
# Run test data (Mixure of 30 PBMC smaples)
script="./Test_data/HLA-Splitter_test.sh"
chmod +x $script
bash $script
```



------

## Prepare HLA Database

HLA-Splitter requires the **IPD-IMGT/HLA database** to build a reference.

1. **Download**: [hla_nuc.fasta](https://www.ebi.ac.uk/ipd/imgt/hla/) (The nucleotide coding sequences/CDS).
2. **Version**: Compatible with version **3.60.0+**.
3. **Tip**: Ensure this file is updated periodically from the source to maintain accuracy.

```bash
mkdir -p ./hladb && cd ./hladb
# Download file via GitHub Repository
wget https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/fasta/hla_nuc.fasta
```



------

## Usage

You can run `HLA-Splitter` directly from your terminal:

```Bash
HLA-Splitter \
    --bam /path/to/possorted_genome_bam.bam \
    --barcodes /path/to/barcodes.tsv.gz \
    --hladb /path/to/hla_nuc.fasta \
    --hlalist /path/to/samples_hla_alleles.csv \
    --outdir /path/to/output_directory \
    --threads 16
```

### Parameters Detail

| **Argument** | **Full Name** | **Required** | **Description**                                              |
| ------------ | ------------- | ------------ | ------------------------------------------------------------ |
| `-B`         | `--bam`       | Yes          | Indexed BAM file (e.g., Cell Ranger's `possorted_genome_bam.bam`). |
| `-b`         | `--barcodes`  | Yes          | Cell barcodes file (e.g., `barcodes.tsv.gz`).                |
| `-H`         | `--hladb`     | Yes          | The IMGT/HLA nucleotide FASTA file (`hla_nuc.fasta`).        |
| `-L`         | `--hlalist`   | Yes          | CSV file listing HLA alleles for each donor (see format below). |
| `-o`         | `--outdir`    | Yes          | Path to the output directory.                                |
| `-t`         | `--threads`   | No           | Number of threads (default: 8).                              |
|              | --graph-corr  | No           | Enable graph-based correlation for demultiplexing refinement (default: False). |
|              | --version     | No           | Get current version.                                         |

------

## Input Requirements: `--hlalist`

The `--hlalist` CSV is crucial for accurate demultiplexing.

- **Precision**: HLA genotypes must be at **2-field precision** (e.g., `A*11:01`). Low precision may lead to poor resolution.

- **Format**: Sample names must be in the header row.

  **Example `hlalist.csv` Structure:**



|         | **Donor1** | **Donor2** | **Donor3** | **Donor4** |
| ------- | ---------- | ---------- | ---------- | ---------- |
| **A1**  | A*11:01    | A*11:01    | A*02:06    | A*11:01    |
| **A2**  | A*24:02    | A*24:02    | A*33:03    | A*24:02    |
| **B1**  | B*15:02    | B*40:01    | B*48:01    | B*13:01    |
| **...** | ...        | ...        | ...        | ...        |

------

## Output Structure

After execution, the following directory structure is generated:

Plaintext

```
output_directory/
├── sc_output/             # Alignment results (Kallisto/Bustools)
│   └── filter_counts/     # Filtered count matrices
│       ├── HLA_type.mtx   # Cell count matrix of HLA alleles
│       └── output.barcodes.txt
├── sorted_genotypes.txt   # Sorted HLA alleles list
├── HLA_demultiplex.csv    # Final demultiplexing results (Barcode -> Donor)
├── HLA_TSNE_Adjusted_Prediction.pdf # Visualization of cell clusters
└── Countplot_Adjusted_Prediction.pdf # Visualization of numbers of sample labels
```

**Example of `HLA_demultiplex.csv`:**

| **Cell barcode** | **Sample_prediction** | **Certainty** |
| ---------------- | --------------------- | ------------- |
| AAACCTGAGAAGGCCT | Donor30               | 0.91          |
| AAACCTGAGCCAGTAG | Donor28               | 0.93          |
| AAACCTGCAAGCCGCT | Donor15               | 0.20          |

------

## How it Works

1. **Read Extraction**: Uses `samtools` to extract reads mapped to the HLA region (chr6).
2. **Pseudo-alignment**: Leverages `kallisto` and `bustools` to align reads against the donor-specific HLA reference.
3. **Bayesian Assignment**: `analysis.py` calculates the probability of each cell belonging to a donor based on allele-specific counts.
4. **Logging**: Full progress and error tracking are handled via the Python logging module.

------

### Citation

If you use HLA-Splitter in your research, please cite:

*Chen Lab, Zhejiang University (2024)*.

