- **HLA-Splitter**

- **Description:** a computational framework leveraging HLA genotypes demultiplex scRNA-seq data.

- **External Dependencies:**

  - `samtools` (>\= v1.21)

  - `cellranger` (>\= v7.0.0) - Mention it's from 10x Genomics.

  - `kallisto` (>\= v0.46.1)

  - `bustools` (>\= v0.43.1)

  - Recommend installation via Conda for these tools.

- **Python Dependencies:**&#x20;

  - Listed in `pyproject.toml`. Installation via `pip` will handle these.

- \*\*HLA database: \*\*IPD-IMGT/HLA database which is is a specialist sequence database for sequences of the human histocompatibility complex. We use the file “[hla_nuc.fasta](https://github.com/ANHIG/IMGTHLA/blob/Latest/hla_nuc.fasta 'hla_nuc.fasta')” which contain the nucleotide coding sequences (CDS) to build HLA reference. The attached file version is 3.60.0. It is best to check and update this file before use.

- **Installation:**

  Bash

  ```
  git clone https://github.com/Yufei/HLA-Splitter_v1.git
  cd HLA-Splitter_v1
  pip install .
  # Or for development: pip install -e .

  ```

- **Usage:**

  Bash

  ```
  HLA-Splitter \
      --bam /path/to/your/possorted_genome_bam.bam \
      --barcodes /path/to/your/filtered_feature_bc_matrix/barcodes.tsv.gz \
      --hladb /path/to/IMGTHLA_latest/hla_nuc.fasta \
      --hlalist /path/to/your/samples_hla_alleles.csv \
      --outdir /path/to/output_directory \
      --threads 16

  ```

- **Input Data:**

  - -B, --bam BAM: required BAM (indexed) file (e.g., possorted_genome_bam.bam from Cell Ranger).

  - -b, --barcodes BARCODES: required cell barcodes file (e.g., filtered_feature_bc_matrix/barcodes.tsv.gz from Cell Ranger).

  - -H , --hladb HLADBE: required IMGT/HLA nucleotide FASTA file (e.g., IMGTHLA_latest/hla_nuc.fasta).

  - -L , --hlalist HLALIST: required CSV file listing HLA alleles for each expected sample (header row \= sample names). The HLA genotype needs to be of 2-filed precision. Low-precision HLA genes may lead to poor resolution results.&#x20;

  - HLA_demultiplex.csv will look like

    ```HLALIST
    	Donor1	Donor2	Donor3	Donor4
    A1	A*11:01	A*11:01	A*02:06	A*11:01
    A2	A*24:02	A*24:02	A*33:03	A*24:02
    B1	B*15:02	B*40:01	B*48:01	B*13:01
    B2	B*51:01	B*55:02	B*58:01	B*15:01
    C1	C*08:01	C*07:02	C*03:02	C*03:04
    C2	C*14:02	C*12:03	C*08:03	C*07:02
    DRB11	DRB1*04:05	DRB1*11:01	DRB1*12:01	DRB1*09:01
    DRB12	DRB1*12:02	DRB1*15:01	DRB1*13:02	DRB1*15:01
    DQB11	DQB1*03:01	DQB1*06:01	DQB1*06:09	DQB1*06:01
    DQB12	DQB1*04:01	DQB1*03:01	DQB1*03:01	DQB1*03:03

    ```

  - -o , --outdir OUTDIR: Path to output directory.

  - -t , --threads THREADS: Number of threads to use for compatible tools (default: 8).

- **Output:** Describe the main output files (e.g., `HLA_demultiplex.csv`, count matrices).The HLA genotype needs to be of two-dimensional precision. Low-precision HLA genes may lead to poor resolution results

  ```output_directory
  output_directory/
  ├── fastq \Extracting chr6 reads and converting to FASTQ
  ├── sc_output \Alignment using Kallisto
  │   ├── filter_counts \filtered
  │       ├──HLA_type.mtx \cell count matrix of hla_alleles
  │       ├──output.barcodes.txt \cell barcodes
  ├── sorted_genotypes.txt \hla_alleles
  ├── HLA_demultiplex.csv \demultiplex result csv
  ├── HLA_TSNE_Adjusted_Prediction.pdf \TSNE plot of demultiplex cells
  ├── Countplot_Adjusted_Prediction.pdf \Count plot of demultiplex cells
  │ ...
  ```

- HLA_demultiplex.csv will look like

  ```HLA_demultiplex.csv
  Cell barcode	Prediction	Adjusted_prediction	Certainty
  AAACCTGAGAAGGCCT	Donor30	Donor30	0.91
  AAACCTGAGACCACGA	Donor1	Donor1	0.98
  AAACCTGAGACTGTAA	Donor18	Donor18	0.95
  AAACCTGAGATCGGGT	Donor3	Donor3	0.99
  AAACCTGAGATCTGCT	Donor23	Donor23	0.99
  AAACCTGAGCCAGTAG	Donor28	Donor28	0.93
  AAACCTGAGCGACGTA	Donor1	Donor1	0.91
  AAACCTGAGCTACCGC	Donor2	Donor2	0.88
  AAACCTGAGCTGAACG	Donor25	Donor2	0.83
  AAACCTGAGTTACCCA	Donor18	Donor18	0.92
  AAACCTGCAAGACGTG	Donor10	Donor10	0.9
  AAACCTGCAAGCCGCT	Donor15	Undefined	0.2
  AAACCTGCACAGCCCA	Donor29	Donor29	0.98
  AAACCTGCACTAGTAC	Donor2	Donor2	0.82
  AAACCTGCACTGAAGG	Donor30	Donor30	0.93
  AAACCTGGTCTCACCT	Donor15	Donor15	0.87
  AAACCTGGTGAGGCTA	Donor3	Donor3	0.98
  AAACCTGGTGTAATGA	Donor4	Donor4	0.89
  AAACCTGGTTAAGGGC	Donor8	Donor8	0.91
  AAACCTGGTTCAACCA	Donor5	Donor5	0.96
  AAACCTGGTTGGACCC	Donor10	Donor10	0.75

  ```

- HLA_TSNE_Adjusted_Prediction.pdf will look like

![image](https://github.com/user-attachments/assets/b1476531-e4ae-4379-bb3b-7972e3a8e36f)


**How it Works:**

- The user installs the package using `pip`.

- This makes the HLA-Splitter command available (defined in `pyproject.toml`'s `[project.scripts]`).

- Running HLA-Splitter executes the `main()` function in `hla_demultiplex/cli.py`.

- `cli.py` parses arguments and calls functions in `steps.py` and `analysis.py`.

- `steps.py` uses `utils.run_command` to execute the necessary external bioinformatics tools (`samtools`, `kallisto`, etc.).

- `analysis.py` contains the Python logic originally in the separate `.py` scripts.

- Logging provides progress information and error details.

######
