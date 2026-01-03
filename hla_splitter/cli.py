import argparse
import os
import sys
import logging
from hla_splitter import analysis
from hla_splitter import steps
from hla_splitter import __version__ # Import version

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
log = logging.getLogger(__name__)

def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description=f"HLA-Splitter Pipeline v{__version__}. Processes single-cell BAM files for HLA typing and sample demultiplexing."
    )
    parser.add_argument("-B", "--bam", required=True, help="Path to input BAM file (indexed).")
    parser.add_argument("-b", "--barcodes", required=True, help="Path to cell barcodes file (e.g., filtered_feature_bc_matrix/barcodes.tsv.gz from Cell Ranger).")
    parser.add_argument("-H", "--hladb", required=True, help="Path to IMGT/HLA nucleotide FASTA file (e.g., hla_nuc.fasta).")
    parser.add_argument("-L", "--hlalist", required=True, help="Path to CSV file listing HLA alleles for each expected sample (header row = sample names).")
    parser.add_argument("-o", "--outdir", required=True, help="Path to output directory.")
    parser.add_argument("-t", "--threads", type=int, default=8, help="Number of threads to use for compatible tools (default: 8).")
    parser.add_argument('--version', action='version', version=f'%(prog)s {__version__}')

    return parser.parse_args()

def main():
    """Main entry point for the HLA-Splitter pipeline."""
    args = get_args()

    log.info(f"Starting HLA-Splitter Pipeline v{__version__}")
    log.info(f"Output directory: {args.outdir}")
    log.info(f"Threads: {args.threads}")

    # --- Input Validation ---
    if not os.path.exists(args.bam):
        log.error(f"Input BAM file not found: {args.bam}")
        sys.exit(1)
    # Add index check? samtools requires .bai or .csi
    bam_index1 = args.bam + ".bai"
    bam_index2 = args.bam.replace(".bam", ".bai") # Handle cases where index might have different name convention
    if not (os.path.exists(bam_index1) or os.path.exists(bam_index2)):
         log.warning(f"BAM index (.bai) not found for {args.bam}. Indexing is required.")
         # Optionally run samtools index here, or just exit
         # sys.exit(1)

    if not os.path.exists(args.barcodes):
        log.error(f"Barcodes file not found: {args.barcodes}")
        sys.exit(1)
    if not os.path.exists(args.hladb):
        log.error(f"HLA FASTA database not found: {args.hladb}")
        sys.exit(1)
    if not os.path.exists(args.hlalist):
        log.error(f"HLA list CSV file not found: {args.hlalist}")
        sys.exit(1)

    # Create output directory
    try:
        os.makedirs(args.outdir, exist_ok=True)
        log.info(f"Ensured output directory exists: {args.outdir}")
    except OSError as e:
        log.error(f"Could not create output directory {args.outdir}: {e}")
        sys.exit(1)


    # --- Pipeline Steps ---
    try:
        log.info("Step 1: Processing HLA sample list...")
        sorted_genotypes_file = analysis.prepare_hla_inputs(args.hlalist, args.outdir)

        log.info("Step 2: Extracting chr6 reads and converting to FASTQ...")
        fastq_dir = steps.extract_chr6_reads(args.bam, args.outdir, args.threads)

        log.info("Step 3: Building HLA reference index...")
        kallisto_index = steps.build_hla_reference(args.hladb, sorted_genotypes_file, args.outdir)

        log.info("Step 4: Running alignment and counting...")
        count_dir = steps.run_alignment_and_counting(args.barcodes, kallisto_index, fastq_dir, args.outdir, args.threads)

        log.info("Step 5: Summarizing counts per HLA type...")
        hla_matrix_file = analysis.summarize_counts(args.outdir) # Uses files within count_dir

        log.info("Step 6: Running demultiplexing analysis...")
        demux_results_file = analysis.run_demultiplex(args.hlalist, args.outdir) # Uses hla_matrix_file

        log.info("--- Pipeline Finished Successfully ---")
        log.info(f"Final demultiplexing results: {demux_results_file}")

    except FileNotFoundError:
         log.error("A required file was not found during pipeline execution. Please check logs.")
         sys.exit(1)
    except subprocess.CalledProcessError:
         log.error("A command-line tool failed. Please check logs for details.")
         sys.exit(1)
    except Exception as e:
         log.error(f"An unexpected error occurred during the pipeline: {e}")
         import traceback
         traceback.print_exc()
         sys.exit(1)


if __name__ == "__main__":
    # This allows running python -m hla_splitter.cli if needed
    main()