import os
import sys
import logging
import gzip
from .utils import run_command # Use the helper function

log = logging.getLogger(__name__)

def extract_chr6_reads(bam_path, out_dir, threads=8):
    """
    Extracts chr6 reads, sorts, and converts to FASTQ.
    Equivalent to ExtractReads.sh.
    Requires samtools and cellranger.
    """
    log.info(f"Extracting chr6 reads from {bam_path}")
    hla_bam = os.path.join(out_dir, "hla.bam")
    hla_sorted_bam = os.path.join(out_dir, "hla.sorted.bam")
    fastq_dir = os.path.join(out_dir, "fastq")

    # Check for cellranger
    try:
        run_command("cellranger -h", check=False) # Check if cellranger is callable
    except FileNotFoundError:
         log.error("`cellranger` command not found. Please ensure 10x Cell Ranger is installed and in PATH.")
         sys.exit(1)

    # Extract chr6
    cmd1 = f"samtools view -@{threads} -b {bam_path} chr6 -o {hla_bam}"
    run_command(cmd1)

    # bamtofastq
    cmd2 = f"cellranger bamtofastq --nthreads={threads} {hla_bam} {fastq_dir}"
    run_command(cmd2) # Cellranger might output useful info to stderr

    log.info(f"FASTQ files generated in {fastq_dir}")
    return fastq_dir

def build_hla_reference(hla_nuc_fasta_path, sorted_genotypes_path, out_dir):
    """
    Builds kallisto index for selected HLA alleles.
    Equivalent to BuildReference.sh.
    Requires grep, samtools, sed, kallisto.
    """
    log.info(f"Building HLA reference index using {hla_nuc_fasta_path}")
    tmp_allele_file = os.path.join(out_dir, "tmpallele.txt")
    cds_fasta = os.path.join(out_dir, "cds.fasta")
    kallisto_index = os.path.join(out_dir, "hla_kallisto_index") # Changed name

    # 1. Grep selected alleles headers (Handle potential errors)
    log.info("Extracting allele headers...")
    # Ensure sorted_genotypes_path exists
    if not os.path.exists(sorted_genotypes_path):
        log.error(f"Sorted genotypes file not found: {sorted_genotypes_path}")
        sys.exit(1)
    # Using Python to avoid complex grep loop
    selected_headers = []
    # Store full headers indexed by their short ID
    full_header_map = {} 

    with open(sorted_genotypes_path, 'r') as sg:
        patterns = [line.strip() for line in sg if line.strip()]
    if not patterns:
         log.error("No genotypes found in sorted_genotypes.txt")
         sys.exit(1)
    # Instead of grep, read fasta header and check if pattern matches start
    # This is safer than direct grep if patterns have special characters
    with open(hla_nuc_fasta_path, 'r') as hla_fasta, open(tmp_allele_file, 'w') as tmp_out:
         current_header = None
         for line in hla_fasta:
             if line.startswith('>'):
                 current_header = line.strip()
                 # Check if any pattern matches the start of the header's allele part
                 header_allele_part = current_header.split(' ')[1][0:] # Get part of HLA allele
                 if (":" in header_allele_part):
                     header_allele_part = header_allele_part.split(":")[0][0:]+":"+header_allele_part.split(":")[1][0:] # 2 field HLA-genes
                 if any(header_allele_part==p for p in patterns):
                     tmp_out.write(current_header + '\n')
                     selected_headers.append(current_header)
                     short_id = current_header[1:].split(' ')[0]
                     full_header_map[short_id] = current_header # Map short ID to full header
             # No need to write sequence lines here
    log.info(f"Found {len(selected_headers)} matching allele headers.")
    if not selected_headers:
        log.error("No headers matching genotypes found in HLA FASTA. Check genotype format and FASTA.")
        sys.exit(1)


    # 2. Extract sequences using samtools faidx
    log.info("Extracting allele sequences using samtools faidx...")
    # Get just the IDs (>ID part) for faidx
    ids_for_faidx = [h.split(' ')[0][1:] for h in selected_headers]
    ids_str = " ".join(ids_for_faidx)
    cmd_faidx = f"samtools faidx {hla_nuc_fasta_path} {ids_str} -o {cds_fasta}"
    run_command(cmd_faidx)

    # 3. Re-add full description to headers in cds.fasta
    log.info("Modifying headers in cds.fasta to include full descriptions...")
    temp_cds_fasta = cds_fasta + ".tmp"
    with open(cds_fasta, 'r') as infile, open(temp_cds_fasta, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                # samtools faidx truncates headers, so line will look like ">HLA:A*01:01:01"
                # We need to extract this short ID to look up the original full header.
                short_id_from_faidx_line = line.strip()[1:] # Remove '>'
                
                # Retrieve the full original header from the map we built earlier
                full_original_header = full_header_map.get(short_id_from_faidx_line)
                
                if full_original_header:
                    outfile.write(full_original_header + '\n') # Write the full original header
                else:
                    # Fallback: if for some reason the short ID isn't in our map, write the line as is.
                    # This should ideally not happen if everything is consistent.
                    log.warning(f"Full original header not found for short ID '{short_id_from_faidx_line}'. Writing original header from samtools output.")
                    outfile.write(line)
            else:
                outfile.write(line) # Write sequence lines as is

    # Replace the original cds.fasta with the modified one
    os.replace(temp_cds_fasta, cds_fasta)
    log.info(f"Headers in {cds_fasta} modified successfully.")

    # 4. Build kallisto index
    log.info(f"Building kallisto index at {kallisto_index}")
    cmd_kallisto = f"kallisto index -i {kallisto_index} {cds_fasta}"
    run_command(cmd_kallisto)

    log.info("HLA reference index built successfully.")
    return kallisto_index

def get_all_file_paths(folder_path):
    """
    Gets the full path to all files in the specified folder and its subfolders.
    Args:
        folder_path (str): Path of the folder to search.
    Returns:
        list: A list containing the full paths to all files.
    """
    all_file_paths = []
    for root, dirs, files in os.walk(folder_path):
        for filename in files:
            file_path = os.path.join(root, filename)
            all_file_paths.append(file_path)
    return all_file_paths

def run_alignment_and_counting(barcodes_10x_path, kallisto_index_path, fastq_dir, out_dir, threads=8):
    """
    Runs kallisto bus and bustools commands.
    Equivalent to Allignment.sh.
    Requires kallisto, bustools, gunzip, sed (or Python equivalent).
    """
    log.info("Starting alignment (kallisto bus) and counting (bustools)...")
    sc_output_dir = os.path.join(out_dir, "sc_output")
    os.makedirs(sc_output_dir, exist_ok=True) # Ensure output dir exists

    # 1. Integrate FASTQ files (handle potential multiple lanes/files)
    r1_files = sorted([f for f in get_all_file_paths(out_dir) if '_R1_' in f and f.endswith('.fastq.gz')])
    r2_files = sorted([f for f in get_all_file_paths(out_dir) if '_R2_' in f and f.endswith('.fastq.gz')])

    if not r1_files or not r2_files:
        log.error(f"No R1 or R2 FASTQ files found in {fastq_dir}")
        sys.exit(1)

    # We merge multiple files directly to kallisto bus
    log.info("Merging R1 FASTQ files...")
    merged_r1_path = os.path.join(out_dir,"fastq/","merged_r1.fastq.gz")
    cat_r1_cmd = f"cat {' '.join(r1_files)} > {merged_r1_path}"
    run_command(cat_r1_cmd, shell=True)
    log.info(f"Merged R1 files into: {merged_r1_path}")

    log.info("Merging R2 FASTQ files...")
    merged_r2_path = os.path.join(out_dir,"fastq/","merged_r2.fastq.gz")
    cat_r2_cmd = f"cat {' '.join(r2_files)} > {merged_r2_path}"
    run_command(cat_r2_cmd, shell=True) # shell=True
    log.info(f"Merged R2 files into: {merged_r2_path}")

    # 2. kallisto bus
    kallisto_out_bus = os.path.join(sc_output_dir, "output.bus")
    cmd_kb = f"kallisto bus -i {kallisto_index_path} -o {sc_output_dir} -x 10xv2 -t {threads} {merged_r1_path} {merged_r2_path}"
    run_command(cmd_kb)

    # 3. bustools sort
    sorted_bus = os.path.join(sc_output_dir, "sorted.bus")
    cmd_sort = f"bustools sort -t {threads} -o {sorted_bus} {kallisto_out_bus}"
    run_command(cmd_sort)

    # 4. Prepare corrected barcodes (remove '-1' suffix)
    corrected_barcodes_tsv = os.path.join(out_dir, "corrected_barcodes.tsv")
    log.info(f"Preparing corrected barcodes file: {corrected_barcodes_tsv}")
    try:
        # Use Python for robust handling
        count = 0
        with gzip.open(barcodes_10x_path, 'rt', encoding='utf-8') as infile, open(corrected_barcodes_tsv, 'w') as outfile:
             for line in infile:
                 bc = line.strip().split('-')[0] # Remove suffix
                 if bc: # Avoid empty lines
                     outfile.write(bc + '\n')
                     count += 1
        log.info(f"Wrote {count} corrected barcodes.")
        # Original script used gunzip | sed. This is safer.
        # cmd_bc = f"gunzip -c {barcodes_10x_path} | sed 's/-[0-9]$//g' > {corrected_barcodes_tsv}"
        # run_command(cmd_bc, shell=True) # shell=True needed for pipe
    except Exception as e:
        log.error(f"Failed to process barcode file {barcodes_10x_path}: {e}")
        sys.exit(1)


    # 5. Create t2g file (Transcript to Gene)
    # Assumes cds.fasta headers are >AlleleID Description... or just >AlleleID
    # We need AlleleID Gene (where Gene might be derived from AlleleID, e.g., HLA-A*01:01:01:01 -> HLA-A)
    # Or maybe bustools just needs TranscriptID TranscriptID (if --genecounts is used?)
    # The original script used tmpallele.txt (>ID Description) -> t2g.txt (ID Description)
    # Let's replicate that structure for bustools count -g
    t2g_file = os.path.join(sc_output_dir, "t2g.txt")
    tmp_allele_file = os.path.join(out_dir, "tmpallele.txt") # Created by build_hla_reference
    log.info(f"Creating transcript-to-gene map (t2g): {t2g_file}")
    if not os.path.exists(tmp_allele_file):
        log.error(f"tmpallele.txt not found in {out_dir}. It should be created by build_hla_reference.")
        sys.exit(1)
    try:
        with open(tmp_allele_file, 'r') as infile, open(t2g_file, 'w') as outfile:
            for line in infile:
                if line.startswith('>'):
                    # Format: ID Gene (often ID == Gene for this type of analysis)
                    transcript_id = line[1:].strip().split(' ')[0][0:]
                    # Use transcript ID as gene ID if no description, or per specific requirement
                    gene_id = line[1:].strip().split(' ')[1][0:] # Modify if gene grouping is needed
                    gene_id = gene_id.split(":")[0][0:]+":"+gene_id.split(":")[1][0:] # 2 field HLA-genes
                    outfile.write(f"{transcript_id} {gene_id}\n") # Tab separated usually preferred
    except Exception as e:
        log.error(f"Failed to create t2g file: {e}")
        sys.exit(1)


    # 6. bustools capture (using corrected barcodes)
    captured_bus = os.path.join(sc_output_dir, "filtered.bus")
    cmd_capture = f"bustools capture -b -o {captured_bus} -c {corrected_barcodes_tsv} {sorted_bus}"
    run_command(cmd_capture)

    # 7. bustools count
    count_out_dir = os.path.join(sc_output_dir, "filter_counts")
    # Ensure count output dir exists AFTER capture success
    os.makedirs(count_out_dir, exist_ok=True)
    # We created t2g, so we use -g
    cmd_count = (f"bustools count -o {count_out_dir}/output -g {t2g_file} "
                 f"-e {os.path.join(sc_output_dir, 'matrix.ec')} "
                 f"-t {os.path.join(sc_output_dir, 'transcripts.txt')} "
                 f"-m --genecounts {captured_bus}") # Using --genecounts based on t2g
    run_command(cmd_count)

    log.info(f"Counting finished. Output matrices in: {count_out_dir}")
    return count_out_dir