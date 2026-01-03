#!/bin/bash

BAM="/Disk3/LabData/scRNAseq/ALZ_PBMC/extract_file/40samples/40samples/outs/possorted_genome_bam.bam"
Barcodes="/Disk3/LabData/scRNAseq/ALZ_PBMC/extract_file/40samples/40samples/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
HLA_list="/Disk3/LabData/scRNAseq/ALZ_PBMC/extract_file/40samples/HLA-genotypes.csv"
HLADB="/Disk3/LabData/scRNAseq/Source/HLA-demultiplex/IMGTHLA-Latest/hla_nuc.fasta"
out_dir="/Disk3/LabData/scRNAseq/HLA_Splitter¡ªtest/"


BAM="/Disk3/LabData/scRNAseq/Urine_PBMC/mRNA/run_count_GEM/outs/possorted_genome_bam.bam"
Barcodes="/Disk3/LabData/scRNAseq/Urine_PBMC/mRNA/run_count_GEM/outs/filtered_feature_bc_matrix/barcodes.tsv.gz"
HLA_list="/Disk3/LabData/scRNAseq/Urine-2/HLA-genotypes_all.csv"
out_dir="/Disk3/LabData/scRNAseq/HLA_Splitter¡ªtest/"


mkdir "/Disk3/LabData/scRNAseq/HLA_Splitter¡ªtest/"

HLA-Splitter -B $BAM -b $Barcodes \
-L $HLA_list -H $HLADB \
-o $out_dir -t 10
