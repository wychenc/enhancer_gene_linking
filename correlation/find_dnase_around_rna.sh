#!/bin/bash
bedtools window -a /srv/scratch/wychen66/output/rna_output/rna_counts_gene-by-cellline_unsorted_row_labels_tss.txt -b /srv/scratch/wychen66/output/dnase_output/idr_files_merged.bed -w 1000000 > dnase_within_1MB_of_rna_tss.bed
