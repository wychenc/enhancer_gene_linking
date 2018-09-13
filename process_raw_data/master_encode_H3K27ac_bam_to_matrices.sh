#!/bin/bash

# list of encode H3K27ac bam file urls
all_bam_urls_f='/oak/stanford/groups/akundaje/projects/enhancer_gene_linking/data/encode_H3K27ac_bam/all_bam_urls.txt'
tss_10KB_window_gtf_4cols_unique_f='/oak/stanford/groups/akundaje/projects/enhancer_gene_linking/master_bam_2_matrices/tss_10KB_window_gtf_4cols_unique.bed'
bigwig_from_bam_folder='/oak/stanford/groups/akundaje/projects/enhancer_gene_linking/data/encode_H3K27ac_bam_to_bw/'
bigWigAverageOverBed_c='/oak/stanford/groups/akundaje/projects/enhancer_gene_linking/software/bigWigAverageOverBed'

#mkdir correct_tss_tabs
#mkdir correct_bam_tabs

# download list of H3K27ac bam files 
input=$1
if [[ $input == '1' ]]; then
	echo $input
	wget -i $all_bam_urls_f
fi

# convert .bam to .bigwig using bamCoverage

# bigWigAverageOverBed to get .tab from bigwig
### get bed file of tss's with 10KB window from gtf file
# tss_10KB_window_gtf_4cols_unique.bed

## break $all_bam_urls_f into smaller lists
## 0_200url.txt, in creasing steps of 200, to 1400_1500url.txt
## for each txt file, do:

if [[ $input == '2' ]]; then
	echo $input
	while read line
	do
	    bw=$bigwig_from_bam_folder$(basename $line).coverage.bw
	    echo $bw
	    $bigWigAverageOverBed_c $bw $tss_10KB_window_gtf_4cols_unique_f $(basename $line).tab
	done < $all_bam_urls_f
fi


# get correct .tab 
if [[ $input == '3' ]]; then
	## move correct tss .tab into directory
	mkdir correct_tss_tabs 
	for file in *bam.tab; do
	    echo $file
	    len=$(awk 'END{print NR}' $file)
	    echo $len
	    if [ "$len" = "57730" ]; then
	         mv $file correct_tss_tabs/
	    fi
	done

	## move correct enh .tab into directory 
	mkdir correct_bam_tabs
	for file in *bam.tab; do
	    echo $file
	    len=$(awk 'END{print NR}' $file)
	    echo $len
	    if [ "$len" = "845537" ]; then
	         mv $file correct_bam_tabs/
	    fi
	done
fi


if [[ $input == '4' ]]; then

	## list correct enh tabs
	touch correct_bam_tabs/correct_bam_tab_list.txt
	for file in correct_bam_tabs/*bam.tab; do
	    echo $(basename $file) >> correct_bam_tabs/correct_bam_tab_list.txt
	done

	## make enh matrix from correct enh tabs
	touch correct_bam_tabs/bam_enh_matrix.txt
	while read line; do
	    echo $line
  	    cut -f6 $line > col.txt
	    paste bam_enh_matrix.txt col.txt > enh.txt
	    mv enh.txt correct_bam_tabs/bam_enh_matrix.txt
	done < correct_bam_tabs/correct_bam_tab_list.txt
	rm col.txt

	## list correct tss tabs
	touch correct_tss_tabs/correct_tss_tab_list.txt
	for file in correct_tss_tabs/*bam.tab; do
	    echo $(basename $file) >> correct_tss_tabs/correct_tss_tab_list.txt
	done

	## make tss matrix from correct tss tabs
	touch correct_tss_tabs/bam_tss_matrix.txt
	while read line; do
	    echo $line
	    cut -f6 $line > col.txt
	    paste bam_tss_matrix.txt col.txt > enh.txt
	    mv enh.txt correct_tss_tabs/bam_tss_matrix.txt
	done < correct_tss_tabs/correct_tss_tab_list.txt
	rm col.txt
fi


