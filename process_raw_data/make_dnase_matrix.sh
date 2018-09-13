#!/bin/bash

step="last"
#module load bedtools
#output_path = "/srv/scratch/wychen66/code_output/dnase_output/2018-3-22/"

if [ "${step}" = "presort" ]
then
    python make_dnase_matrix_pre_shell_script.py /srv/scratch/wychen66/code_output/dnase_output/2018-3-22/
    while read line
    do
    echo $line
    zcat ${line} | cut -f1-3 >> ${output_path}presort.bed  
    done < idr_files.bed 
fi

if [ "${step}" = "sort" ]
then
    sort -k1,1 -k2,2n presort.bed > postsort.bed 
    #rm presort.bed
fi

if [ "${step}" = "merge" ]
then
    bedtools merge -i postsort.bed > idr_files_merged.bed 
    #rm postsort.bed
fi

if [ "${step}" = "eval" ]
then
    zcat -f idr_files_merged.bed | awk '{print ($3-$2)}' | awk '{sum+=$1} END {print sum}'
fi

if [ "${step}" = "add_4th_col" ]
then
    awk '{print $1":"$2"-"$3}' idr_files_merged.bed > idr_files_merged_4th.bed
    paste idr_files_merged.bed idr_files_merged_4th.bed > idr_files_merged_4.bed 
    #rm idr_files_merged.bed
    #rm idr_files_merged_4th.bed
fi

if [ "${step}" = "last" ]
then
    touch dnase_file.txt

    while read line; do
    echo $line
    zcat -f ${line} | cut -f1-3 | bedtools genomecov -i stdin -g /mnt/data/annotations/by_organism/human/hg19.GRCh37/hg19.chrom.sizes -bg > unsorted.bedGraph
    sort -k1,1 -k2,2n unsorted.bedGraph > sorted.bedGraph
    /srv/scratch/wychen66/software/bedGraphToBigWig sorted.bedGraph /mnt/data/annotations/by_organism/human/hg19.GRCh37/hg19.chrom.sizes sorted.bw
    /srv/scratch/wychen66/software/bigWigAverageOverBed sorted.bw idr_files_merged_4.bed dnase_out.tab 

    cut -f6 dnase_out.tab > col.txt
    paste dnase_file.txt col.txt > dnase_final.txt
    mv dnase_final.txt dnase_file.txt
    done < tag_files.bed # 122 lines

    mv dnase_file.txt final_dnase_matrix.txt 
#    rm unsorted.bedGraph
#    rm sorted.bedGraph
#    rm sorted.bw
#    rm dnase_out.tab
#    rm col.txt
fi

if [ "${step}" = "label" ]
then
    while read line; do
    echo $line >> /srv/scratch/wychen66/output/dnase_output/labels.txt
    done < /srv/scratch/wychen66/output/dnase_output/tag_files.bed # 122 lines

    python make_dnase_matrix_post_shell_script.py /srv/scratch/wychen66/output/dnase_output/
fi




