#! usr/env/bin bash

# problem-set-3 answers
datasets='/vol3/home/pezoas/data-sets'

# Q1 Identify the size of the largest overlap between CTCF and H3K4me3 
tfbs=$datasets/bed/encode.tfbs.chr22.bed.gz
histone=$datasets/bed/encode.h3k4me3.hela.chr22.bed.gz

zcat $tfbs | awk '$4 == "CTCF"' > ctcf-peaks.bed

a_1=$(bedtools intersect -a ctcf-peaks.bed -b $histone -wo \
	 | awk '{print $NF}' \
	 | sort -nr | head -n1)

echo "answer-1: $a_1"

# Q2 Use bedtools to calculate GC content on nucleotides 19,000,000 to 19,000,500 on chr22
# bedtools nuc build bedfile with coordinates echo -e echo -e "this\that" 
chr22=$datasets/fasta/hg19.chr22.fa
echo -e 'chr22\t19000000\t19000500' > interval.bed
# prints columns with nucleotide percentages want fifth column and bottom line 
a_2=$(bedtools nuc -fi $chr22 -bed interval.bed | cut -f5 | tail -n1)
echo "answer-2: $a_2"

# Q3 identify length of CTCF chip-seq peak in that has largest mean signal 
# map interval calculate signal in intervals prints columns with intervals need to subtract 3 from 2 last column is signal period means no signal? sort signal in last column for largest then substract columns
ctcf=$datasets/bedtools/ctcf.hela.chr22.bg
a_3=$(bedtools map -a ctcf-peaks.bed -b $ctcf -c 4 -o mean | sort -k5n | tail -n1 | awk '{print $3 - $2}')  
echo "answer-3: $a_3"

# Q4 Identify the gene promoter 1000 bp up of TSS with highest median signal in ctcf.hela.
# Report gene name to find promoters use bedtools flank compare intervals to signal with map summary stat is median
tss=$datasets/bed/tss.hg19.chr22.bed
genome=$datasets/genome/hg19.genome 
a_4=$(bedtools flank -i $tss -g $genome -l 1000 -r 0 -s | bedtools sort -i \
	| bedtools map -a - -b $ctcf -c 4 -o median \
	| sort -k7n | tail -n1 | awk '{print $4}')
echo "answer-4: $a_4"

# Q5 Identify the longest interval on chr22 not covered by genes.hg19.bed use betools complement
genes=$datasets/bed/genes.hg19.bed_sorted
hg19_sorted=$datasets/genome/hg19.genome_sorted
a_5=$(bedtools complement -i $genes -g $hg19_sorted | awk '$1 == "chr22" {print $1,$2,$3-$2}' | head -n1 | awk 'BEGIN {OFS=""} {print $1,":",$2,"-",$3}')
echo "answer-5: $a_5" 


  
