#!/bin/bash
#gff2bed   < ../TriTrypDB-56_TbruceiTREU927.gff  | grep -P "\tCDS\t" > CDS.bed
#gff2bed   < ../TriTrypDB-56_TbruceiTREU927.gff  | grep -P "\tmRNA\t" > mRNA.bed
#perl  perl gff2bed.pl

THREADS=28
INPUT_FILE=/media/logen/sdc/trypanosoma/CPSF/long_reads/Customer_BMK191104-W104-0101/BMK191104-W104-0101-cleandata/W104-01-N04.fq.gz
OUT_HEAD=bmk_cpsf30

source /home/logen/anaconda3/etc/profile.d/conda.sh
conda activate base

porechop -i $INPUT_FILE  -o $OUT_HEAD.fq.gz --threads $THREADS
minimap2 -ax map-ont -Y -k14 -t $THREADS --sam-hit-only --secondary=no ../TriTrypDB-56_TbruceiTREU927_Genome.fasta $OUT_HEAD.fq.gz > $OUT_HEAD.sam
perl match_pA.pl $OUT_HEAD.sam

samtools view -@ $THREADS -hu ${OUT_HEAD}_perl.sam | samtools sort -@ $THREADS -o ${OUT_HEAD}_perl.bam - 
samtools index ${OUT_HEAD}_perl.bam
bamToBed -i ${OUT_HEAD}_perl.bam > ${OUT_HEAD}.bed

OUT_HEAD=bmk_cpsf30RNAi
bedtools intersect -s -a $OUT_HEAD.bed -b mRNAcds.bed -wa -wb > ${OUT_HEAD}_inter.bed
