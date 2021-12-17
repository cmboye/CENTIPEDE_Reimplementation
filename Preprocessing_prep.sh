#!/bin/bash
#This is a script to preprocess data for centiPEDE. It includes the PWMscan, and PhastCons
#Required Input Data: PWM for motif of interest (MEME format), genome of interest (.fa), DNAse-Seq Data (.narrowPeak.bed); Optional: Conservation data (.bw)

#module load meme
#module load bedtools
#module load ucscbrowser/2019-05-25

#Set variables
PWM= 
DNASE= 
CONS=
GENOME=

zcat ${DNASE}.bed.gz | awk '{if ($7 > 8) print}' > ${DNASE}_filtered.bed #Filter for p<1e-8 DNase peaks; note that narrowPeak uses -log10(p)
bedtools getfasta -fi ${GENOME}.fa -bed ${DNASE}_filtered.bed -fo > ${DNASE}_FASTA.fa #Find sequences within the regions in the peaks
#Note: Will need to convert to MEME format if not provided
fimo --text --parse-genomic-coord ${PWM}.meme ${DNASE}_FASTA.fa > ${DNASE}_${PWM}.fimo.txt 

###PhastCons-- Conservation data from UCSC
#First we will use this publicly available script to define the bigWigRegions function that extracts bigWig data from certain regions indicated by a bed file
bigWigRegions() {
  bw="$1"
  bed="$2"
  IFS=$'\t'
  while read chrom beg end rest; do
    # Write temporary files to RAM.
    out="/dev/shm/bigWigRegions_${USER}_${$}_${chrom}_${beg}_${end}"
    bigWigToBedGraph -chrom=$chrom -start=$beg -end=$end "$bw" "$out"
  done < "$bed"
  # Print the temporary files to stdout and then delete them.
  cat /dev/shm/bigWigRegions_${USER}_${$}_* | sort -k1V -k2n -k3n
  rm -f /dev/shm/bigWigRegions_${USER}_${$}_*
}
bigWigRegions ${CONS}.bw ${DNASE}_filtered.bed > ${CONS}_${DNASE}_raw.bed
awk '{print $2 "\t" $3 "\t" $4}' ${DNASE}_${PWM}.fimo.txt | tail -n +2 > ${DNASE}_prep.bed
#Need to sort before using bedtools map
sort -k1,1 -k2,2n ${DNASE}_prep.bed -o ${DNASE}_prep.bed
sort -k1,1 -k2,2n ${CONS}_${DNASE}_raw.bed -o ${CONS}_${DNASE}_raw.bed
bedtools map -a ${DNASE}_prep.bed -b ${CONS}_${DNASE}_raw.bed -c 4 -o mean > ${CONS}_${DNASE}_mean.txt
