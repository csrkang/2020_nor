#!/bin/bash
wget -P fasta/ https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/752/655/GCA_002752655.1_ASM275265v1/GCA_002752655.1_ASM275265v1_genomic.fna.gz
gunzip fasta/GCA_002752655.1_ASM275265v1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/752/655/GCA_002752655.1_ASM275265v1/GCA_002752655.1_ASM275265v1_genomic.gtf.gz
gunzip GCA_002752655.1_ASM275265v1_genomic.gtf
bwa index fasta/GCA_002752655.1_ASM275265v1_genomic.fna
fileName='SRX_Acc_List.txt'
while IFS= read -r line
do
  prefetch ${line}
  fasterq-dump "$line" -O fastq/
  sickle se -f fastq/${line}.fastq -t sanger -o fastq/${line}_trimmed.fastq
  fastqc -o fastqc/ fastq/${line}_trimmed.fastq
  bwa mem fasta/GCA_002752655.1_ASM275265v1_genomic.fna fastq/${line}_trimmed.fastq |\
          samtools view -Sbh - > bam/${line}_aligned.bam
  rm fastq/*.fastq
  featureCounts -t gene -g gene_id -a GCA_002752655.1_ASM275265v1_genomic.gtf -o counts/${line}_counts.txt\
                 bam/${line}_aligned.bam
  rm bam/*.bam
done < "$fileName"
R CMD BATCH counts_to_FPKM.R
