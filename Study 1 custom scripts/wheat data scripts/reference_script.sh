#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 1
#SBATCH --mem=32g
#SBATCH -N 1
# ----------------Load Modules--------------------
module load SAMtools/0.1.20-IGB-gcc-4.9.4
# ----------------Commands------------------------

gunzip /home/n-z/wittney2/Reference_genome_wheat/Taestivum_296_v2.fa.gz
samtools faidx /home/n-z/wittney2/Reference_genome_wheat/Taestivum_296_v2.fa
gzip /home/n-z/wittney2/Reference_genome_wheat/Taestivum_296_v2.fa
cat /home/n-z/wittney2/Reference_genome_wheat/Taestivum_296_v2.fa.fai | cut -f 1  > ./chrom.list
