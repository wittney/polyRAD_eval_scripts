#!/bin/bash
#SBATCH -p normal
#SBATCH --mem=100g
#SBATCH -N 2
#SBATCH -n 4
#SBATCH -D /home/n-z/wittney2
module load radinitio/1.1.1-IGB-gcc-4.9.4-Python-3.6.1

radinitio  -g /home/n-z/wittney2/Reference_genome_wheat/Taestivum_296_v2_m.fa  -l /home/n-z/wittney2/new_chrom_file2 -o /home/n-z/wittney2/Simulated_seq_reads/ -e PstI -d MspI -b ddRAD -c 9 -p 1 -s 400 -n 2500  
