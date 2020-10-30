#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 2
#SBATCH --mem=41g
#SBATCH -N 2
# ----------------Load Modules--------------------
module load Bowtie2/2.3.5.1-IGB-gcc-4.9.4

bowtie2 -k 4 --very-sensitive -U /home/n-z/wittney2/Wheat_Data/wheat_tassel_output/input_bowtie_PstI_11232019.fa -x /home/n-z/wittney2/Wheat_Data/wheat_tassel_output/wheat_genome_index/wheat_genome_index -S /home/n-z/wittney2/Wheat_Data/wheat_tassel_output/wheat_output__PstI_11232019.sam
