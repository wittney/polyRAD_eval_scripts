#!/bin/bash
# ----------------SLURM Parameters----------------
#SBATCH -p normal
#SBATCH -n 4
#SBATCH --mem=64G
#SBATCH -N 2
#SBATCH -D /home/n-z/wittney2
# ----------------Load Modules--------------------
module load TASSEL/5.2.28-Java-1.8.0_121
# ----------------Commands------------------------
run_pipeline.pl -Xmx60G -fork1 -GBSSeqToTagDBPlugin -e PstI -i /home/n-z/wittney2/Wheat_Data/wheat_fastq_files/ -db /home/n-z/wittney2/Wheat_Data/wheat_tassel_output/wheat_database_PstI_11222019.db -k /home/n-z/wittney2/Wheat_Data/wheat_PstI_key_All.txt -kmerLength 64 -minKmerL 20 -mxKmerNum 100000000 -endPlugin -runfork1 
