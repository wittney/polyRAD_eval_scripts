#!/bin/bash

for file in msp*
do
zcat $file |  awk 'BEGIN {RS = ">" ; FS = "\n"} NR > 1 {print "@"$1"\n"$2"\n+"; for(c=0;c<length($2);c++) printf "H"; printf "\n"}' | gzip > /media/HDD/wittney2/Simulated_sequence_reads/rad_reads/fastq_files/${file}.fq.gz 
done
