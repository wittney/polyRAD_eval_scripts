#!/bin/bash

cat f_out_fastq_3.fq | awk '{if(NR%4==2) {printf("CAGTCC%s\n",$0)} else if(NR%4 == 0) {printf("IIIIIII%s\n", $0)} else {print $0}}' | gzip > C8B56ANXX4_3_.modified.fq.gz
zcat C8B56ANXX4_3_.modified.fq.gz | sed -e 's/msp.*$//g' -e '/^@.*:/s/:/:56:C8B56ANXX1:1:3:1101:1203:20371:N:0/g' > fC8B56ANXX4_3_.modified.fq.gz
