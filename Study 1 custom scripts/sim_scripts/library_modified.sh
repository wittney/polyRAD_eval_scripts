#!/bin/bash


zcat C8B56ANXX1_3_.modified.fq.gz | sed -e 's/msp.*$//g' -e '/^@.*:/s/:/:56:C8B56ANXX1:1:3:1101:1203:20371:N:0/g' > fC8B56ANXX1_3_.modified.fq.gz
