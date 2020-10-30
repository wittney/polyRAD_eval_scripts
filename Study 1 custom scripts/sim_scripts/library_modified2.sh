#!/bin/bash

for file in C8B56*
do
cat $file  | sed -e 's/msp.*$//g' -e '/^@.*:/s/:/:56:C8B56ANXX1:1:3:1101:1203:20371:N:0/g' | gzip f{$file}.fq.gz
