#!/bin/bash

for file in fC8B56*
do
    mv -i "${file}" "${file/3/1}"
done
