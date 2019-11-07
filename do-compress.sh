#!/bin/bash

fname="test-16-100.vcf"
comp_fname="${fname}.vcfc"
decomp_fname="${fname}.decompressed"

./main compress $fname $comp_fname 2>&1 | tee compress.log

hexdump $comp_fname -C | tee "$comp_fname.hexdump"

./main decompress $comp_fname $decomp_fname 2>&1 | tee decompress.log


