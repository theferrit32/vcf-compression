#!/bin/bash

fname="test-16-100.vcf"

./main compress $fname "$fname.vcfc" 2>&1 | tee compress.log

hexdump "$fname.vcfc" -C | tee "$fname.vcfc.hexdump"

./main decompress "$fname.vcfc" "$fname.decompressed" 2>&1 | tee decompress.log
