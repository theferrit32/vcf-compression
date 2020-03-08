#!/bin/bash
set -e -x
dir="vcf-compression"
#node_addr='c0'
# rsync -av -e 'ssh -J "krferrit@arc.csc.ncsu.edu"' ../${dir}/  krferrit@${node_addr}:${dir}/

rsync -av ../${dir}/ krferrit@arc:${dir}/
