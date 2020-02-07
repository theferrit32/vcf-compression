#!/bin/bash
set -e -x
node_addr='c29'
dir="vcf-compression"

rsync -av -e 'ssh -J "krferrit@arc.csc.ncsu.edu"' ../${dir}/  krferrit@${node_addr}:${dir}/
