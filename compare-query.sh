#!/bin/bash
set -e
set -x
flush_cache() {
	echo 3 | sudo tee /proc/sys/vm/drop_caches
}

vcf="/mnt/ext/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf"
BCR_COOR_TABIX="22:23521891-23660224"
BCR_COOR_VCFC="22 23521891 23660224"

flush_cache

time htslib/tabix "${vcf}.gz" ${BCR_COOR_TABIX} > BCR-tabix.vcf

flush_cache

time vcf-compression/main_release query "${vcf}.vcfc" ${BCR_COOR_VCFC} > BCR-vcfc.vcf
