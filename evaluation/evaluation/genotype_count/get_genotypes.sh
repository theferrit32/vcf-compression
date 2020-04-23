#!/bin/bash
set -xe

directory=/mnt/beegfs/krferrit/1000Genomes-vcf-orig
outfile=/mnt/ext4/genotypes-new.gz
count_outfile=/mnt/ext4/genotype-counts.txt

if [[ -f $outfile ]]; then
  echo "$outfile already exists, remove?"
  read -n 1 -p "Press any key to remove and continue"
  rm $outfile
fi

# Use compressed versions of files, and just open a gzip decompress
# stream on them for reading. Better for I/O performance.
for f in $directory/*.vcf.gz; do
  echo "Scanning $f"
  # Exclude header/meta lines, exclude 9 variant columns
  # Replace tabs with newlines
  gzip -d -c $f | grep -vP '^#.*' | cut -f 10- | tr '\t' '\n' | gzip -c - >> $outfile
done

echo "Done recording genotypes in to $outfile"

cat $outfile | tr '\t' '\n' | sort -S 10G --parallel=16 | uniq -c > $count_outfile

echo "Done counting"

