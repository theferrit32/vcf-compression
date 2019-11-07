# vcf-compression

## Overview
Compression strategy for VCF files.

1. Leave VCF metadata lines (starting with '##') as-is. A few KiB could be compressed out of these, but since these savings are constant and do not scale up with larger files, they are a low priority.
2. Leave TSV header line as-is. This line occurs exactly once, begins with a '#', and each column is unique so there is little benefit to going through the effort to compress it, with similar reasoning as for the metadata lines.
3. For each VCF data line:
    - Leave the first 8 columns uncompressed. These are mandatory columns each of which are somewhat unique, and each is relatively short, besides `INFO`.  The `INFO` column can be arbitrarily long and in many cases is very long, on the order of hundreds of bytes.
    - If there is another column after the initial 8, it is the `FORMAT` column specifying the format of values in each sample genotype column. (See VCFv4.3 spec section 1.6.2).
    - For the remaining 0 or more columns, see section *VCF Genotype Compression*.

## VCF Genotype Compression

In the usual case of VCF data lines, there are 9 mandatory TSV columns (including `FORMAT`), followed by a variable number of columns, each representing genotype information about a particular DNA sample. The number of these sample columns is pre-computable by looking at the number of sample columns in the TSV header row, as each row must have the same number of columns.

The variant call process is to compare individual samples with a reference genome. By construction, this reference genome seeks to reflect the most common