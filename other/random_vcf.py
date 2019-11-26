#!/bin/env python3
import random, math

sample_count = 10
variant_count = 50
alt_count = 2 # valid values: 1, 2, 3
alt_indexes = range(0, alt_count+1)
random.seed(5)
output_filename = 'test-%d-%d.vcf' % (sample_count, variant_count)
bases = ['A', 'T', 'G', 'C']
pos_counter = 10000
tsv_header = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']

def cumulative_dist(seq):
    c = 0
    dist = []
    for s in seq:
        c += s
        dist.append(c)
    return dist

def choice_with_probability(seq, probs):
    assert(len(seq) == len(probs))
    s = sum(probs)
    cdist = cumulative_dist(probs)
    r = random.random() * s
    print('s: %s, r: %.4f' % (s, r))
    for i in range(0, len(seq)):
        if r < cdist[i]:
            return seq[i]

def fwrite(f, s):
    f.write(s.encode('utf-8'))

with open(output_filename, 'wb') as f:
    # write headers
    fwrite(f, '##fileformat=VCFv4.1\n')
    fwrite(f, '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    fwrite(f, '##fileDate=20150218\n')
    fwrite(f, '#' + '\t'.join(tsv_header))
    digits = int(math.ceil(math.log10(sample_count)))
    sample_name_fmt = 'HG%0' + str(digits) + 'd'
    for j in range(0, sample_count):
        fwrite(f, '\t')
        fwrite(f, sample_name_fmt % j)
    fwrite(f, '\n')


    # write variants
    for i in range(0, variant_count):
        chrom = '1'
        pos = str(pos_counter)
        pos_counter += 2
        identifier = 'var' + str(i)
        ref = random.choice(bases)
        alts = [b for b in bases if b != ref]
        random.shuffle(alts)
        alts = alts[:alt_count]
        qual = '100'
        filter_string = 'PASS'
        info = 'INFO'
        fmt = 'GT'
        line_terms = [chrom, pos, identifier, ref, ','.join(alts),
                qual, filter_string, info, fmt]
        for j in range(0, sample_count):
            alt_vals = [0, 1, 2]
            alt_probs = [0.90, 0.08, 0.02]
            allele1 = choice_with_probability(alt_vals, alt_probs)
            allele2 = choice_with_probability(alt_vals, alt_probs)
            line_terms.append(str(allele1) + '|' + str(allele2))

        fwrite(f, '\t'.join(line_terms))
        fwrite(f, '\n')

    print('finished writing %s' % output_filename)
