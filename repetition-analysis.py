#!/bin/env python3
import re, os
import matplotlib.pyplot as plt
font = {'family' : 'normal',
        #'weight' : 'bold',
        'size'   : 18}
plt.rc('font', **font)

# vcf_file = '/home/me/dev/stanford/1000genomes/' \
#    + 'ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'

################
# parameters
it = 100000
original_filename = 'chr16.%d.vcf' % it
filename = '16.%d.out' % it
allele_gt = '1|1'
seq_len = {
    '0|0': 128,
    '1|1': 32,
    '0|1': 32,
    '1|0': 32
}.get(allele_gt, 32)
# end parameters
################

with open(filename, 'r') as f:
    lines = f.readlines()

p = re.compile('sample %s repeated ([0-9]+) times' % (allele_gt.replace('|', '\\|')))
d = {}
for line in lines:
    m = p.match(line)
    if m:
        num = m.group(1)
        if num not in d:
            d[num] = 0
        d[num] = d[num] + 1
print(d)

x = sorted(list(d.keys()))
x = sorted([int(dk) for dk in d.keys()])
y = [d[str(xk)] for xk in x]
print(x)
print(y)

def calc_savings(d):
    savings = 0
    for k,v in d.items():
        b_orig = 4 # 3 bytes to save sequence '0|0' followed by tab
        b_new = 1
        savings += int(v) * int(k) * (b_orig - b_new)
    return savings

def calc_ratio(d, original_filename):
    savings = float(calc_savings(d))
    original_filesize = os.path.getsize(original_filename)
    print('savings: %f, original_filesize: %f' % (savings, original_filesize))
    return (original_filesize - savings) / original_filesize

print('savings: %s' % calc_savings(d))
reduction_ratio = 1.0 - calc_ratio(d, original_filename)
print('reduction ratio: %s' % reduction_ratio)

plt.bar(x, y)
plt.xticks([xt for xt in range(min(x), max(x)+1, 4)])
plt.xlabel('Repetition count (max sequence length %d)' % seq_len)
plt.ylabel('Instances of repetition sequence')
plt.title('Run lengths of %s allele alts of 100k samples on CHR16 (reduction: %.8f)' % (allele_gt, reduction_ratio))
plt.show()
