#!/bin/env python3
import os
import math
import matplotlib.pyplot as plt
#fname = '/home/me/dev/stanford/1000genomes/' +\
#    'ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf'

def load_file(position_file_path):
    positions = []
    uncompressed_lengths = []
    compressed_lengths = []

    with open(position_file_path) as fin:
        for line in fin:
            line = line.strip()
            terms = line.split(' ')
            #print(terms)
            positions.append(int(terms[0]))
            uncompressed_lengths.append(int(terms[1]))
            compressed_lengths.append(int(terms[2]))

    return positions, uncompressed_lengths, compressed_lengths


def main():
    fname = 'start-positions.txt'
    positions, uncompressed_lengths, compressed_lengths = load_file(fname)
    print('Finished reading file ' + fname)

    x = []
    y = []
    items = sorted(zip(positions, compressed_lengths), key=lambda a: a[0])
    prev_pos = 0
    d = {}
    gap_lengths = []
    pos_gaps = []
    block_counts = []
    internal_fragment_lengths = []
    block_length = 4096
    for pos, compressed_length in items:
        pos_gap = (pos - prev_pos)
        pos_gaps.append(pos_gap)
        gap_length = (pos - prev_pos) - compressed_length
        if gap_length not in d:
            d[gap_length] = 0
        d[gap_length] += 1
        gap_lengths.append(gap_length)

        # compute block alignment
        if compressed_length <= block_length:
            blocks = 1
            internal_fragment_length = block_length - compressed_length
        else:
            blocks = int(math.ceil(compressed_length / block_length))
            internal_fragment_length = block_length - (compressed_length % block_length)
        block_counts.append(blocks)
        internal_fragment_lengths.append(internal_fragment_length)


    # for gap_length, count in sorted(d.items(), key=lambda a: a[0]):
    #     x.append(gap_length)
    #     y.append(count)

    fig, axes = plt.subplots(1, 3)

    #axes[0].hist(gap_lengths, bins=100)
    axes[0].hist(compressed_lengths, bins=8)
    axes[1].hist(internal_fragment_lengths, bins=100)
    axes[2].hist(block_counts, bins=100)
    axes[0].set_xlabel('Compressed line lengths')
    axes[1].set_xlabel('Sparse Internal FS Block Fragmentation')
    axes[2].set_xlabel('%d byte Blocks Required Per Line' % block_length)

    axes[0].set_ylabel('Occurences')
    axes[1].set_ylabel('Occurences')
    axes[2].set_ylabel('Occurences')

    axes[0].set_yscale('log')
    axes[1].set_yscale('log')
    axes[2].set_yscale('log')

    #total_internal_fragmentation = sum([block_length - fl for fl in internal_fragment_lengths])
    total_internal_fragmentation = sum(internal_fragment_lengths)
    print('Total internal fragmentation: ' + str(total_internal_fragmentation))
    print('Mininum position gap: %d' % min(pos_gaps))
    print('Maximum compressed line block length: %d' % max(block_counts))

    print(sorted(compressed_lengths, reverse=True)[:50])
    plt.show()
    #plt.hist(gap_lengths, bins=100)
    #plt.yscale('log')
    # plt.xlabel('Compressed line lengths')
    # plt.ylabel('Number of occurences')
    # plt.title('Number of gap length occurrences (1000 Genomes Phase 3 CHR 22 VCF)')
    # plt.show()


if __name__ == '__main__':
    main()

