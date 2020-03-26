#!/bin/env python3
import os, time, json, sys
import subprocess
import matplotlib.pyplot as plt

#cmd = './main_release sparse-query /mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.100k.vcf.vcfc.sparse 22 16050075 19757157'
#cmd = './tabix /mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 22 16050075 19757157'

from .command import run_vcfc_sparse_query
from .config import Config

def measure():
    vcfc_cmd = '/home/krferrit/vcf-compression/main_release'
    tabix_cmd = '/home/krferrit/vcf-compression/tabix'
    bgzip_cmd = 'home/krferrit/vcf-compression/bgzip'

    data = {}
    sparse_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.vcfc.sparse'
    bgzip_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'

    config = Config(tabix_cmd, vcfc_cmd, bgzip_cmd)

    reference_name = '22'
    min_pos = 16050075
    max_pos = 19757157
    max_pos = 51244237
    step = 100000

    # measure percentage of search space for each, not whole file
    iter_max = int((max_pos - min_pos) / step * 0.05)

    # VCFC measurements
    # Iterate through ranges starting at front of file
    sparse_values_from_front = []
    iter = 0
    for end_pos in range(min_pos, max_pos, step):
        if iter >= iter_max:
            break
        iter += 1
        # cmd = '{vcfc_cmd} sparse-query {filename} {ref}:{start}-{end}'.format(
        #     vcfc_cmd=vcfc_cmd, filename=sparse_filename,
        #     ref=reference_name, start=min_pos, end=end_pos)
        duration = run_vcfc_sparse_query(config, sparse_filename, reference_name, min_pos, end_pos)
        sparse_values_from_front.append((min_pos, end_pos, duration))
    data['sparse_values_from_front'] = sparse_values_from_front

    # Iterate through ranges starting at back of file
    sparse_values_from_back = []
    iter = 0
    for start_pos in range(max_pos, min_pos, -step):
        if iter >= iter_max:
            break
        iter += 1
        # cmd = '/home/krferrit/vcf-compression/main_release sparse-query %s %s %d %d' % (
        #     sparse_filename, reference_name, min_pos, end_pos)
        cmd = '{vcfc_cmd} sparse-query {filename} {ref}:{start}-{end}'.format(
            vcfc_cmd=vcfc_cmd, filename=sparse_filename,
            ref=reference_name, start=start_pos, end=max_pos)
        print(cmd)
        start_time = time.time()
        proc = subprocess.Popen(cmd.split(), stdin=None, stdout=subprocess.DEVNULL)
        proc.wait()
        duration = round(time.time() - start_time, 6)
        sparse_values_from_back.append((start_pos, max_pos, duration))
    data['sparse_values_from_back'] = sparse_values_from_back

    # Iterate through ranges starting at mid of file
    sparse_values_from_mid = []
    mid_pos = int((min_pos + max_pos + 1)/2)
    iter = 0
    for start_pos in range(mid_pos, 0, -int(step/2)):
        if iter >= iter_max:
            break
        iter += 1
        end_pos = (mid_pos-start_pos) + mid_pos
        cmd = '{vcfc_cmd} sparse-query {filename} {ref}:{start}-{end}'.format(
            vcfc_cmd=vcfc_cmd, filename=sparse_filename,
            ref=reference_name, start=start_pos, end=end_pos)
        print(cmd)
        start_time = time.time()
        proc = subprocess.Popen(cmd.split(), stdin=None, stdout=subprocess.DEVNULL)
        proc.wait()
        duration = round(time.time() - start_time, 6)
        sparse_values_from_mid.append((start_pos, end_pos, duration))
    data['sparse_values_from_mid'] = sparse_values_from_mid


    # Tabix+BGZIP measurements
    tabix_bgzip_values_from_front = []
    iter = 0
    for end_pos in range(min_pos, max_pos, step):
        if iter >= iter_max:
            break
        iter += 1
        cmd = '{tabix_cmd} {filename} {ref}:{start}-{end}'.format(
            tabix_cmd=tabix_cmd, filename=bgzip_filename,
            ref=reference_name, start=min_pos, end=end_pos)
        print(cmd)
        start_time = time.time()
        proc = subprocess.Popen(cmd.split(), stdin=None, stdout=subprocess.DEVNULL)
        proc.wait()
        duration = round(time.time() - start_time, 6)
        tabix_bgzip_values_from_front.append((min_pos, end_pos, duration))
    data['tabix_bgzip_values_from_front'] = tabix_bgzip_values_from_front

    tabix_bgzip_values_from_back = []
    iter = 0
    for start_pos in range(max_pos, min_pos, -step):
        if iter >= iter_max:
            break
        iter += 1
        # cmd = '/home/krferrit/vcf-compression/main_release sparse-query %s %s %d %d' % (
        #     sparse_filename, reference_name, min_pos, end_pos)
        cmd = '{tabix_cmd} {filename} {ref}:{start}-{end}'.format(
            tabix_cmd=tabix_cmd, filename=bgzip_filename,
            ref=reference_name, start=start_pos, end=max_pos)
        print(cmd)
        start_time = time.time()
        proc = subprocess.Popen(cmd.split(), stdin=None, stdout=subprocess.DEVNULL)
        proc.wait()
        duration = round(time.time() - start_time, 6)
        tabix_bgzip_values_from_back.append((start_pos, max_pos, duration))
    data['tabix_bgzip_values_from_back'] = tabix_bgzip_values_from_back

    # Iterate through ranges starting at mid of file
    tabix_bgzip_values_from_mid = []
    mid_pos = int((min_pos + max_pos + 1)/2)
    iter = 0
    for start_pos in range(mid_pos, 0, -int(step/2)):
        if iter >= iter_max:
            break
        iter += 1
        end_pos = (mid_pos-start_pos) + mid_pos
        cmd = '{tabix_cmd} {filename} {ref}:{start}-{end}'.format(
            tabix_cmd=tabix_cmd, filename=bgzip_filename,
            ref=reference_name, start=start_pos, end=mid_pos)
        print(cmd)
        start_time = time.time()
        proc = subprocess.Popen(cmd.split(), stdin=None, stdout=subprocess.DEVNULL)
        proc.wait()
        duration = round(time.time() - start_time, 6)
        tabix_bgzip_values_from_mid.append((start_pos, end_pos, duration))
    data['tabix_bgzip_values_from_mid'] = tabix_bgzip_values_from_mid


    return data

def graph(measurements):
    if True:
        sparse_values_from_front_x = []
        sparse_values_from_front_y = []
        for val in measurements['sparse_values_from_front']:
            sparse_values_from_front_x.append(val[1] - val[0])    # range size
            sparse_values_from_front_y.append(val[2])             # duration

        tabix_bgzip_values_from_front_x = []
        tabix_bgzip_values_from_front_y = []
        for val in measurements['tabix_bgzip_values_from_front']:
            tabix_bgzip_values_from_front_x.append(val[1] - val[0])
            tabix_bgzip_values_from_front_y.append(val[2])

        sparse_values_from_back_x = []
        sparse_values_from_back_y = []
        for val in measurements['sparse_values_from_back']:
            sparse_values_from_back_x.append(val[1] - val[0])
            sparse_values_from_back_y.append(val[2])

        tabix_bgzip_values_from_back_x = []
        tabix_bgzip_values_from_back_y = []
        for val in measurements['tabix_bgzip_values_from_back']:
            tabix_bgzip_values_from_back_x.append(val[1] - val[0])
            tabix_bgzip_values_from_back_y.append(val[2])

        sparse_values_from_mid_x = []
        sparse_values_from_mid_y = []
        for val in measurements['sparse_values_from_mid']:
            sparse_values_from_mid_x.append(val[1] - val[0])
            sparse_values_from_mid_y.append(val[2])

        tabix_bgzip_values_from_mid_x = []
        tabix_bgzip_values_from_mid_y = []
        for val in measurements['tabix_bgzip_values_from_mid']:
            tabix_bgzip_values_from_mid_x.append(val[1] - val[0])
            tabix_bgzip_values_from_mid_y.append(val[2])

    plt.xlabel('Query Range')
    plt.ylabel('Time (seconds)')
    plt.plot(sparse_values_from_front_x, sparse_values_from_front_y, label='VCFC Sparse Offset Index (from start of file)')
    plt.plot(sparse_values_from_back_x, sparse_values_from_back_y, label='VCFC Sparse Offset Index (from end of file)')
    plt.plot(sparse_values_from_mid_x, sparse_values_from_mid_y, label='VCFC Sparse Offset Index (from mid of file)')

    plt.plot(tabix_bgzip_values_from_front_x, tabix_bgzip_values_from_front_y, label='BGZIP+Tabix Index (from start of file)')
    plt.plot(tabix_bgzip_values_from_back_x, tabix_bgzip_values_from_back_y, label='BGZIP+Tabix Index (from end of file)')
    plt.plot(tabix_bgzip_values_from_mid_x, tabix_bgzip_values_from_mid_y, label='BGZIP+Tabix Index (from mid of file)')

    plt.legend()
    plt.show()

def main(args):
    if args[0] == 'measure':
        measurements = measure()
        with open('range-measurements.json', 'w') as f:
            json.dump(measurements, f)
    elif args[0] == 'graph':
        with open('range-measurements.json', 'r') as f:
            measurements = json.load(f)
        graph(measurements)


if __name__ == '__main__':
    main(sys.argv[1:])
