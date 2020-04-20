#!/bin/env python3
import os, time, json, sys
import re
import copy
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import humanize

from evaluation.command import (
    run_tabix,
    run_vcfc_sparse_query,
    create_binned_index,
    run_vcfc_binned_index_query,
    run_vcfc_binned_index_timing_profile
)
from evaluation.config import Config

plt.rcParams.update({
    'figure.dpi': 200,
    'savefig.dpi': 200,
    'figure.figsize': (10.0, 5.0),
    'axes.labelsize': 'medium',
    'font.size': 10
})

colors = [
    # Primary data color, secondary data feature color

    # ('blue', 'aqua'),
    ('blue', 'dodgerblue'),
    ('orange', 'coral'),
    # ('darkgreen', 'lightgreen'),
    ('darkgreen', 'limegreen'),
    ('purple', 'orchid'),
    ('red', 'tomato'),
]

vcfc_dir = '/home/krferrit/vcf-compression/'
tabix_cmd = '/home/krferrit/vcf-compression/tabix'
bgzip_cmd = '/home/krferrit/vcf-compression/bgzip'

# Chromosome 22
vcfc_filename = '/mnt/ext4/backup/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.vcfc'
sparse_filename = vcfc_filename + '.sparse'
bgzip_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
bcf_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf'
reference_name = '22'
min_pos = 16050075
max_pos = 51244237

# Chromosome 1
# vcfc_filename = '/mnt/ext4/backup/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.vcfc'
# sparse_filename = vcfc_filename + '.sparse'
# bgzip_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
# bcf_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf'
# reference_name = '1'
# min_pos = 16050075
# max_pos = 51244237



test_runs = 10

def measure_binned_index_creation_time():
    data = {}
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    bin_sizes = [
        *range(5, 50, 5),
        *range(50, 200, 25),
        *range(200, 2000+1, 100)
    ]

    creation_times = []
    for bin_size in bin_sizes:
        times = []
        for _ in range(test_runs):
            print('Creating binned index with bin size %d' % bin_size)
            bin_index_creation_time = create_binned_index(config, vcfc_filename, bin_size)
            print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))
            times.append(bin_index_creation_time)
        creation_times.append({
            'bin_size': bin_size,
            'time': sum(times)/len(times),
            'stddev': np.std(times)
        })

    data['vcfc_binned_index_creation_time'] = {
        'data': creation_times,
        'label': 'VCFC Binned Index Creation Time'
    }

    return {
        'data': data,
        'title': 'VCFC Binned Index Creation Time by Bin Size',
        'name': 'binned-index-creation-time',
        'xlabel': 'Bin Size',
        'ylabel': 'Time (seconds)'
    }

def graph_binned_index_creation_time(measurements):
    measurement = measurements['data']['vcfc_binned_index_creation_time']
    xvals = [m['bin_size'] for m in measurement]
    yvals = [m['time'] for m in measurement]
    errs = [m['stddev'] for m in measurement]

    plt.errorbar(xvals, yvals, errs)

    plt.xlabel(measurements['xlabel'])
    plt.ylabel(measurements['ylabel'])
    plt.title(measurements['title'])
    plt.legend()

    # Tighten layout
    plt.margins(0.01)
    plt.tight_layout()

    if 'name' in measurements:
        filename = measurements['name'] + '.png'
    else:
        filename = measurements['title'].replace(' ', '_') + '.png'
    plt.savefig(filename, dpi=200)
    print('Saved figure: %s' % filename)



def measure_binned_index_single_exhaustive(queries:int=200):
    data = {}

    step = int((max_pos - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    bin_size = 150
    print('Creating binned index with bin size %d' % bin_size)
    bin_index_creation_time = create_binned_index(config, vcfc_filename, bin_size)
    print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))

    vcfc_binned_durations = []
    print('Running vcfc_binned_exhaustive')
    for pos in range(min_pos, max_pos+1, step):
        print('vcfc_binned_exhaustive: %d' % pos)
        durs = []
        for _ in range(test_runs):
            durs.append(run_vcfc_binned_index_query(config, vcfc_filename, '22', pos, pos))
        # vcfc_binned_durations.append((pos, pos, sum(durs)/len(durs)))
        vcfc_binned_durations.append({
            'start_position': pos,
            'end_position': pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })

    data['vcfc_binned_exhaustive_%d' % bin_size] = {
        'data': vcfc_binned_durations,
        'label': 'VCFC Binned Index (bin size=%d)' % bin_size
    }


    bgzip_durations = []
    print('Running bgzip_exhaustive')
    for pos in range(min_pos, max_pos+1, step):
        print('bgzip_exhaustive: %d' % pos)
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bgzip_filename, '22', pos, pos))
        # bgzip_durations.append((pos, pos, sum(durs)/len(durs)))
        bgzip_durations.append({
            'start_position': pos,
            'end_position': pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })

    data['bgzip_exhaustive'] = {
        'data': bgzip_durations,
        'label': 'BGZIP + Tabix'
    }

    bcf_durations = []
    print('Running bcf_exhaustive')
    for pos in range(min_pos, max_pos+1, step):
        print('bcf_exhaustive: %d' % pos)
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bcf_filename, '22', pos, pos))
        # bcf_durations.append((pos, pos, sum(durs)/len(durs)))
        bcf_durations.append({
            'start_position': pos,
            'end_position': pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })

    data['bcf_exhaustive'] = {
        'data': bcf_durations,
        'label': 'BCF + Tabix'
    }

    return {
        'data': data,
        'title': 'VCFC Binned Index, Single Variant Lookup',
        'name': 'binned-exhaustive-single',
        'xlabel': 'Variant Position',
        'ylabel': 'Time (seconds)'
    }

def graph_binned_index_exhaustive(measurements):
    graph_measurements(
        measurements,
        lambda val: val['start_position'],
        lambda val: val['time'])


def measure_binned_index_range_queries(query_range:int=10000, queries:int=100):
    data = {}
    assert queries > 0, 'queries > 0'
    # Set step to fit `queries` number of positions into the range(min_pos, max_pos)
    step = int((max_pos - query_range - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    bin_size = 150
    print('Creating binned index with bin size %d' % bin_size)
    bin_index_creation_time = create_binned_index(config, vcfc_filename, bin_size)
    print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))

    # VCFC measurements
    # Iterate through ranges starting at front of file
    vcfc_durations = []
    for pos in range(min_pos, (min_pos+step*queries)+1, step):
        print('binned_query: %d' % pos)
        end_pos = pos + query_range
        durs = []
        for _ in range(test_runs):
            durs.append(run_vcfc_binned_index_query(
                config, vcfc_filename, reference_name, pos, end_pos))

        vcfc_durations.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['vcfc_binned_exhaustive'] = {
        'data': vcfc_durations,
        'label': 'Sparse Offset-as-Index'
    }


    # Tabix+BGZIP measurements
    tabix_bgzip_durations = []
    for pos in range(min_pos, (min_pos+step*queries)+1, step):
        print('bgzf_query: %d' % pos)
        end_pos = pos + query_range
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bgzip_filename, reference_name, pos, end_pos))
        tabix_bgzip_durations.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['tabix_bgzip_exhaustive'] = {
        'data': tabix_bgzip_durations,
        'label': 'BGZF + Tabix'
    }


    # Tabix+BCF measurements
    tabix_bcf_durations = []
    for pos in range(min_pos, (min_pos+step*queries)+1, step):
        print('bcf_query: %d' % pos)
        end_pos = pos + query_range
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bcf_filename, reference_name, pos, end_pos))
        tabix_bcf_durations.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['tabix_bcf_exhaustive'] = {
        'data': tabix_bcf_durations,
        'label': 'BCF + Tabix'
    }


    return {
        'data': data,
        'title': 'VCFC Binned Index, Query Range %d' % query_range,
        'name': 'binned-exhaustive-range',
        'xlabel': 'Query Range Start Position',
        'ylabel': 'Time (seconds)'
    }


def graph_measurements(measurements, xfunc, yfunc, stddev_label=False):
    coloridx = 0
    plots = []
    legend_labels = []

    max_yval = 0

    for k in measurements['data']:
        measurement = measurements['data'][k]
        print('measurement:\n' + json.dumps(measurement))
        secondarycolor, primarycolor = colors[coloridx]
        # secondarycolor = colors[coloridx][0]
        coloridx += 1
        xvals = []
        yvals = []
        for val in measurement['data']:
            xvals.append(xfunc(val))
            yvals.append(yfunc(val))
        print('Plotting %s' % measurement['label'])

        # Update max_yval
        if max(yvals) > max_yval:
            max_yval = max(yvals)
        # print('xvals:\n' + str(xvals))
        # print('yvals:\n' + str(yvals))

        marker_size = 1
        p1 = plt.scatter(xvals, yvals, s=marker_size, label=measurement['label'], color=primarycolor)
        plots.append(p1)
        legend_labels.append(measurement['label'])

        z = np.polyfit(xvals, yvals, 1)
        p = np.poly1d(z)
        p2 = plt.plot(xvals, p(xvals), 'r--', linewidth=2, label=measurement['label'], color=secondarycolor)
        plots.append(p2)
        legend_labels.append(measurement['label'])

    plt.xlabel(measurements['xlabel'])
    plt.ylabel(measurements['ylabel'])
    # plt.xticks(rotation=90)
    plt.ylim(bottom=0, top=max_yval)
    plt.legend()
    plt.title(measurements['title'])

    # Tighten layout
    plt.margins(0.01)
    plt.tight_layout()

    if 'name' in measurements:
        filename = measurements['name'] + '.png'
    else:
        filename = measurements['title'].replace(' ', '_') + '.png'
    plt.savefig(filename, dpi=200)
    print('Saved figure: %s' % filename)


# TODO use new graph_measurements function
def measure_sparse_single_variant_exhaustive(queries:int=200):
    data = {}
    step = int((max_pos - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    sparse_durations = []
    print('Running vcfc_sparse_exhaustive')
    for pos in range(min_pos, max_pos+1, step):
        print('vcfc_sparse_exhaustive: %d' % pos)
        durs = []
        for _ in range(test_runs):
            durs.append(run_vcfc_sparse_query(config, sparse_filename, '22', pos, pos))
        sparse_durations.append({
            'start_position': pos,
            'end_position': pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['vcfc_sparse_exhaustive'] = {
        'data': sparse_durations,
        'label': 'VCFC Sparse Offset-as-Index'
    }

    bgzip_durations = []
    print('Running bgzip_exhaustive')
    for pos in range(min_pos, max_pos+1, step):
        print('bgzip_exhaustive: %d' % pos)
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bgzip_filename, '22', pos, pos))
        bgzip_durations.append({
            'start_position': pos,
            'end_position': pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['bgzip_exhaustive'] = {
        'data': bgzip_durations,
        'label': 'BGZF + Tabix'
    }

    bcf_durations = []
    print('Running bcf_exhaustive')
    for pos in range(min_pos, max_pos+1, step):
        print('bcf_exhaustive: %d' % pos)
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bcf_filename, '22', pos, pos))
        bcf_durations.append({
            'start_position': pos,
            'end_position': pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['bcf_exhaustive'] = {
        'data': bcf_durations,
        'label': 'BCF + Tabix'
    }

    return {
        'data': data,
        'title': 'VCFC Sparse Offset-as-Index, Single Variant Lookup',
        'name': 'sparse-exhaustive-single',
        'xlabel': 'Variant query location',
        'ylabel': 'Time (seconds)'
    }


def graph_sparse_single_variant_exhaustive(measurements):
    graph_measurements(measurements, lambda val:val['start_position'], lambda val:val['time'])


def measure_sparse_range_queries(query_range:int=10000, queries:int=100):
    data = {}
    assert queries > 0, 'queries > 0'
    # Set step to fit `queries` number of positions into the range(min_pos, max_pos)
    step = int((max_pos - query_range - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    # measure percentage of search space for each, not whole file
    # iter_max = int((max_pos - min_pos) / step * 0.05)

    # VCFC measurements
    # Iterate through ranges starting at front of file
    sparse_durations = []
    for pos in range(min_pos, (min_pos+step*queries)+1, step):
        print('sparse_query: %d' % pos)
        end_pos = pos + query_range
        durs = []
        for _ in range(test_runs):
            durs.append(run_vcfc_sparse_query(config, sparse_filename, reference_name, pos, end_pos))

        sparse_durations.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['vcfc_sparse_exhaustive'] = {
        'data': sparse_durations,
        'label': 'Sparse Offset-as-Index'
    }


    # Tabix+BGZIP measurements
    tabix_bgzip_durations = []
    for pos in range(min_pos, (min_pos+step*queries)+1, step):
        print('bgzf_query: %d' % pos)
        end_pos = pos + query_range
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bgzip_filename, reference_name, pos, end_pos))
        tabix_bgzip_durations.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['tabix_bgzip_exhaustive'] = {
        'data': tabix_bgzip_durations,
        'label': 'BGZF + Tabix'
    }


    # Tabix+BCF measurements
    tabix_bcf_durations = []
    for pos in range(min_pos, (min_pos+step*queries)+1, step):
        print('bcf_query: %d' % pos)
        end_pos = pos + query_range
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bcf_filename, reference_name, pos, end_pos))
        tabix_bcf_durations.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['tabix_bcf_exhaustive'] = {
        'data': tabix_bcf_durations,
        'label': 'BCF + Tabix'
    }


    return {
        'data': data,
        'title': 'VCFC Sparse Offset-as-Index, Query Range %d' % query_range,
        'name': 'sparse-exhaustive-range',
        'xlabel': 'Query Range Start Position',
        'ylabel': 'Time (seconds)'
    }


def measure_sparse_range_queries_old(queries:int=200):
    data = {}
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)
    # step = 100000
    step = int((max_pos - min_pos) / queries)

    # measure percentage of search space for each, not whole file
    iter_max = int((max_pos - min_pos) / step * 0.05)

    # VCFC measurements
    # Iterate through ranges starting at front of file
    sparse_values_from_front = []
    iteration = 0
    for end_pos in range(min_pos, max_pos+1, step):
        if iteration >= iter_max:
            break
        iteration += 1
        durs = []
        for _ in range(test_runs):
            durs.append(run_vcfc_sparse_query(config, sparse_filename, reference_name, min_pos, end_pos))

        sparse_values_from_front.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['sparse_values_from_front'] = {
        'data': sparse_values_from_front,
        'label': 'Sparse Offset-as-Index, front of file'
    }


    # Iterate through ranges starting at back of file
    sparse_values_from_back = []
    iteration = 0
    for start_pos in range(max_pos, min_pos, -step):
        if iteration >= iter_max:
            break
        iteration += 1
        durs = []
        for _ in range(test_runs):
            durs.append(run_vcfc_sparse_query(config, sparse_filename, reference_name, start_pos, max_pos))
        sparse_values_from_back.append({
            'start_position': start_pos,
            'end_position': max_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['sparse_values_from_back'] = {
        'data': sparse_values_from_back,
        'label': 'Sparse Offset-as-Index, back of file'
    }


    # Iterate through ranges starting at mid of file
    sparse_values_from_mid = []
    mid_pos = int((min_pos + max_pos + 1)/2)
    iteration = 0
    for start_pos in range(mid_pos, 0, -int(step/2)):
        if iteration >= iter_max:
            break
        iteration += 1
        end_pos = (mid_pos-start_pos) + mid_pos
        durs = []
        for _ in range(test_runs):
            durs.append(run_vcfc_sparse_query(config, sparse_filename, reference_name, start_pos, end_pos))
        sparse_values_from_mid.append({
            'start_position': start_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['sparse_values_from_mid'] = {
        'data': sparse_values_from_mid,
        'label': 'Sparse Offset-as-Index, center of file'
    }


    # Tabix+BGZIP measurements
    tabix_bgzip_values_from_front = []
    iteration = 0
    for end_pos in range(min_pos, max_pos, step):
        if iteration >= iter_max:
            break
        iteration += 1
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bgzip_filename, reference_name, min_pos, end_pos))
        tabix_bgzip_values_from_front.append({
            'start_position': min_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['tabix_bgzip_values_from_front'] = {
        'data': tabix_bgzip_values_from_front,
        'label': 'BGZF + Tabix, back of file'
    }


    tabix_bgzip_values_from_back = []
    iteration = 0
    for start_pos in range(max_pos, min_pos, -step):
        if iteration >= iter_max:
            break
        iteration += 1
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bgzip_filename, reference_name, start_pos, max_pos))
        tabix_bgzip_values_from_back.append({
            'start_position': start_pos,
            'end_position': max_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['tabix_bgzip_values_from_back'] = {
        'data': tabix_bgzip_values_from_back,
        'label': 'BGZF + Tabix, back of file'
    }


    # Iterate through ranges starting at mid of file
    tabix_bgzip_values_from_mid = []
    mid_pos = int((min_pos + max_pos + 1)/2)
    iteration = 0
    for start_pos in range(mid_pos, 0, -int(step/2)):
        if iteration >= iter_max:
            break
        iteration += 1
        end_pos = (mid_pos-start_pos) + mid_pos
        durs = []
        for _ in range(test_runs):
            durs.append(run_tabix(config, bgzip_filename, reference_name, start_pos, end_pos))
        tabix_bgzip_values_from_mid.append({
            'start_position': start_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    data['tabix_bgzip_values_from_mid'] = {
        'data': tabix_bgzip_values_from_mid,
        'label': 'BGZF + Tabix, center of file'
    }


    return {
        'data': data,
        'title': 'VCFC Binned Index, Range Lookup',
        'name': 'sparse-exhaustive-range',
        'xlabel': 'Query Range',
        'ylabel': 'Time (seconds)'
    }


def graph_sparse_range_queries(measurements):
    graph_measurements(measurements, lambda val: val['end_position']-val['start_position'], lambda val: val['time'])


def measure_binned_index_time_profile(queries:int=100):
    data = {}
    # step = 50000
    step = int((max_pos - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    bin_sizes = [
        *range(5, 50, 5),
        *range(50, 200, 10),
        *range(200, 2000+1, 25)
    ]

    # for bin_size in range(min_bin_size, max_bin_size+bin_step+1, bin_step):
    for bin_size in bin_sizes:
        bin_profile = {}
        print('Creating binned index with bin size %d' % bin_size)
        bin_index_creation_time = create_binned_index(config, vcfc_filename, bin_size)
        print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))

        # Run regular queries, aggregate profile after each

        # vcfc_binned_durations = []
        print('Running exhaustive binned queries')
        test_count = 0
        for pos in range(min_pos, max_pos+1, step):
            test_count += 1
            print('vcfc_binned_exhaustive: %d' % pos)

            profiles = []
            for _ in range(test_runs):
                profiles.append(run_vcfc_binned_index_timing_profile(config, vcfc_filename, reference_name, pos, pos))

            # Merge timing profiles
            for k in profiles[0]:
                if k not in bin_profile:
                    bin_profile[k] = 0
                for p in profiles:
                    bin_profile[k] += p[k] / test_runs

        for k in bin_profile:
            bin_profile[k] /= test_count

        data['vcfc_binned_index_%d' % bin_size] = {
            'data': bin_profile,
            'label': 'Bin Size %d' % bin_size
        }


    return {
        'data': data,
        'title': 'VCFC Binned Index Query Phase Time Profile, Single Variant Lookup',
        'name': 'binned-timing-profile',
        'xlabel': 'Bin Size (not linear scale)',
        'ylabel': 'Time (seconds)'
    }

def measure_binned_index_time_profile_range(query_range:int=10000, queries:int=100):
    data = {}
    assert queries > 0, 'queries > 0'
    # Set step to fit `queries` number of positions into the range(min_pos, max_pos)
    step = int((max_pos - query_range - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    bin_sizes = [
        *range(5, 50, 5),
        *range(50, 200, 10),
        *range(200, 2000+1, 25)
    ]

    # for bin_size in range(min_bin_size, max_bin_size+bin_step+1, bin_step):
    for bin_size in bin_sizes:
        bin_profile = {}
        print('Creating binned index with bin size %d' % bin_size)
        bin_index_creation_time = create_binned_index(config, vcfc_filename, bin_size)
        print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))

        # Run regular queries, aggregate profile after each
        print('Running %d exhaustive binned queries of size %d')
        # get exactly `queries` loops, account for rounding
        test_count = 0
        for pos in range(min_pos, (min_pos+step*queries)+1, step):
            test_count += 1
            print('vcfc_binned_exhaustive: %d' % pos)
            # profile = run_vcfc_binned_index_timing_profile(config, vcfc_filename, reference_name, pos, pos+step)
            profiles = []
            for test_run in range(test_runs):
                print('%d ' % test_run, end='')
                profiles.append(run_vcfc_binned_index_timing_profile(config, vcfc_filename, reference_name, pos, pos+query_range))

            # Merge timing profiles
            for k in profiles[0]:
                if k not in bin_profile:
                    bin_profile[k] = 0
                for p in profiles:
                    bin_profile[k] += p[k] / test_runs

        for k in bin_profile:
            bin_profile[k] /= test_count

        data['vcfc_binned_index_%d' % bin_size] = {
            'data': bin_profile,
            'label': 'Bin Size %d' % bin_size
        }

    return {
        'data': data,
        'title': 'VCFC Binned Index Interval, Query Range %d' % query_range,
        'name': 'binned-timing-profile-range',
        'xlabel': 'Bin Size (not linear scale)',
        'ylabel': 'Time (seconds)'
    }

def graph_binned_index_time_profile(measurements):
    '''
    See measure_binned_index_time_profile and binned-timing-profile.json for measurements schema.
    '''

    data = measurements['data']

    time_categories = [
        {
            'flags': ['decompress2_metadata_headers', 'decompress2_metadata_headers_fd'],
            'label': 'Reading Headers'
        },
        {
            'flags': ['decompress_iteration'],
            'label': 'Line Decompression'
        },
        {
            'flags': ['index_search'],
            'label': 'Index Search'
        },
        # {
        #     'flags': ['get_or_set_decompressed_cache'],
        #     'label': 'Run Decompression or Cache Retrieval'
        # },
        {
            'flags': ['decompress_seeking'],
            'label': 'Seeking From Start of Bin'
        },
    ]

    xvals = []
    yvals = [] # objects matching above

    for entry_key in data:
        entry = data[entry_key]
        bin_size = re.compile(r'.*_(\d+)$').match(entry_key).group(1)

        xvals.append(bin_size)
        yval = {}
        for time_category in time_categories:
            category_time = 0
            for flag in time_category['flags']:
                category_time += entry['data'].get(flag, 0)
            yval[time_category['label']] = category_time
        yvals.append(yval)

    print(yvals)

    combined_yvals = []
    yval_labels = []

    for time_category in time_categories:
        category_yvals = [v[time_category['label']] for v in yvals]
        combined_yvals.append(category_yvals)
        yval_labels.append(time_category['label'])

    # yval_labels = list(reversed(yval_labels))

    plots = []
    cumulative_yvals = None
    max_cumulative_yval = 0
    color_idx = 0
    for i in range(0, len(combined_yvals)):
        if i == 0:
            p = plt.bar(xvals, combined_yvals[i], color=colors[color_idx][1])
            cumulative_yvals = combined_yvals[i]
        else:
            # Stack on previous
            p = plt.bar(xvals, combined_yvals[i], bottom=copy.copy(cumulative_yvals), color=colors[color_idx][1])
            # Add to cumulative
            cumulative_yvals = [cumulative_yvals[j] + combined_yvals[i][j] for j in range(len(cumulative_yvals))]

        color_idx += 1

        if max(cumulative_yvals) > max_cumulative_yval:
            max_cumulative_yval = max(cumulative_yvals)

        plots.append(p)

    xticks = [xvals[i] for i in range(0, len(xvals), 2)]

    plt.xlabel(measurements['xlabel'])
    plt.ylabel(measurements['ylabel'])
    plt.xticks(xticks, rotation=270)
    plt.ylim(top=max_cumulative_yval*1.1)
    # plt.legend()
    plt.legend(reversed(plots), reversed(yval_labels))
    plt.title(measurements['title'])
    if 'name' in measurements:
        filename = measurements['name'] + '.png'
    else:
        filename = measurements['title'].replace(' ', '_') + '.png'

    # Tighten layout
    plt.margins(0.01)
    plt.tight_layout()

    plt.savefig(filename, dpi=200)
    print('Saved figure: %s' % filename)
    # plt.show()


def main(args):
    ops = [
        ['sparse-exhaustive-range', measure_sparse_range_queries, graph_sparse_range_queries],
        ['sparse-exhaustive-single', measure_sparse_single_variant_exhaustive, graph_sparse_single_variant_exhaustive],

        ['binned-exhaustive-single', measure_binned_index_single_exhaustive, graph_binned_index_exhaustive],
        ['binned-exhaustive-range', measure_binned_index_range_queries, graph_binned_index_exhaustive],

        ['binned-timing-profile', measure_binned_index_time_profile, graph_binned_index_time_profile],
        ['binned-timing-profile-range', measure_binned_index_time_profile_range, graph_binned_index_time_profile],

        ['binned-index-creation-time', measure_binned_index_creation_time, graph_binned_index_creation_time]
    ]

    def print_usage():
        print('Usage: range_query.py opname [measure|graph]')
        print('opnames:\n' + '\n'.join(o[0] for o in ops))

    if len(args) < 2:
        print_usage()
        return -1

    if args[0] not in [o[0] for o in ops]:
        print_usage()
        raise RuntimeError('Unrecognized subcommand')

    for op in ops:
        opname = op[0]
        measure_func = op[1]
        graph_func = op[2]
        if args[0] == opname:
            filename = opname + '.json'
            if args[1] == 'measure':
                data = measure_func()
                with open(filename, 'w') as f:
                    json.dump(data, f)
                print('Wrote file: ' + filename)
            elif args[1] == 'graph':
                with open(filename) as f:
                    data = json.load(f)
                graph_func(data)


if __name__ == '__main__':
    exit(main(sys.argv[1:]))
