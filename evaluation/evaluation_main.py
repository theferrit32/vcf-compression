#!/bin/env python3
import os, time, json, sys
import re
import copy
import subprocess
import numpy as np
import humanize

from evaluation.plt_common import plt
from evaluation.plt_common import colors

from evaluation.config import Config
from evaluation.command import (
    run_tabix,
    create_tabix_index,

    run_vcfc_sparse_query,

    create_vcfc_sparse_external_index,
    run_vcfc_sparse_external_index_query,

    create_binned_index,
    run_vcfc_binned_index_query,
    run_vcfc_binned_index_timing_profile,
)

vcfc_dir = '/home/krferrit/vcf-compression/'
tabix_cmd = '/home/krferrit/vcf-compression/tabix'
bgzip_cmd = '/home/krferrit/vcf-compression/bgzip'

default_binsize = 150
test_runs = 10

# Chromosome 22
vcfc_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.vcfc'
sparse_filename = vcfc_filename + '.sparse'
bgzip_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
bcf_filename = '/mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf'
reference_name = '22'
min_pos = 16050075
max_pos = 51244237

# Chromosome 1
# file_directory = '/mnt/ext4'
# vcfc_filename = file_directory + '/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.vcfc'
# sparse_filename = vcfc_filename + '.sparse'
# bgzip_filename = file_directory + '/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
# bcf_filename = file_directory + '/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.bcf'
# reference_name = '1'
# min_pos = 10177
# max_pos = 249240543



def measure_indexing_times():
    data = {}
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    bin_size = 150
    vcfc_binned_index_times = []
    for _ in range(test_runs):
        vcfc_binned_index_times.append(
            create_binned_index(config, vcfc_filename, bin_size))
    data['vcfc_binned_index'] = {
        'times': vcfc_binned_index_times,
        'mean': np.average(vcfc_binned_index_times),
        'stddev': np.std(vcfc_binned_index_times)
    }

    vcfc_sparse_external_index_times = []
    for _ in range(test_runs):
        vcfc_sparse_external_index_times.append(
            create_vcfc_sparse_external_index(config, vcfc_filename))
    data['vcfc_sparse_external_index'] = {
        'times': vcfc_sparse_external_index_times,
        'mean': np.average(vcfc_sparse_external_index_times),
        'stddev': np.std(vcfc_sparse_external_index_times)
    }


    bgzip_index_times = []
    for _ in range(test_runs):
        bgzip_index_times.append(
            create_tabix_index(config, bgzip_filename))
    data['bgzip_index'] = {
        'times': bgzip_index_times,
        'mean': np.average(bgzip_index_times),
        'stddev': np.std(bgzip_index_times)
    }


    bcf_index_times = []
    for _ in range(test_runs):
        bcf_index_times.append(
            create_tabix_index(config, bcf_filename))
    data['bcf_index'] = {
        'times': bcf_index_times,
        'mean': np.average(bcf_index_times),
        'stddev': np.std(bcf_index_times)
    }
    return data

def graph_indexing_times(data):
    fmt = 'o'
    plt.errorbar(['VCFC Binned'],
            [data['vcfc_binned_index']['mean']],
            # vcfc_binned_index_results['mean'],
            yerr=data['vcfc_binned_index']['stddev'],
            label='VCFC Binned External Index (bin_size=150)',
            fmt=fmt)
    plt.errorbar(['VCFC Sparse External'],
            [data['vcfc_sparse_external_index']['mean']],
            yerr=data['vcfc_sparse_external_index']['stddev'],
            label='VCFC Sparse External Index',
            fmt=fmt)
    plt.errorbar(['BGZIP + Tabix'],
            [data['bgzip_index']['mean']],
            yerr=data['bgzip_index']['stddev'],
            label='BGZIP Compression with Tabix Index',
            fmt=fmt)
    plt.errorbar(['BCF + Tabix'],
            [data['bcf_index']['mean']],
            yerr=data['bcf_index']['stddev'],
            label='BCF Compression with Tabix Index',
            fmt=fmt)

    plt.ylabel('Time (seconds)')
    plt.xlabel('Compression/Indexing Strategy')
    plt.title('Index Creation Time')
    plt.legend()

    plt.xticks(rotation=30)
    plt.grid(axis='y')

    # Tighten layout
    plt.margins(0.01)
    plt.tight_layout()

    filename = 'index-creation-times.png'
    plt.savefig(filename)
    print('Saved figure: %s' % filename)



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
    measurement = measurements['data']['vcfc_binned_index_creation_time']['data']
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


def graph_measurements_single(measurements):
    return graph_measurements(
        measurements,
        lambda val:val['start_position'],
        lambda val:val['time'])

def graph_measurements_range(measurements):
    return graph_measurements(
        measurements,
        lambda val: val['end_position']-val['start_position'],
        lambda val: val['time'])

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


"""
Function to perform a single run against a qurey function with a file and lists of start/end positions
"""
def _run_variant_query_exhaustive(
        start_positions:list,
        end_positions:list,
        test_name:str,
        config:Config,
        function,
        filename:str):
    assert len(start_positions) == len(end_positions)

    total_durations = []
    print('Running %s' % test_name)
    for (start_pos, end_pos) in zip(start_positions, end_positions):
        print('%s: %d-%d' % (test_name, start_pos, end_pos))
        durs = []
        for _ in range(test_runs):
            durs.append(function(config, filename, '22', start_pos, end_pos))
        total_durations.append({
            'start_position': start_pos,
            'end_position': end_pos,
            'time': sum(durs)/len(durs),
            'stddev': np.std(durs)
        })
    return total_durations

def measure_all_range_variant(queries:int=200, query_range:int=10000):
    data = {}
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)
    # step = int((max_pos - min_pos) / queries)

    # Set step to fit `queries` number of positions into the range(min_pos, max_pos)
    step = int((max_pos - query_range - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    print('Creating binned index for binned index queries')
    bin_index_creation_time = create_binned_index(config, vcfc_filename, default_binsize)
    print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))

    # for pos in range(min_pos, (min_pos+step*queries)+1, step):
        # end_pos = pos + query_range

    start_positions = list(range(min_pos, (min_pos+step*queries)+1, step))
    end_positions = [s + query_range for s in start_positions]

    # VCFC Sparse External Index
    sparse_external_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'vcfc-sparse-external-exhaustive',
        config,
        run_vcfc_sparse_external_index_query,
        vcfc_filename
    )
    data['vcfc_sparse_external_exhaustive'] = {
        'data': sparse_external_durations,
        'label': 'VCFC Sparse External Index'
    }

    # VCFC Sparse Offset-as-Index
    sparse_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'vcfc-sparse-exhaustive',
        config,
        run_vcfc_sparse_query,
        sparse_filename
    )
    data['vcfc_sparse_exhaustive'] = {
        'data': sparse_durations,
        'label': 'VCFC Sparse Offset-as-Index'
    }

    # VCFC Binned External Index
    create_binned_index(config, vcfc_filename, default_binsize)
    binned_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'vcfc-binned-external-exhaustive',
        config,
        run_vcfc_binned_index_query,
        vcfc_filename
    )
    data['vcfc_sparse_external_exhaustive'] = {
        'data': binned_durations,
        'label': 'VCFC Binned External Index (Bin Size %d)' % default_binsize
    }

    # BGZIP with Tabix Index
    bgzip_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'bgzip-tabix-external-exhaustive',
        config,
        run_tabix,
        bgzip_filename
    )
    data['bgzip_tabix_exhaustive'] = {
        'data': bgzip_durations,
        'label': 'BGZIP + Tabix Index'
    }

    # BCF with Tabix Index
    bcf_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'bcf-tabix-external-exhaustive',
        config,
        run_tabix,
        bcf_filename
    )
    data['bcf_tabix_exhaustive'] = {
        'data': bcf_durations,
        'label': 'BCF + Tabix Index'
    }

    return {
        'data': data,
        'title': 'Variant Range Query Time (Query Range %d)' % query_range,
        'name': 'all-exhaustive-range',
        'xlabel': 'Range Start Position',
        'ylabel': 'Time (seconds)'
    }


def measure_all_single_variant(queries:int=200):
    data = {}
    step = int((max_pos - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    print('Creating binned index for binned index queries')
    bin_index_creation_time = create_binned_index(config, vcfc_filename, default_binsize)
    print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))

    start_positions = list(range(min_pos, max_pos+1, step))
    end_positions = start_positions


    # VCFC Sparse External Index
    sparse_external_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'vcfc-sparse-external-exhaustive',
        config,
        run_vcfc_sparse_external_index_query,
        vcfc_filename
    )
    data['vcfc_sparse_external_exhaustive'] = {
        'data': sparse_external_durations,
        'label': 'VCFC Sparse External Index'
    }

    # VCFC Sparse Offset-as-Index
    sparse_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'vcfc-sparse-exhaustive',
        config,
        run_vcfc_sparse_query,
        sparse_filename
    )
    data['vcfc_sparse_exhaustive'] = {
        'data': sparse_durations,
        'label': 'VCFC Sparse Offset-as-Index'
    }

    # VCFC Binned External Index
    create_binned_index(config, vcfc_filename, default_binsize)
    binned_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'vcfc-binned-external-exhaustive',
        config,
        run_vcfc_binned_index_query,
        vcfc_filename
    )
    data['vcfc_sparse_external_exhaustive'] = {
        'data': binned_durations,
        'label': 'VCFC Binned External Index (Bin Size %d)' % default_binsize
    }

    # BGZIP with Tabix Index
    bgzip_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'bgzip-tabix-external-exhaustive',
        config,
        run_tabix,
        bgzip_filename
    )
    data['bgzip_tabix_exhaustive'] = {
        'data': bgzip_durations,
        'label': 'BGZIP + Tabix Index'
    }

    # BCF with Tabix Index
    bcf_durations = _run_variant_query_exhaustive(
        start_positions, end_positions,
        'bcf-tabix-external-exhaustive',
        config,
        run_tabix,
        bcf_filename
    )
    data['bcf_tabix_exhaustive'] = {
        'data': bcf_durations,
        'label': 'BCF + Tabix Index'
    }

    return {
        'data': data,
        'title': 'Single Variant Query Time',
        'name': 'all-exhaustive-single',
        'xlabel': 'Variant Position',
        'ylabel': 'Time (seconds)'
    }


def measure_binned_index_time_profile(queries:int=100):
    data = {}
    step = int((max_pos - min_pos) / queries)
    config = Config(tabix_cmd, vcfc_dir, bgzip_cmd)

    bin_sizes = [
        *range(5, 50, 5),
        *range(50, 200, 10),
        *range(200, 2000+1, 50)
    ]

    for bin_size in bin_sizes:
        bin_profile = {}
        print('Creating binned index with bin size %d' % bin_size)
        bin_index_creation_time = create_binned_index(config, vcfc_filename, bin_size)
        print('Finished creating binned index, took %f seconds' % (bin_index_creation_time))

        # Run regular queries, aggregate profile after each
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
                    bin_profile[k] += p[k] / test_count

        for k in bin_profile:
            bin_profile[k] /= test_count

        data['vcfc_binned_index_%d' % bin_size] = {
            'data': bin_profile,
            'label': 'Bin Size %d' % bin_size
        }


    return {
        'data': data,
        'title': 'VCFC Binned Index Query Phase Time Profile, Single Variant Lookup',
        'name': 'binned-timing-profile-single',
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
        *range(200, 2000+1, 50)
    ]

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
            endpos = pos + query_range
            profiles = []
            for test_run in range(test_runs):
                print('%d ' % test_run, end='')
                profiles.append(
                    run_vcfc_binned_index_timing_profile(config, vcfc_filename, reference_name, pos, endpos))

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

    plt.savefig(filename)
    print('Saved figure: %s' % filename)


def main(args):
    ops = [
        ['all-exhaustive-single', measure_all_single_variant, graph_measurements_single],
        ['all-exhaustive-range', measure_all_range_variant, graph_measurements_range],

        ['binned-timing-profile-single', measure_binned_index_time_profile, graph_binned_index_time_profile],
        ['binned-timing-profile-range', measure_binned_index_time_profile_range, graph_binned_index_time_profile],

        ['binned-index-creation-time', measure_binned_index_creation_time, graph_binned_index_creation_time],

        ['all-indexing-times', measure_indexing_times, graph_indexing_times]
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
