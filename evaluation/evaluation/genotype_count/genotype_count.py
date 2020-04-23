import sys, os, re
import json
import gzip
import matplotlib.pyplot as plt

plt.rcParams.update({
    'figure.dpi': 200,
    'savefig.dpi': 200,
    'figure.figsize': (10.0, 5.0),
    'axes.labelsize': 'small',
    'font.size': 9
})

def main(args):
    if len(args) < 2:
        raise RuntimeError('Need command and filename')
    command = args[0]
    filename = args[1]

    if not os.path.exists(filename) or not os.path.isfile(filename):
        raise RuntimeError('Path is not a file: ' + str(filename))

    if command == 'graph-runs-by-line':
        graph_genotype_runs(filename)
    elif command == 'counts-by-line':
        # graph_genotype_counts_by_line(filename)
        pass
    elif command == 'graph_genotype_counts_aggregate':
        graph_genotype_counts_aggregate(filename)
    elif command == 'collect_genotype_counts_old':
        collect_genotype_counts_old(filename)
    elif command == 'graph_genotype_counts_old':
        graph_genotype_counts_old()
    else:
        raise RuntimeError('Unknown command: ' + command)



def collect_genotype_counts_old(file_name):
    pattern = re.compile(r'.*.vcf$')

    counts = {}
    in_data = False

    #for file_name in os.listdir(directory_name):
    if pattern.match(file_name) is None:
        return

    print('Reading file: ' + file_name)
    if file_name.startswith('/'):
        path = file_name
    else:
        path = os.path.join(os.getcwd(), file_name)
    # path = '%s/%s' % (directory_name, file_name)
    # if os.path.isfile(path):
    with open(path) as f:
        for line in f:
            if in_data or not line.startswith('#'):
                in_data = True
            elif in_data and line.startswith('#'):
                raise RuntimeError('Got header line after data line')
            elif not in_data and line.startswith('#'):
                continue

            if in_data:
                terms = line.strip().split('\t')
                if len(terms) <= 9:
                    raise RuntimeError('line didn\'t have enough terms: ' + line)
                if terms[8] != 'GT':
                    print('FORMAT was not GT: ' + terms[8])
                    # raise RuntimeError('FORMAT was not GT: ' + terms[8])
                for sample_i in range(9, len(terms)):
                    sample_gt = terms[sample_i]
                    if sample_gt not in counts:
                        counts[sample_gt] = 0
                    counts[sample_gt] += 1

    with open('genotype_counts.json', 'w') as f_out:
        json.dump(counts, f_out)



def graph_genotype_counts_old():
    with open('genotype_counts.json') as f:
        counts = json.load(f)

    items = list(counts.items())
    items = sorted(items, key=lambda it: it[0]) # sort based on key, genotype value

    xvals = [it[0] for it in items]
    yvals = [it[1] for it in items]

    def gt_numeric(gt):
        terms = gt.split('|')
        return int(terms[0]) + int(terms[1])

    filtered = [(x,y) for x,y in zip(xvals, yvals) if y > 3]
    filtered = sorted(filtered, key=lambda v:gt_numeric(v[0]))
    xvals = [f[0] for f in filtered]
    yvals = [f[1] for f in filtered]

    plt.bar(xvals, yvals)
    plt.yscale('log')

    plt.xlabel('Sample Variant Genotype')
    plt.ylabel('Number of occurrences in Chromosome 1 VCF')
    plt.title('Occurrences of Sample Genotypes Within the 1000 Genomes VCF dataset')

    # Tighten layout
    plt.margins(0.01)
    plt.tight_layout()

    plt.xticks(rotation=270)
    plt.grid(axis='y')

    image_file_name = 'genotype_counts.png'
    plt.savefig(image_file_name)
    print('Saved file: ' + image_file_name)


def graph_genotype_counts_aggregate(filename):
    # sample gt -> count
    gt_counts = {}
    sample_pattern = re.compile(r'\d+\|\d+')

    with open(filename) as f:
        for line in f:
            terms = line.split(' ')
            assert len(terms) == 2, 'line did not have 2 terms %s' % line
            sample = terms[0]
            count = int(terms[1])
            if sample_pattern.match(terms[0]):
                if sample not in gt_counts:
                    gt_counts[sample] = 0
                gt_counts[sample] += count

    # xvals = [v[0] for v in vals]
    # yvals = [v[1] for v in vals]
    xvals = gt_counts.keys()
    yvals = [gt_counts[x] for x in xvals]

    plt.bar(xvals, yvals)
    image_file_name = 'genotype-counts-aggregate.png'
    plt.savefig(image_file_name)
    print('Saved file: ' + image_file_name)
    os.system('xdg-open ' + image_file_name)


def graph_genotype_runs(filename):
    # sample gt -> run lengths -> run length count
    sample_run_lengths = {}
    sample_pattern = re.compile(r'\d+\|\d+')

    line_counter = 0
    line_counter_print_interval = 100000

    with open(filename) as f:
        for line in f:
            if line_counter % line_counter_print_interval == 0:
                print('line_counter = %d' % line_counter)
            line_counter += 1

            terms = line.split(' ')
            if len(terms) != 2:
                raise RuntimeError('line did not have 2 terms %s' % line)
            sample = terms[0]
            run_length = int(terms[1])
            if sample_pattern.match(sample):
                # Initialize sample and or run length counter
                if sample not in sample_run_lengths:
                    sample_run_lengths[sample] = {}
                if run_length not in sample_run_lengths[sample]:
                    sample_run_lengths[sample][run_length] = 0

                # Increment run length for this sample gt
                sample_run_lengths[sample][run_length] += 1


    xvals = []
    yvals = []
    marker_sizes = []
    max_marker_size = 0

    yval_bin_size = 25

    for sample in sample_run_lengths:
        print(sample)

        run_lengths = sample_run_lengths[sample].keys() # run lengths

        # Purge any with runs of only 1 length (not really a "run")
        min_run_length_threshold = 2
        run_lengths = [rl for rl in run_lengths if rl >= min_run_length_threshold]
        if len(run_lengths) == 0:
            continue

        run_marker_sizes = [sample_run_lengths[sample][rl] for rl in run_lengths]

        merged = merge_key_bins(run_lengths, run_marker_sizes, yval_bin_size)

        run_lengths = [m[0] for m in merged]
        run_marker_sizes = [m[1] for m in merged]

        xvals += [sample] * len(run_lengths)
        yvals += run_lengths
        marker_sizes += run_marker_sizes

    max_marker_size = max(marker_sizes)

    print(list(zip(xvals, yvals, marker_sizes)))

    # scale to N as max marker size
    N = 100 # uses base 10
    minyval = min(yvals)
    maxyval = max(yvals)
    # import numpy as np
    # logvals = np.logspace(np.log10(1), np.log10(N), max_marker_size+1)
    # # logvals = np.geomspace(1, N, max_marker_size+1)
    # max_logval = max(logvals)
    # logvals = [(max_logval - i) for i in logvals]
    # # print('logvals:\n' + str(logvals))
    # marker_sizes = [logvals[y] for y in marker_sizes]
    # print(marker_sizes)

    # alphas = [(1-(y/maxyval)) for y in yvals]
    # print(alphas)

    # limit_marker_to = 2000
    import numpy as np
    markersize_min = 5
    markersize_max = 1000
    markersize_map = np.linspace(markersize_min, markersize_max, max(marker_sizes)+1)
    # markersize_map = np.geomspace(markersize_min, markersize_max, max(marker_sizes)+1)

    # marker_sizes = [limit_marker_to * (m/max_marker_size) for m in marker_sizes]
    adjusted_marker_sizes = [markersize_map[m] for m in marker_sizes]

    # set true zeros to zeros
    for i in range(len(yvals)):
        if marker_sizes[i] == 0:
            adjusted_marker_sizes[i] = 0

    print(list(zip(xvals, yvals, marker_sizes)))

    plt.scatter(xvals, yvals, s=adjusted_marker_sizes, alpha=0.8, marker='o')
    ytick_step = int((maxyval-minyval)/25)
    yticks = [y for y in range(minyval, maxyval+2*ytick_step+1, ytick_step)]
    # print(yticks)
    plt.yticks(yticks)

    plt.xlabel('Sample Variant Genotype')
    plt.ylabel('Run Length, Binned by %d-length intervals (size is frequency of run length)' % yval_bin_size)
    plt.title('Run Lengths by Sample Genotype Within the 1000 Genomes VCF dataset')

    # Tighten layout
    plt.margins(0.05)
    plt.tight_layout()

    image_file_name = 'genotype-run-lengths.png'
    plt.savefig(image_file_name, dpi=300)
    print('Saved file: ' + image_file_name)
    # os.system('xdg-open ' + image_file_name)

    # plt.show()


def graph_genotype_runs_linegraph(filename):
    # sample gt -> run lengths -> run length count
    sample_run_lengths = {}
    sample_pattern = re.compile(r'\d+\|\d+')

    line_counter = 0
    line_counter_print_interval = 100000

    with open(filename) as f:
        for line in f:
            if line_counter % line_counter_print_interval == 0:
                print('line_counter = %d' % line_counter)
            line_counter += 1

            terms = line.split(' ')
            assert len(terms) == 2, 'line did not have 2 terms %s' % line
            sample = terms[0]
            run_length = int(terms[1])
            if sample_pattern.match(sample):
                # Initialize sample and or run length counter
                if sample not in sample_run_lengths:
                    sample_run_lengths[sample] = {}
                if run_length not in sample_run_lengths[sample]:
                    sample_run_lengths[sample][run_length] = 0

                # Increment run length for this sample gt
                sample_run_lengths[sample][run_length] += 1


    marker_sizes = []
    max_marker_size = 0
    for sample in sample_run_lengths:
        xvals = []
        yvals = []

        for run_length in sample_run_lengths[sample]:
            run_length_count = sample_run_lengths[sample][run_length]
            xvals.append(run_length)
            yvals.append(run_length_count)

        plt.scatter(xvals, yvals, label=sample, s=1)



    # plt.scatter(xvals, yvals, s=marker_sizes, alpha=0.8, marker='_')
    # ytick_step = int((maxyval-minyval)/25)
    # yticks = [y for y in range(minyval, maxyval+2*ytick_step+1, ytick_step)]
    # print(yticks)
    # plt.yticks(yticks)

    # Tighten layout
    plt.margins(0.01)
    plt.tight_layout()


    plt.legend()
    plt.yscale('log')

    image_file_name = 'genotype-run-lengths.png'
    plt.savefig(image_file_name, dpi=200)
    print('Saved file: ' + image_file_name)
    # os.system('xdg-open ' + image_file_name)

    plt.show()


def merge_key_bins(keys, values, bin_size:int=10):
    items = sorted(zip(keys, values), key=lambda v:v[0])
    output_list = []
    items_idx = 0
    for bin_start in range(0, max(keys)+1, bin_size):
        # print(bin_start)
        bin_total_value = 0
        while items_idx < len(items):
            if items[items_idx][0] < bin_start + bin_size:
                bin_total_value += items[items_idx][1]
                # print(items[items_idx][1])
                items_idx += 1
            else:
                break
        output_list.append((bin_start, bin_total_value))
    return output_list


def chunk_list(input_list, chunk_size):
    '''
    This function takes the input list and returns a list with adjacent groups of size
    chunk_list combined into tuples. The last tuple may be smaller than chunk_size.

    Example: input_list = [a, b, c, d, e, f, g], chunk_size = 3
        returns = [(a, b, c), (d, e, f), (g,)]
    '''
    '''
    for i in range(0, int(len(bin_id_list)/chunk_size + 2)):
        row_list = []
        for ci in range(0, chunk_size):
            idx = i*chunk_size + ci
            if idx >= len(bin_id_list):
                break
            row_list.append(bin_id_list[i*chunk_size + ci])
        chunked_bin_ids.append(tuple(row_list))
    '''
    output = []
    for i in range(0, int(len(input_list)/chunk_size + 2)):
        row_list = []
        for ci in range(0, chunk_size):
            idx = i*chunk_size + ci
            if idx >= len(input_list):
                break
            row_list.append(input_list[idx])
        if len(row_list) > 0:
            output.append(tuple(row_list))
    return output


if __name__ == '__main__':
    exit(main(sys.argv[1:]))
