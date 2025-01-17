import os, re, json, time, shutil
import subprocess
import requests
import sqlite3
import logging
import matplotlib.pyplot as plt
#cmd = tabix ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 16:60000-61000
# vcf.gz 1-22, X, Y = 17007956 KiB, 16.22 GiB
logging.basicConfig(level=logging.DEBUG)

bgzip_dir = '/home/me/dev/stanford/1000genomes/'
bcf_dir = '/home/me/dev/stanford/1000genomes/bcf/'
tabix_cmd = '/home/krferrit/vcf-compression/tabix'
vcfc_cmd = '/home/krferrit/vcf-compression/main_release'
sparse_vcfc_ext4_dir = '/mnt/ext4/'
sparse_vcfc_xfs_dir = '/mnt/xfs/'

ext4_workdir = '/mnt/ext4/'

def get_bgzip_file_path(reference_name):
    if reference_name == 'Y':
        return os.path.join(bgzip_dir, 'ALL.chrY.phase3_integrated_v2a.20130502.genotypes.vcf.gz')
    elif reference_name == 'X':
        return os.path.join(bgzip_dir, 'ALL.chrX.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz')
    else:
        basename = 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
        return os.path.join(bgzip_dir, basename % reference_name)

def get_bcf_file_path(reference_name):
    if reference_name in ('Y', 'X'):
        basename = 'ALL.chr%s.phase3_shapeit2_mvncall_integrated.20130502.genotypes.bcf'
        return os.path.join(bcf_dir, basename % reference_name)
    else:
        basename = 'ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.bcf'
        return os.path.join(bcf_dir, basename % reference_name)

def basename_in_dir(file_path, directory_path):
    if not file_path.startswith('/'):
        raise RuntimeError('file_path must be absolute: ' + str(file_path))
    if not directory_path.startswith('/'):
        raise RuntimeError('directory_path must be absolute: ' + str(directory_path))
    if file_path.startswith(directory_path):
        basename = os.path.basename(file_path)
        if os.path.join(directory_path, basename) == file_path:
            return True
    return False

def tabix(path, chrom, start, end):
    start = int(start)
    end = int(end)
    cmd = [
        tabix_cmd,
        path,
        '%s:%d-%d' % (chrom, start, end)
    ]
    s_time = time.time()
    proc = subprocess.Popen(cmd, stdin=None, stdout=subprocess.DEVNULL)
    proc.wait()
    duration = round(time.time() - s_time, 6)
    return duration

def sparse_vcfc(path, chrom, start, end):
    assert(os.path.exists(path))
    start = int(start)
    end = int(end)
    cmd = [
        vcfc_cmd,
        path,
        '%s %d %d' % (chrom, start, end)
    ]
    s_time = time.time()
    proc = subprocess.Popen(cmd, stdin=None, stdout=subprocess.DEVNULL)
    proc.wait()
    duration = round(time.time() - s_time, 6)
    return duration


####
# This returns the coordinates for the gene_symbol passed in.
# If not cached in local sqlite db, it queries the Ensembl REST API.
def get_gene_coordinates(gene_symbol):
    def getconn():
        db_path = 'gene_coordinates_hg19.db'
        if not os.path.exists(db_path):
            conn = sqlite3.connect(db_path)
            c = conn.cursor()
            c.execute('create table coordinates(gene_name text, reference_name text, start integer, end integer)')
            conn.commit()
        return sqlite3.connect(db_path)
    def get_from_db(gene_symbol):
        c = getconn().cursor()
        c.execute(
            'select gene_name, reference_name, start, end from coordinates where gene_name = ?',
            (gene_symbol,))
        rs = c.fetchall()
        if len(rs) > 1:
            logging.warn('More than one db record for %s' % gene_symbol)
        elif len(rs) == 0:
            return None
        return {
            'gene_name': rs[0][0],
            'reference_name': rs[0][1],
            'start': rs[0][2],
            'end': rs[0][3]
        }
    def store_in_db(gene_symbol, reference_name, start, end):
        logging.debug('Storing %s in db' % gene_symbol)
        conn = getconn()
        c = conn.cursor()
        c.execute(
            'insert into coordinates(gene_name, reference_name, start, end) values (?,?,?,?)',
            (gene_symbol, reference_name, start, end))
        conn.commit()
    def fetch_from_ensembl(gene_symbol):
        # HowTo:
        # https://grch37.rest.ensembl.org/documentation/info/symbol_lookup
        url = 'https://grch37.rest.ensembl.org/lookup/symbol/human/%s' % gene_symbol
        url += '?content-type=application/json'
        resp = requests.get(url)
        content = resp.content.decode('utf-8')
        if resp.status_code != 200:
            logging.error('Error looking up gene %s: %s' % gene_symbol, content)
            return False
        j = json.loads(content)
        ref = j['seq_region_name']
        start = j['start']
        end = j['end']
        store_in_db(gene_symbol, ref, start, end)
        return True

    coords = get_from_db(gene_symbol)
    if coords is None:
        if not fetch_from_ensembl(gene_symbol):
            raise RuntimeError('Failed to lookup gene symbol: %s' % gene_symbol)
        coords = get_from_db(gene_symbol)
    return coords

def bzgip_find_first_startpos_atleast(bgzip_file_path, min_start_position):
    pass

def vcfc_find_first_startpos_atleast(vcfc_file_path, min_start_position):
    pass

def main():
    coords = [
        get_gene_coordinates('TP53'),
        get_gene_coordinates('ALDH2'),
        get_gene_coordinates('BRCA1'),
        get_gene_coordinates('BRCA2'),

        # https://www.nature.com/articles/d41586-017-07291-9
        get_gene_coordinates('TNF'),
        get_gene_coordinates('EGFR'),
        get_gene_coordinates('VEGFA'),
        get_gene_coordinates('APOE'),
        get_gene_coordinates('IL6'),
        get_gene_coordinates('TGFB1'),
        get_gene_coordinates('MTHFR'),
        get_gene_coordinates('ESR1'),
        get_gene_coordinates('AKT1')
    ]

    iterations = 3
    gene_benchmarks = {}

    def add_benchmark(gene_name, benchmark_type, benchmark_time):
        if gene_name not in gene_benchmarks:
            gene_benchmarks[gene_name] = {}
        print('setting %s benchmark for %s to %s' % (gene_name, benchmark_type, benchmark_time))
        gene_benchmarks[gene_name][benchmark_type] = benchmark_time

    def ensure_file_in_workdir(path, workdir):
        assert(os.path.exists(path))
        if not basename_in_dir(path, workdir):
            shutil.copy(path, workdir)

    # BGZIP + TBI location benchmark
    print('Benchmark: BGZIP + TBI')
    for coord in coords:
        chromosome = coord['reference_name']
        file_path = get_bgzip_file_path(chromosome)
        start = coord['start']
        end = coord['end']
        duration = 0.0
        for _ in range(0, iterations):
            duration += tabix(file_path, chromosome, start, end)
        avg_duration = duration / iterations
        add_benchmark(coord['gene_name'], 'BGZIP+TBI', avg_duration)
        print('gene: %6s, chromosome: %2s, avg_duration: %f' % (
            coord['gene_name'], chromosome, avg_duration))

    # BCF + CSI benchmark
    print('Benchmark: BCF + CSI')
    for coord in coords:
        chromosome = coord['reference_name']
        file_path = get_bcf_file_path(chromosome)
        start = coord['start']
        end = coord['end']
        duration = 0.0
        for _ in range(0, iterations):
            duration += tabix(file_path, chromosome, start, end)
        avg_duration = duration / iterations
        add_benchmark(coord['gene_name'], 'BCF+CSI', avg_duration)
        print('gene: %6s, chromosome: %2s, avg_duration: %f' % (
            coord['gene_name'], chromosome, avg_duration))

    # Sparse VCFC benchmark EXT4
    print('Benchmark: Sparse VCFC EXT4 (block size=4096)')
    for coord in coords:
        chromosome = coord['reference_name']
        file_path = get_bcf_file_path(chromosome)
        start = coord['start']
        end = coord['end']
        duration = 0.0
        for _ in range(0, iterations):
            duration += sparse_vcfc(file_path, chromosome, start, end)
        avg_duration = duration / iterations
        add_benchmark(coord['gene_name'], 'Sparse EXT4(4096)', avg_duration)
        print('gene: %6s, chromosome: %2s, avg_duration: %f' % (
            coord['gene_name'], chromosome, avg_duration))

    print(gene_benchmarks)

    show_graph = True
    if show_graph:
        x = []
        bgzip_y = []
        bcf_y = []
        sparse_vcfc_ext4_4096 = []

        for gene_name in gene_benchmarks:
            benchmark = gene_benchmarks[gene_name]
            print('%8s, BGZIP+TBI: %.10f, BCF+CSI: %.10f' % (
                gene_name, benchmark['BGZIP+TBI'], benchmark['BCF+CSI']
            ))
            x.append(gene_name)
            bgzip_y.append(benchmark['BGZIP+TBI'])
            bcf_y.append(benchmark['BCF+CSI'])
            sparse_vcfc_ext4_4096.append(benchmark['Sparse EXT4(4096)'])

        plt.xlabel('Scheme')
        plt.ylabel('Time (seconds)')
        plt.show()

if __name__ == '__main__':
    exit(main())