import time, os
import subprocess
import re
from evaluation.config import Config

#cmd = './main_release sparse-query /mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.100k.vcf.vcfc.sparse 22 16050075 19757157'
#cmd = './tabix /mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 22 16050075 19757157'

def flush_cache():
    # Relies on passwordless sudo
    cmd_args = ['sudo', '/usr/local/sbin/flush-cache']
    proc = subprocess.Popen(cmd_args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = proc.communicate()
    if len(stderr) > 0:
        raise RuntimeError(stderr.decode('utf-8'))

def time_cmd(cmd_args, do_flush_cache=True) -> float:
    if do_flush_cache:
        flush_cache()
    print(cmd_args)
    start_time = time.time()
    # TODO see if capturing stderr in PIPE has any measurable impact. Should only
    # store content on error, in which case we don't care as much about accurate time anyways.
    proc = subprocess.Popen(cmd_args, stdin=None, stdout=subprocess.DEVNULL)
    proc.wait()
    end_time = time.time()
    duration = round(end_time - start_time, 6)
    if proc.returncode != 0:
        raise RuntimeError('cmd: %s failed with status %d' % (cmd_args, proc.returncode))
    return duration

def create_tabix_index(config:Config, filename:str) -> float:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    cmd_args = [
        config.tabix_command_path,
        '-f',
        filename
    ]
    return time_cmd(cmd_args)


def run_tabix(config:Config, filename:str, ref:str, start:int, end:int) -> float:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    cmd_args = [
        config.tabix_command_path,
        filename,
        '%s:%d-%d' % (ref, start, end)
    ]
    return time_cmd(cmd_args)


def create_vcfc_sparse_file(config:Config, filename:str) -> float:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    cmd_args = [
        config.get_vcfc_release_cmd(),

    ]


def run_vcfc_sparse_query(config:Config, filename:str, ref:str, start:int, end:int) -> float:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    cmd_args = [
        config.get_vcfc_release_cmd(),
        'sparse-query',
        filename,
        '%s:%d-%d' % (ref, start, end)
    ]
    return time_cmd(cmd_args)


def create_binned_index(config:Config, filename:str, bin_size:int) -> float:
    cmd_args = [
        config.get_vcfc_release_cmd(),
        'create-binned-index',
        str(bin_size),
        filename
    ]
    return time_cmd(cmd_args)

def run_vcfc_binned_index_query(config:Config, filename:str, ref:str, start:int, end:int) -> float:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    if not os.path.exists(filename + '.vcfci'):
        raise RuntimeError('%s does not exist' % (filename + '.vcfci'))

    cmd_args = [
        config.get_vcfc_release_cmd(),
        'query-binned-index',
        filename,
        '%s:%d-%d' % (ref, start, end)
    ]
    return time_cmd(cmd_args)


def create_vcfc_sparse_external_index(config:Config, filename:str) -> float:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    cmd_args = [
        config.get_vcfc_release_cmd(),
        'create-sparse-index',
        filename,
        # '%s:%d-%d' % (ref, start, end)
    ]
    return time_cmd(cmd_args)

def run_vcfc_sparse_external_index_query(config:Config, filename:str, ref:str, start:int, end:int) -> float:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    if not os.path.exists(filename + '.vcfci-sparse'):
        raise RuntimeError('%s does not exist' % (filename + '.vcfci-sparse'))

    cmd_args = [
        config.get_vcfc_release_cmd(),
        'query-sparse-index',
        filename,
        '%s:%d-%d' % (ref, start, end)
    ]
    return time_cmd(cmd_args)


def construct_timing_profile(output:list, convert_to_seconds=True) -> dict:
    '''
    Expects `output` to be a byte array, like that returned from Popen.communicate
    '''

    output = output.decode('utf-8')
    profile = {}
    pattern = re.compile(r'TIMING (\w+): (\d+)')
    for line in output.split('\n'):
        match = pattern.match(line)
        if match:
            component = match.group(1)
            nanoseconds = int(match.group(2))
            if component not in profile:
                profile[component] = 0
            profile[component] += nanoseconds

    if convert_to_seconds:
        for label in profile:
            profile[label] = profile[label] / 1e9

    return profile


def run_vcfc_binned_index_timing_profile(config:Config, filename:str, ref:str, start:int, end:int) -> dict:
    if not os.path.exists(filename):
        raise RuntimeError('%s does not exist' % filename)
    if not os.path.exists(filename + '.vcfci'):
        raise RuntimeError('%s does not exist' % (filename + '.vcfci'))

    cmd_args = [
        config.get_vcfc_timing_cmd(),
        'query-binned-index',
        filename,
        '%s:%d-%d' % (ref, start, end)
    ]

    # Load file into host cache
    # os.system('cat %s > /dev/null' % filename)

    # flush_cache()

    start_time = time.time()
    proc = subprocess.Popen(cmd_args, stdin=None, stdout=subprocess.PIPE)
    # proc.wait() # cannot use this with stdout=subprocess.PIPE
    stdout, stderr = proc.communicate()
    end_time = time.time()
    duration = round(end_time - start_time, 6)

    # print('Constructing profile')
    profile = construct_timing_profile(stdout)
    return profile
