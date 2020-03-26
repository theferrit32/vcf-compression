import time
import subprocess

from .config import Config

#cmd = './main_release sparse-query /mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.100k.vcf.vcfc.sparse 22 16050075 19757157'
#cmd = './tabix /mnt/ext4/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz 22 16050075 19757157'

def time_cmd(cmd_args) -> float:
    start_time = time.time()
    proc = subprocess.Popen(cmd_args, stdin=None, stdout=subprocess.PIPE)
    proc.wait()
    end_time = time.time()
    duration = round(end_time - start_time, 6)
    return duration

def run_tabix(config:Config, ref:str, start:int, end:int) -> float:
    cmd_args = [
        config.tabix_command_path,
        ref, str(start), str(end)
    ]
    return time_cmd(cmd_args)

def run_vcfc_sparse_query(config:Config, filename:str, ref:str, start:int, end:int) -> float:
    cmd_args = [
        config.vcfc_command_path,
        'sparse-query',
        filename,
        ref, str(start), str(end)
    ]
    return time_cmd(cmd_args)
