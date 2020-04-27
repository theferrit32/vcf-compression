#!/bin/bash
set -xe

date
start=$(date +%s)
if [[ -z "$1" ]]; then
    echo "Must provide a node name, like 'c[80]'"
    exit 1
fi

if [[ -z "$2" ]]; then
    echo "Must provide a filesystem name, like ext4 or xfs"
    exit 1
fi

eval_root=~/vcf-compression/evaluation
node="$1"
fs="$2"

function queue {
    if [[ -z "$1" ]]; then
        echo "Must provide a test name"
        exit 1
    fi
    sbatch -p max -w "${node}" bash ${eval_root}/timing.sh "$1" "${fs}"
}
function run {
    if [[ -z "$1" ]]; then
        echo "Must provide a test name"
        exit 1
    fi
    bash ${eval_root}/timing.sh "$1" "${fs}"
}
op=run

$op all-exhaustive-single
$op all-exhaustive-range

$op binned-timing-profile-single
$op binned-timing-profile-range

$op binned-index-creation-time
$op all-indexing-times

date
end=$(date +%s)
dur=$((end-start))
echo "Finished all runs in ${dur} seconds"
