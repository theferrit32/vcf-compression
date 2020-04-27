#!/bin/bash
set -xe

date

eval_root=~/vcf-compression/evaluation
# node="$1"

function run {
    if [[ -z "$1" ]]; then
        echo "Must provide a test name"
        exit 1
    fi
    python3 ${eval_root}/evaluation_main.py "$1" graph
}
op=run

$op all-exhaustive-single
$op all-exhaustive-range

$op binned-timing-profile-single
$op binned-timing-profile-range

$op binned-index-creation-time
$op all-indexing-times


echo "Finished all runs"
