#!/bin/bash
set -ex
#python3 ~/vcf-compression/evaluation/range_query.py measure-binned-timing-profile-range

evaluation_root=~/vcf-compression/evaluation

function run_measure {
  if [[ -z "$1" ]]; then
    echo "Must pass a test name"
    exit 1
  fi
  python3 -u ${evaluation_root}/evaluation_main.py "$1" measure | tee "${1}.log"
}

if [[ -z "$1" ]]; then
  echo "Must pass a test name"
  exit 1
fi

echo "Running $1"

run_measure "$1"

# run_measure sparse-exhaustive-range
# run_measure sparse-exhaustive-single
# run_measure binned-exhaustive-single
# run_measure binned-timing-profile
# run_measure binned-timing-profile-range

echo "Finished run"
