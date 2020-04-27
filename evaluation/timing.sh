#!/bin/bash
set -ex
#python3 ~/vcf-compression/evaluation/range_query.py measure-binned-timing-profile-range

evaluation_root=~/vcf-compression/evaluation

function run_measure {
  if [[ -z "$1" ]]; then
    echo "Must pass a test name"
    exit 1
  fi
  if [[ -z "$2" ]]; then
    echo "Must pass a filesystem name (ext4|xfs)"
    exit 1
  fi
  python3 -u ${evaluation_root}/evaluation_main.py "$1" measure --filesystem "$2" | tee "${1}.log"
}

if [[ -z "$1" ]]; then
  echo "Must pass a test name"
  exit 1
fi

if [[ -z "$2" ]]; then
  echo "Must pass a filesystem name (ext4|xfs)"
  exit 1
fi

testname="$1"
fs="$2"

echo "Running $testname"

run_measure "$testname" "$fs"

# run_measure sparse-exhaustive-range "$fs"
# run_measure sparse-exhaustive-single "$fs"
# run_measure binned-exhaustive-single "$fs"
# run_measure binned-timing-profile "$fs"
# run_measure binned-timing-profile-range "$fs"

echo "Finished run"
