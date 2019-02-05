#!/bin/bash

get_abs_parentdir() {
  # $1 : relative filename
  filename=$1
  parentdir=$(dirname "${filename}")
  echo "$(cd "${parentdir}" && pwd)"
}

test_bin=$1

if [[ ! -f $test_bin ]]; then
  echo "$#"
  echo "test file \"$test_bin\" not found!!!"
  exit 1
fi

ref_genome=./data/test.mmi
fastq1=./data/test_1.fastq.gz
fastq2=./data/test_2.fastq.gz

echo "GLOG_v=3 GLOG_logtostderr=1 $test_bin $ref_genome $fastq1 $fastq2"

GLOG_v=3 GLOG_logtostderr=1 $test_bin $ref_genome $fastq1 $fastq2