#!/bin/bash
cd /home/sxy/project/data/raw

fastq-dump --split-files --gzip SRR25386819 &
fastq-dump --split-files --gzip SRR25386815 &
fastq-dump --split-files --gzip SRR25386807 &
fastq-dump --split-files --gzip SRR25386803 &

jobs

wait

echo "download is completed"
