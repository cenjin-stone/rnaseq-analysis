#!/bin/bash

cd /home/sxy/project/data/raw

declare -A name
name[SRR25386819]=control1
name[SRR25386815]=control2
name[SRR25386807]=treat1
name[SRR25386803]=treat2

for old in "${!name[@]}";do
	new="${name[$old]}"
	mv ${old}_1.fastq.gz ${new}_R1.fastq.gz
	mv ${old}_2.fastq.gz ${new}_R2.fastq.gz
done
echo "$old --- $new"
