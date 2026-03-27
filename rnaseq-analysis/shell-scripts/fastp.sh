#!/bin/bash

cd /home/sxy/project

sample_raw=("control1" "control2" "treat1" "treat2")

for sample in "${sample_raw[@]}";do
	fastp	-i data/raw/${sample}_R1.fastq.gz \
		-I data/raw/${sample}_R2.fastq.gz \
                 -o data/clean/${sample}_R1_clean.fastq.gz \
        	 -O data/clean/${sample}_R2_clean.fastq.gz \
		--detect_adapter_for_pe \
		--cut_front --cut_tail \
		--cut_window_size 4 \
		--cut_mean_quality 20 \
		--length_required 40 \
		 --n_base_limit 5 \
         	 --average_qual 20 \
        	  --thread 4 \
        	  -h results/qc/${sample}_report.html \
        	  -j results/qc/${sample}_report.json \
        	  2>&1 | tee logs/${sample}_fastp.log
    
    echo "$sample is completed"
done
