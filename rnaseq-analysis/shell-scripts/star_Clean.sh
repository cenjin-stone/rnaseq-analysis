#!/bin/bash
cd /home/sxy/project
clean_date=(control1 control2 treat1 treat2)

for i in "${clean_date[@]}";do
	STAR	--genomeDir reference/star_index \
		--readFilesIn data/clean/${i}_R1_clean.fastq.gz data/clean/${i}_R2_clean.fastq.gz \
		--readFilesCommand zcat \
		--outFileNamePrefix results/aligned/${i}_ \
		--outSAMtype BAM SortedByCoordinate \
		--runThreadN 4 2>&1 | tee logs/star_${i}.log
done
echo "all is complete"
