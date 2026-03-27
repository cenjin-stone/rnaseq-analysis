#!/bin/bash

cd /home/sxy/project

samples=(control1 control2 treat1 treat2)

for group in "${samples[@]}";do
	featureCounts -T 4 \
	-a reference/homo.gtf \
	-o results/counts/${group}_counts.txt \
	-t exon \
	-g gene_id \
	-p \
	results/aligned/${group}_Aligned.sortedByCoord.out.bam 2>&1 |tee logs/feature_${group}.log
	echo "${group} is completed"
done

echo "all is completed"
