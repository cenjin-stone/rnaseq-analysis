#!/bin/bash
#only ch1 index
cd /home/sxy/project
mkdir -p reference/star_index
STAR --runMode genomeGenerate \
     --genomeDir reference/star_index \
     --genomeFastaFiles reference/ch1.fa \
     --sjdbGTFfile reference/homo.gtf \
     --sjdbOverhang 74 \
     --genomeSAindexNbases 12 \
     --runThreadN 4 2>&1 | tee logs/star_index.log
