#!/bin/bash
cd /home/sxy/project/results/counts

group=(control1 control2 treat1 treat2)
name="${group[0]}_counts.txt"

for file in "${group[@]:1}";do
	name="${name} ${file}_counts.txt"
done	
echo "$name"

#paste $name|head -2 |tail -1

paste $name|awk 'BEGIN{OFS="\t"} NR==2{print "Gene_id","control1","control2","treat1","treat2"} NR>2{print $1,$7,$14,$21,$28}' > counts_matrix.txt

cat counts_matrix.txt |head -3 


