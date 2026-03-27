rm(list = ls())#环境清空
setwd("F:/生物数据/myproject")
Ecoli <- read.table("Ecoli_counts_CDS.txt",
                    header = TRUE,#第一行设置为列名
                    skip = 1,#跳过第一行
                    row.names = 1)#第一列设置为行名
colnames(Ecoli)
colnames(Ecoli)[6] <- "Counts"
exp <- Ecoli[,"Counts",drop=FALSE]#选择Counts列，drop=F保持数据框格式
top_gene <- exp[order(-exp$Counts),,drop=FALSE]#按照行排序，order函数是排序
head(top_gene,10)

top10_genes <- rownames(head(exp[order(-exp$Counts), , drop = FALSE], 10))#获取前10基因的名称
print(top10_genes)
summary(exp$Counts)
#直方图
hist(exp$Counts,
     breaks = 50,
     main = "Gene Expression Distribution",
     xlab = "Counts",
     col = "steelblue")#break 柱子数量

high_exp <- exp[exp$Counts > 5000, , drop = FALSE]
cat("表达量 > 5000 的基因数:", nrow(high_exp), "\n")
low_exp <- exp[exp$Counts < 100, , drop = FALSE]
cat("表达量 < 100 的基因数:", nrow(low_exp),"\n")
zero_exp <- exp[exp$Counts == 0, , drop = FALSE]
cat("表达量为 0 的基因数:", nrow(zero_exp), "\n")
#箱线图
boxplot(exp$Counts, main = "Gene Expression Distribution", 
        ylab = "Counts", col = "lightblue")
#密度图
plot(density(exp$Counts), main = "Density Plot of Expression", 
     xlab = "Counts", col = "blue")
#数据保存
write.csv(exp,"Ecoli_all_gene.csv")
write.csv(high_exp,"Ecoli_high_exp.csv")
