rm(list = ls())#环境清空
setwd("F:/生物数据/myproject")

#只做了人类第一条染色体的STAR数据，后续数据不代表整个基因组的真实数据。仅为代码参考
ch1_exp <- read.table("counts_matrix.txt",
                    header = TRUE,#第一行设置为列名
                    row.names = 1)#第一列设置为行名
table(is.na(ch1_exp))#统计下NA
exp[is.na(ch1_exp)] <- 0#将NA值转换成数值0
###1.分组----
symbol_id <- colnames(ch1_exp)
group <- rep(c("C","T"),each=2)
group
rt <- data.frame(row.names = symbol_id,group=group)
identical(rownames(rt),colnames(ch1_exp))#确定rt：行名和exp：列名是否一致

###2.DESeq2分析----
library(tidyverse)
library(DESeq2)
rt$group <- as.factor(rt$group)#变成因子
str(rt$group)
####2.1 构建对象
dds <- DESeqDataSetFromMatrix(countData = ch1_exp,colData = rt,design = ~group)
####2.2 筛选表达量
select_exp <- rowSums(counts(dds))>1
dds <- dds[select_exp,]

###3.差异表达分析----
dds <- DESeq(dds)
levels(dds$group)
res <- results(dds,contrast = c('group','T','C'))
#30行代码格式为：c('因子名'，'实验组'，'对照组')

###4.提取结果----
res1 <- as.data.frame(res)
###5.结果筛选（padj < 0.05 & log2FoldChange >1）可更改
degs <- filter(res1,padj <=0.05 & abs(log2FoldChange)>=1)

###绘图----
# 绘制离散图
plotDispEsts(dds)
# P值直方图，看P值显著基因有多少
hist(res1$padj)

###保存数据----
write.table(degs,file = "ch1_DEG.txt",sep = "\t",row.names = T,col.names = NA,quote = F)#保存行名，第一列NA，纯净数据

###6.找到上下调基因----
DEG <- res1
logFC_cutoff <- 1#做一个分割
type1=(DEG$padj < 0.05)&(DEG$log2FoldChange < -logFC_cutoff)#down
type2=(DEG$padj < 0.05)&(DEG$log2FoldChange > logFC_cutoff)#up
DEG$change <- ifelse(type2,"Up",ifelse(type1,"Down","NOT"))#标出上下调
table(DEG$change)#查看一下

###绘图----
library(ggplot2)
library(cowplot)
library(ggrepel)

DEG <- mutate(DEG,Gene_symbol=rownames(DEG))#mutate函数添加新的列名
UpDEG <- subset(DEG,change=='Up')#数据过少，就这样子吧
DownDEG <- subset(DEG,change=='Down')
DownDEG_5 <- top_n(x=DEG,n=-5,wt=pvalue)#挑选最显著的5个：x数据框；n行数，正数最大，反之最小；wt依据的列

# 基础火山图代码----data=DEG
p <- ggplot(data = DEG, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.5, size = 4.5, aes(color = change)) +
  ylab('-log10(padj)') +
  scale_color_manual(values = c('#1F7784', 'grey', '#FF7F0E')) +
  
  # 添加阈值线
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue", linewidth = 0.8) +
  
  # 调整坐标轴
  scale_x_continuous(breaks = seq(floor(min(DEG$log2FoldChange)), ceiling(max(DEG$log2FoldChange)), by = 1)) +
  
  # 添加标题和图例
  labs(title = "Volcano Plot of Differential Expression",
       subtitle = paste0("Total genes: ", nrow(DEG), 
                         " | Up: ", sum(DEG$change == "Up"，na.rm = TRUE), 
                         " | Down: ", sum(DEG$change == "Down"na.rm = TRUE)),
       x = "log2 Fold Change",
       color = "Expression") +
  
  # 主题设置
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank()
  )

# 标记显著基因（可选）
if(exists("UpDEG") && exists("DownDEG_5")){
  top_genes <- rbind(UpDEG, DownDEG_5)
  p <- p + 
    geom_point(data = top_genes, aes(x = log2FoldChange, y = -log10(padj)), 
               color = "black", size = 3, shape = 1) +
    ggrepel::geom_text_repel(
      data = top_genes,
      aes(label = Gene_symbol),
      size = 4,
      box.padding = 0.5,
      max.overlaps = Inf
    )
}
print(p)
ggsave("Volcano_Plot.png", plot = p, width = 10, height = 8, dpi = 300)

###热图----data=ch1_exp
library(pheatmap)
identical(colnames(ch1_exp),rownames(rt))
diff <- rownames(DEG)[DEG$change !='NOT']#找出所有非NOT的行
exp_diff <- ch1_exp[diff,]#找到diff所有行

diff <- DEG[DEG$change != 'NOT', ]
Up <- diff %>% arrange(desc(log2FoldChange)) %>%head(5)#从大到小，函数desc
Down <- diff %>% arrange(log2FoldChange)%>% head(5)#从小到大

all_genes <- c(rownames(Up), rownames(Down))

exp_diff <- ch1_exp[all_genes, ]

heatmap_colors <- colorRampPalette(c('blue', 'white', 'red'))(100)
pheatmap(exp_diff,
         annotation_col = rt,  # 列注释（样本分组）
         color = heatmap_colors,  # 颜色梯度
         scale = "row",  # 按行标准化（Z-score标准化）
         cluster_cols = TRUE,  # 对列（样本）进行聚类
         cluster_rows = TRUE,  # 对行（基因）进行聚类
         show_rownames = TRUE,  # 不显示所有基因名（通常太多）
         show_colnames = TRUE,  # 显示样本名
         fontsize_row = 8,  # 行名字体大小
         fontsize_col = 10,  # 列名字体大小
         border_color = NA,  # 无边框颜色
         main = "Gene Expression Heatmap ",  # 标题
         annotation_colors = list(  # 自定义注释颜色
           Group = c(Control = "grey", Treatment = "red")  # 假设rt中有Group列
         ),
         clustering_distance_rows = "euclidean",  # 行聚类距离
         clustering_distance_cols = "euclidean",  # 列聚类距离
         clustering_method = "complete"  # 聚类方法
)




