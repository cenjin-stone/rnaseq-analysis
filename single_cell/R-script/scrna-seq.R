#Seurat 模拟数据流程
#数据来源：10X Genomics
rm(list=ls())
setwd("F:/rnaseq-analysis/singleCell/PBMC/matrix")
library(Seurat)
library(tidyverse)
###1.数据读取----
pbmc.data <- Read10X(data.dir = "F:/rnaseq-analysis/singleCell/PBMC/matrix")
#直接读取路径下
dim(pbmc.data)#查看基因数量；细胞数量
pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc5k",min.cells = 3,min.features = 200)
#参数：counts:数据集；project:项目名称自己取；mincell:1个基因至少3个细胞有表达才保留，minifeature一个细胞至少检测到200个基因
#图1
hist(pbmc$nFeature_RNA, breaks = 100, main = "Genes per cell", xlab = "Number of genes")

###2.查看线粒体基因占比----
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
#用百分比特征函数，找到MT开头的基因，放在自己命名的组下（线粒体基因占比高说明细胞存在各种问题）
VlnPlot(pbmc, features = "percent.mt")
pbmc <- subset(pbmc, subset = percent.mt < 10)#挑选这个组小于10的留下
pbmc#查看下剩余基因与细胞数目

###3.细胞标准化----
pbmc <- NormalizeData(pbmc)
pbmc#查看一下

###4.找高变基因----
pbmc <- FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 2000)#选择好方法，筛选出2000的高变基因
top10 <- VariableFeatures(pbmc)[1:10]
p1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = p1, points = top10, repel = TRUE)#添加标签的图p1，添加点
#高变基因图----
plot2

###5.基因PCA分析----
pbmc <- ScaleData(pbmc)#标准化
pbmc <- RunPCA(pbmc)#PCA分析
#观察下pc拐点
ElbowPlot(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)#找到细胞的“邻居”
pbmc <- FindClusters(pbmc, resolution = 0.5)#把相似细胞合成一个集群，值越小分的群就越好，越大分的群就越多越细致
pbmc <- RunUMAP(pbmc, dims = 1:10)#高纬数据转变为二维
#UMAP图----
DimPlot(pbmc, label = TRUE)
DimHeatmap(pbmc, dims = 1:6, cells = 500, balanced = TRUE)

###6.找maker基因----某个细胞高表达，其他细胞低表达，是一种身份象征的基因，
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)#只找positive基因，至少要25%的细胞中表达，
head(pbmc.markers)
DefaultAssay(pbmc)

###6.1选取某一个cluster，例子是0----
cd4_markers <- FindMarkers(pbmc, ident.1 = "CD4+ T", min.pct = 0.25, logfc.threshold = 0.25)#这个因为是后面取名字了所以用这个。不然用数字编号
head(cd4_markers)
cd4_up <- cd4_markers[cd4_markers$p_val_adj<0.05 & cd4_markers$avg_log2FC>0.5, ]
VlnPlot(pbmc, features = rownames(head(cd4_up, 6)), ncol = 3)#feature:直接填基因，c（），每行3个图
FeaturePlot(pbmc, features = rownames(head(cd4_up, 4)))#Umap分布图

###7.找到每个集群的top5基因名称
top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)#把前面得到基因，按照cluster分组，在排序函数，值最大的前5个，根据：a_l分
top5_list <- split(top5$gene, top5$cluster)#要把gene拆分出来，依据cluster来拆
print(top5_list)

#每个cluster命名
new.cluster.ids <- c(
  "CD4+ T", "Monocyte", "CD4+ T", "CD8+ T", "CD8+ T",
  "CD4+ T", "NK cell", "B cell", "Monocyte", "Proliferating",
  "Dendritic cell", "Platelet", "Platelet")
names(new.cluster.ids) <- levels(pbmc)#给每个名字加上编号，levels(pbmc)的值是代码41行分出来的集群
pbmc <- RenameIdents(pbmc, new.cluster.ids)#重新给cluster命名

###8.Umap----
png("Umap_annotated.png",width = 8, height = 6, units = "in", res = 300)
DimPlot(pbmc,label = TRUE,label.size=4,split.by = "orig.ident")#显示标签
dev.off()  

###9.经典maker分布图----
png("featureplot_classic_markers.png", width = 12, height = 10, units = "in", res = 300)
FeaturePlot(pbmc, features = c("CD4", "CD8A", "MS4A1", "CD14", "NKG7"))
dev.off()

###10.点图----
png("dotplot_markers.png", width = 12, height = 8, units = "in", res = 300)
DotPlot(pbmc, features = unique(top5$gene)) + RotatedAxis()#去重，旋转45度
dev.off()

###11.小提琴图----
png("vlnplot_markers.png", width = 12, height = 8, units = "in", res = 300)
VlnPlot(pbmc, features = c("CCR4", "S100A8", "KLRF1", "MS4A1","FCER1A"))#填写感兴趣的基因
dev.off()

###12.热图
png("heatmap_top5_markers.png", width = 20, height = 10, units = "in", res = 300)
DoHeatmap(pbmc, features = top5$gene)
dev.off()


