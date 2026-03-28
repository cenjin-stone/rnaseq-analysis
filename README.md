# RNA-seq & scRNA-seq 分析项目

## 项目一：Bulk RNA-seq 差异表达分析
- **数据**：人类肺癌细胞系 H520，4 样本（2 对照 + 2 处理），GSE237935
- **流程**：fastp → STAR → featureCounts → DESeq2
- **结果**：30 个差异基因（5 上调，25 下调），火山图、热图
- **说明**：本地内存限制，仅使用 1 号染色体验证流程

## 项目二：单细胞转录组分析（PBMC 5k）
- **数据**：10X Genomics 官方 PBMC 5k 数据集
- **流程**：Seurat 标准流程（QC、降维、聚类、注释）
- **结果**：鉴定出 13 个细胞群（CD4+ T、CD8+ T、B 细胞、单核细胞、NK 细胞、树突状细胞、血小板等）
- **图**：UMAP、热图、点图、小提琴图、FeaturePlot

## 环境
- Linux / R / Seurat / DESeq2 / STAR / fastp

## 作者
石新宇
