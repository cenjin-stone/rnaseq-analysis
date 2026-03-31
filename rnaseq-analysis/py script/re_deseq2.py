import pandas as pd     # 读数据
import numpy as np      # 数学计算
import matplotlib.pyplot as plt   # 画图
#from matplotlib.patches import Ellipse
#import matplotlib.transforms as transforms
import seaborn as sbs   # 
from sklearn.decomposition import PCA #pca

print("all pakege is already")

##1.读入
df=pd.read_csv("counts_matrix.txt",sep="\t",index_col=0)#index=0第一列当作行名
print(df.shape)#.shape输出多少行（0），多少列（1）
##清洗
df_clean=df[df.sum(axis=1)>10]#按照列的方向求和
print(df_clean.shape)

#PCA
#df_pca=df_clean.T#转置
df_pca = np.log2(df_clean.T + 1)#用工具压缩对数，+1防报错
print(df_pca.head(5))
#从PCA工具降维
pca=PCA(n_components=2)#降维2个维度
pca_result = pca.fit_transform(df_pca)
print(pca_result)#数据
#pca图
plt.figure(figsize=(6,5))#画布创建6*5
#分散点---定位结果中第一个样品的坐标x=0，y=1，颜色，点大小，图例名
group=[
    {"key":[0,1],
     "color":"red",
     "label":"Control"},
    {"key":[2,3],
     "color":"blue",
     "label":"Treat"}
]
for g in group:
    for k in g["key"]:#一个点循环
        if k==g["key"][0]:
            plt.scatter(pca_result[k,0],pca_result[k,1],c=g["color"],s=120,label=g["label"])
        else:
             plt.scatter(pca_result[k,0],pca_result[k,1],c=g["color"],s=120)

plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%})")#x轴计算
plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%})")#里面()为了好看
#f运行文字变量，0--->plt.xlabel
plt.title("PCA Plot")
plt.legend()#图例显示
plt.show()

#数据in
DEG=pd.read_csv("ch1_res.txt",index_col=0,sep="\t")
DEG=DEG[DEG["padj"].notna()].copy()#clean
print(DEG.head(5))

#火山：Y=-log10(padj)
DEG["-log10(padj)"]=-np.log10(DEG["padj"])#添加列
DEG["group"]="Not"#全标记为not
#搜索器.loc（行）
mask_up = (DEG["padj"] < 0.05) & (DEG["log2FoldChange"] > 1)
DEG.loc[mask_up, "group"] = "Up"
mask_down = (DEG["padj"] < 0.05) & (DEG["log2FoldChange"] < -1)
DEG.loc[mask_down, "group"] = "Down"
print(DEG["group"].value_counts())#频数统计

deg_all=DEG[DEG["group"]!="Not"].copy()
#火山图
plt.figure(figsize=(6,5),dpi=100)
colors={
    "Up":"red",
    "Down":"blue",
    "Not":"gray"
}
for g in colors.keys():
    subset=DEG[DEG["group"]==g]
    #plt.scatter([],[],c=colors[g],label=g)
    #plt.scatter(DEG["log2FoldChange"],DEG["-log10(padj)"])#x轴，Y轴
    plt.scatter(
        subset["log2FoldChange"],subset["-log10(padj)"],
        c=colors[g],
        s=50,#point size
        alpha=0.7,#透明度
        label=g
    )
plt.title("Differential gene volcano plot")
plt.legend()
plt.xlabel("log2FoldChange")
plt.ylabel("-log10(padj)")
plt.grid(alpha=0.3)
plt.show()

