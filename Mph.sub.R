library(qs)
library(Seurat)
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
#library(SingleR)
#library(CCA)
library(clustree)
library(cowplot)
library(monocle)
library(tidyverse)
library(SCpubr)
library(UCell)
library(irGSEA)
library(GSVA)
library(GSEABase)
library(harmony)
library(plyr)
library(randomcoloR)
library(CellChat)
library(cowplot)
library(dplyr)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
library(scCustomize)
library(qs)
library(scRNAtoolVis)
library(jjAnno)
library(plotly)
library(SCP)


.libPaths(c("D:\\R_tools_install\\R-4.4.2\\Seuratobject_V5", .libPaths()))
library(SeuratObject)
packageVersion("SeuratObject")
library(SCP)



setwd("F:\\Sprouts\\CRC.wsl\\3.Mph.sub\\")


sce  <- qread("Mph.qs")

Idents(sce) <- sce@meta.data$celltype
# Fib <- subset(sce,idents=c("Mye"))
# 
# sce = Fib
############################################################################
###                            1.处理数据过滤                           ###
###########################################################################
sce  <- NormalizeData(object = sce , 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
sce<- FindVariableFeatures(object = sce , 
                           selection.method = "vst", nfeatures = 3000)
scale.genes <-  VariableFeatures(sce)
sce  <- ScaleData(sce, features = scale.genes)
#YDL <- ScaleData(object = YDL, features = rownames(YDL))对全部基因标准化
#线性降维（PCA）,默认用高变基因集,但也可通过features参数自己指定；

sce =RunPCA(object= sce,
            npcs = 10,pc.genes=VariableFeatures)     #PCA分析
##整合方法1，单个样本间进行整合
##max.iter.cluster默认是10，决定运行RunHarmony结果再迭代多少次后完成了收敛
#lamba参数，默认是1，决定整合力度，值越小，整合力度越大
ElbowPlot(sce , ndims=50, reduction="pca") 
pc.num=1:4

############################################################################
###                               1.3数据   harmony                      ### 
###########################################################################
sce   <- sce  %>%
RunHarmony("Patient", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sce, 'harmony')

############################################################################
###                               1.4数据   选定分辨率                  ### 
###########################################################################
options(future.globals.maxSize = 512*1000 * 1024^2)

sce=FindNeighbors(sce, dims = pc.num, reduction = "harmony")
for (res in c(0.01, 0.02,0.03,0.05, 0.1, 0.2,0.3, 0.5,1.5,2.0,2.5,3.0)) {
  sce=FindClusters(sce, graph.name = "RNA_snn", 
                   resolution = res, algorithm = 1)}
apply(sce@meta.data[,grep("RNA_snn_res",colnames(sce@meta.data))],2,table)

#p2_tree=clustree(sce@meta.data, prefix = "RNA_snn_res.")
############################################################################
###                               1.5数据降维                           ### 
###########################################################################
#min_dist -- 是低维空间中点之间的最小距离。这个参数控制了UMAP将点紧密地聚集在一起的程度，数值越低，嵌入的点就越紧密。
#更大的min_dist值会使UMAP将点更松散地包装在一起，而专注于保存广泛的拓扑结构。
sce <- sce  %>%
  RunUMAP(reduction = "harmony",
          dims = pc.num,min.dist = 0.4) %>%
  #RunTSNE(reduction = "harmony", dims = pc.num) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution = 0.2) %>%
  identity()
############################################################################
###                               1.5Topgene                           ### 
##########################################################################
Idents(sce) <- sce@meta.data$RNA_snn_res.0.2
sce.markers <- FindAllMarkers(sce, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25)

write.csv(sce.markers,file = 
            "sce.ALL_markers.csv")
write.csv(top20af.markers,file = 
            "top20af.markers.markers.csv")
head(sce.markers)
# get top 10 genes
top20af.markers <- sce.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)


#写出GMT文件
df=as.data.frame(sce.markers[,c("gene","cluster")])
table(df$cluster)
dirnb=(paste0(getwd(),"/D1.afgmt.gmt"))
gmtlist=list()
for (afk in as.character(names(table(df$cluster)))) {
  ddgene=rownames(df)[df$cluster==afk]
  aaname=afk
  gmtlist[[aaname]]=ddgene
}
output_gmt<-function(geneset,file){
  sink(file) 
  lapply(names(geneset),function(i){
    cat(paste(c(i,'NA',geneset[[i]]),collapse='\t')) 
    cat('\n') 
  })
  sink() 
}
output_gmt(gmtlist,dirnb)


Idents(sce) <- sce$celltype
CAFC3 <- subset(sce,ident="15")
cell_names <- colnames(CAFC3)
sce=SetIdent(object = sce,cells = cell_names,value = 'TIMP13_Mph')
sce@meta.data$celltype<- sce@active.ident#
table(sce$celltype)
############################################################################
###                               1.6 看一下效果                         ### 
###########################################################################
library(ggsci)
mycolors = c(pal_d3("category20")(25),
             pal_d3("category20b")(25),
             pal_d3("category20c")(25),
             pal_d3("category20")(20))
mycolors = c("#9ab882","#584f8a","#b02c3e","#d9792b","#C8E3C1","#A4C9DE","#F6CBDE","#A05528","#9bbaca","#b6acc2","#387091")
DimPlot(sce, reduction = "tsne",label = T,
        raster=FALSE, group.by = "celltype",pt.size = 0.05,
        cols = mycolors)+ scale_colour_manual(values = mycolors)


DimPlot(sce, reduction = "umap",label = T,
        raster=FALSE, group.by = "celltype",pt.size = 0.01,
        cols = mycolors)

DimPlot(sce, reduction = "umap",label = T,
        raster=FALSE, group.by = "RNA_snn_res.3",pt.size = 0.1,
        cols = mycolors)+ scale_colour_manual(values = mycolors)




############################################################################
###                         1.6 看一下marker检验                         ### 
###########################################################################



genes_to_check = c(
  "CCL5","NKG7", "CD3D", "CD3E", "CD8A",#0：T cell:
  "PECAM1","KRT19","CLDN4","GKN1","PGC",#1:Epi:
  "IGHG1","CD38",'MZB1', 'SDC1',"JCHAIN","SLAMF7",#2:plasma:
  "C1QC","CSF1R","CD68","LYZ","S100A8",#3:Macro:
  'LUM' , 'FGF7',"DCN","COL1A1","CXCL14",#4:Fib:
  "CD79A",'MS4A1','CD19',"BANK1","BLK",#5:B/plasma:
  "CDH5","VWF","FABP4","PLVAP",#6:Endo:
  "RGS5","TAGLN","ACTA2","MYL9","IGFBP7",#7:Pericyte:
  'CPA3' ,'CST3', 'KIT', 'TPSAB1','TPSB2',#8:Mast cell
  "MKI67","STMN1","PCNA"#epi
)

markers.plasma <- c("MZB1","JCHAIN","TNFRSF17","ITM2C","DERL3","POU2AF1","IGHA1","TXNDC11","CD79A")
markers.naive = c("MS4A1","CD27","TNFRSF13B","TCL1A","CD38","PPP1R14A","PLD4","SOX4","VPREB3","CD72","IL4R","PLPP5","SELL","FCER2")
markers.memory = c("MS4A1","TCL1A","CD27","TNFRSF13B","IGHD","IGHM","LINC01857","COBLL1","CD1C","AP3B1","CRIP1","CRIP2","TAGLN2","ITGB1","S100A4","S100A10","AHNAK")
#pdf("Mph_marker.pdf",width = 4.08,height = 8.16)
DotPlot(sce, group.by = 'RNA_snn_res.0.01',
        features = unique(genes_to_check))+
  coord_flip()+theme_bw()+#去除背景，旋转图片
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))+#文字90度呈现  
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#5FB4C4",'#A89490','#B74638'))#颜色渐变设置  
#dev.off()

############################################################################
###                               注释                                  ### 
###########################################################################
Idents(sce) <- sce@meta.data$celltype
table(sce$RNA_snn_res.0.2)

sce <- RenameIdents(sce, 
                    "0" = "APOE_Mph",
                    "1" = "TIMP1_Mph",
                    "2" = "FCER1A_Mph",
                    "3" = "SPP1_Mph",
                    "4" =  "SEPP1_Mph"
                    
)

sce$celltype <- Idents(sce)
Idents(sce) <- sce$celltype

qsave(sce,"Mph.qs")


Idents(sce) <- sce@meta.data$celltype
library(COSG)
library(dplyr)
# 记录开始时间
starttime <- Sys.time()

# 运行COSG marker基因识别
COSG_markers <- cosg(
  sce,
  groups = 'all',
  assay = 'RNA',
  slot = 'data',
  mu = 1,
  n_genes_user = 100
)

# 记录结束时间
endtime <- Sys.time()

# 计算耗时
timecost <- endtime - starttime
print(timecost)

COSG_markers$names
markers =  as.data.frame(COSG_markers$names)


write.csv(markers,file = "F:\\Sprouts\\CRC.wsl\\9.Mph.KM\\markers.csv")

############################################################################
###                               细胞分布偏好                           ### 
###########################################################################
#devtools::install_github('satijalab/seurat-data')
#install.packages("epitools")
library(epitools)
library(ComplexHeatmap)
library(reshape2)
library(Seurat)
options(stringsAsFactors = FALSE)
#计算不同分组的细胞数量（样本/细胞群等）
test1=as.matrix(as.data.frame.matrix(table(sce$sample, sce$celltype)))
observe.data <- t(test1) 
expected.data <- expected(observe.data) # 生成期望频数表
plot.data <- observe.data/expected.data # 计算实际细胞数目与期望细胞数目的比值
#col.order <- c("PB", "AL", "HCC", "CHC", "ICC") # 绘图数据的列顺序
#plot.data <- plot.data[, col.order] # 对绘图数据列进行排序
plot.data <- plot.data[order(apply(plot.data, 1, which.min), decreasing = T), ] # 对绘图数据行进行排序
# 如果想深色色块在主对角线使用which.min，在副对角线则使用which.max

# 热图颜色设定,阈值根据需求进行更改
cols <- setNames(object = c("#F9FEAF","#EBDA7B","#E2B255","#D6883A","#B6441F","#A51414"),
                 nm = c("1", "1.5", "2", ">3"))

cols <- setNames(object = c("#faf8cb","#f3bd8f","#d78c4d","#cc5a17"),
                 nm = c("1", "1.5", "2", ">3"))


cell_fun = function(j, i, x, y, width, height, fill) {
  if(plot.data[i,j] < as.numeric(names(cols)[1])) {
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[1], fill = cols[1]))
    grid.text("\u00B1", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(cols)[2])) {
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[2], fill = cols[2]))
    grid.text("+", x = x, y = y)
  } else if(plot.data[i,j] < as.numeric(names(cols)[3])) {
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[3], fill = cols[3]))
    grid.text("++", x = x, y = y)
  } else {
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[4], fill = cols[4]))
    grid.text("+++", x = x, y = y)
  }
}

# 生成热图
hm = Heatmap(plot.data,  col = c("#faf8cb","#f3bd8f","#d78c4d","#cc5a17"),
             cluster_rows = F, cluster_columns = F,
             width = unit(5, "cm"), height = unit(3, "cm"), 
             show_heatmap_legend = F)
# 生成图例
lgd = Legend(labels = names(cols), title = expression(R[o/e]), legend_gp = gpar(fill = cols))
pdf(file="Roe2.pdf",width = 3.39,height = 5.0)
draw(hm, annotation_legend_list = lgd)
dev.off()



############################################################################
###                             换一个风格                         ### 
###########################################################################
library(data.table)
observe.data <- as.matrix(as.data.frame.matrix(table(sce$sample, sce$celltype)))

# 计算Ro/e
expected.data <- expected(observe.data)
Ratio <- observe.data/expected.data
# 生成绘图数据
plot.data <- melt(Ratio)
colnames(plot.data) <- c("sample", "celltype", "Ratio")
plot.data$celltype <- factor(as.character(plot.data$celltype),
                             levels = colnames(Ratio)[order(Ratio[1, ])])


table(sce$sample)
mycolors = c("#9ab882","#584f8a","#b02c3e","#d9792b","#C8E3C1","#A4C9DE","#F6CBDE","#A05528","#9bbaca","#b6acc2","#387091")
color <- setNames(object = mycolors, 
                  nm = c("mBrain", "mLN","nLN","nLung",'PE','tL/B','tLung'))
# 画折线图
ggplot(plot.data, aes(x = celltype, y = Ratio, group = sample, color = sample))+
  geom_point(size = 3) +
  geom_line(linetype = "dashed") +
  scale_color_manual(values = color) +
  theme_classic() + 
  labs(x = "", y = "Ro/e", fill = "position") + 
  theme(axis.text.x = element_text(hjust = 1, size = 10, color = "black",angle = 45),
        axis.title.y = element_text(hjust = 0.5, size = 10, color = "black"),
        axis.text.y = element_text(hjust = 0.5, size = 10, color = "black"),
        axis.ticks = element_line(size = 0.2, color="black"),
        axis.ticks.length = unit(0.2, "cm"),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.9, 0.9), # 自己调整，把图例放到空白处
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10))



############################################################################
###                              PROGENy富集分析                       ### 
###########################################################################
library(Seurat)
library(ggplot2)
library(tidyr)#%>% 使用的%>% 来自这个包
library(dplyr)
library(readr)
library(pheatmap)
library(tibble)
library(progeny)
#library(Magrittr)
Idents(sce) <- sce@meta.data$celltype
sce <- progeny(sce, scale=FALSE, organism="Human", top=500, perm=1, 
               return_assay = TRUE)

load('G:\\1.Sprouts\\LUAD\\9.Mye.sub\\ALL.Celltype.RData')


sce <- Seurat::ScaleData(sce, assay = "progeny") 



progeny_scores_df <- 
  as.data.frame(t(GetAssayData(sce, slot = "scale.data", 
                               assay = "progeny"))) %>%
  rownames_to_column("Cell") %>%
  gather(Pathway, Activity, -Cell) 

CellsClusters <- data.frame(Cell = names(Idents(sce)), 
                            CellType = as.character(Idents(sce)),
                            stringsAsFactors = FALSE)
progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)


summarized_progeny_scores <- progeny_scores_df %>% 
  group_by(Pathway, CellType) %>%
  summarise(avg = mean(Activity), std = sd(Activity))


summarized_progeny_scores_df <- summarized_progeny_scores %>%
  dplyr::select(-std) %>%   
  spread(Pathway, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 


paletteLength = 100
myColor = colorRampPalette(c("#38608b", "white","#743232"))(paletteLength)

progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                      length.out=ceiling(paletteLength/2) + 1),
                  seq(max(summarized_progeny_scores_df)/paletteLength, 
                      max(summarized_progeny_scores_df), 
                      length.out=floor(paletteLength/2)))

progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                        fontsize_row = 10, 
                        color=myColor, breaks = progenyBreaks, 
                        main = "PROGENy (500)", angle_col = 45,
                        treeheight_col = 0,  border_color = NA)


hm = Heatmap(summarized_progeny_scores_df[,-1], c("#38608b", "white","#743232"),
            cluster_rows = F, cluster_columns = F,
            width = unit(10, "cm"), height = unit(6, "cm"), 
            show_heatmap_legend = F)

p+hm


############################################################################
###                              PROGENy富集分析                       ### 
###########################################################################
#1. 加载必要包
library(MASS) 
library(fields)
#2. 提取UMAP坐标和分组信息
umap_df <- as.data.frame(Embeddings(sce, "umap")) 
umap_df$Group <- sce$sample
#3. 定义计算富集分数的函数
get_enrich_score <- function(umap_df, group1, group2) {sub1 <- umap_df[umap_df$Group == group1, 1:2, drop=FALSE] 
sub2 <- umap_df[umap_df$Group == group2, 1:2, drop=FALSE] 
dens1 <- kde2d(sub1$UMAP_1, sub1$UMAP_2, n = 100) 
dens2 <- kde2d(sub2$UMAP_1, sub2$UMAP_2, n = 100) 
enrich <- log2(dens2$z + 1e-6) - log2(dens1$z + 1e-6) 
enrich_score <- fields::interp.surface( obj = list(x = dens1$x, y = dens1$y, z = enrich), loc = umap_df[, c("UMAP_1", "UMAP_2")]) 
return(enrich_score)}

table(sce$sample)
#4. 计算并写入meta.data（举例：比较NC与Prei、Prei与OSCC）
umap_df$tLung_vs_nLung <- get_enrich_score(umap_df, "tLung", "nLung") 
umap_df$mLN_vs_nLN <- get_enrich_score(umap_df, "mLN", "nLN") 
umap_df$tLung_vs_mLN <- get_enrich_score(umap_df, "tLung", "mLN") 
umap_df$mLN_vs_nLung <- get_enrich_score(umap_df, "mLN", "nLung") 

sce$tLung_vs_nLung <- umap_df$tLung_vs_nLung 
sce$tLung_vs_mLN <- umap_df$tLung_vs_mLN
sce$mLN_vs_nLung <- umap_df$mLN_vs_nLung
sce$mLN_vs_nLN <- umap_df$mLN_vs_nLN


#5. 可视化富集分数
FeatureDimPlot( srt = sce, 
                features = "tLung_vs_nLung", 
                pt.size = 0.5, reduction = "umap", 
                theme_use = "theme_blank" )



#5. 可视化富集分数
FeatureDimPlot( srt = sce, 
                features = "tLung_vs_mLN", 
                pt.size = 0.5, reduction = "umap", 
                theme_use = "theme_blank" )


#5. 可视化富集分数
FeatureDimPlot( srt = sce, 
                features = "mLN_vs_nLung", 
                pt.size = 0.5, reduction = "umap", 
                theme_use = "theme_blank" )

#5. 可视化富集分数
FeatureDimPlot( srt = sce, 
                features = "mLN_vs_nLN", 
                pt.size = 0.5, reduction = "umap", 
                theme_use = "theme_blank" )


############################################################################
###                               细胞分布偏好                           ### 
###########################################################################
###########################################################################
###                               细胞分布偏好                           ### 
###########################################################################
###定义mono----
scRNA_mye =  sce
FeatureDimPlot(scRNA_mye,features = 'VCAN')

Monocyte <- c("VCAN", "FCN1", "S100A8", "S100A9", "CD14", "FCGR3A", "CX3CR1")
M2_signature <- c("GPNMB", "APOE", "SPP1", "SELENOP", "TREM2", "MRC1", "CCL18", "MSR1", "VEGFA")
M1_signature <- c("CD86", "TNF", "CXCL9", "CXCL10", "CXCL11")

signature_list_cd4 <- list(
  MonocyteNEW = Monocyte,
  M2_signaturenew = M2_signature,
  M1_signaturenew = M1_signature
)

# 2. 打分
DefaultAssay(scRNA_mye) <- 'RNA'
table(scRNA_mye$celltype)
library(Seurat)
scRNA_mye <- AddModuleScore(
  object = scRNA_mye, 
  features = signature_list_cd4, 
  name = c("Monocytenew", "M2new", "M1new")
)

# 3. 提取感兴趣的分组标签和分数
celltype_col <- "celltype"   # 替换成你的分组列
score_cols <- c("Monocytenew1", "M2new2", "M1new3") # AddModuleScore自动生成

df_cd4 <- scRNA_mye@meta.data[, c(celltype_col, score_cols)]
colnames(df_cd4) <- c("celltype", "Monocyte", "M2", "M1")
df_cd4 <- as.data.frame(df_cd4)

# 4. 计算每个celltype的均值
# 计算每个celltype的均值
unique_types <- unique(df_cd4$celltype)
print(unique_types)
library(dplyr)
df_cd4$celltype <- as.character(df_cd4$celltype)
df_mean <- dplyr::summarise(
  dplyr::group_by(df_cd4, celltype),
  Monocyte = mean(Monocyte, na.rm = TRUE),
  M2 = mean(M2, na.rm = TRUE),
  M1 = mean(M1, na.rm = TRUE)
)
print(df_mean)
dim(df_mean)
print(df_mean)
data <- df_mean
# 5. 三角坐标轴投影

# 定义坐标轴角度（0°, 120°, 240°）
angles <- c(0, 90, 210) * pi / 180

# 计算笛卡尔坐标
data$x <- with(data, Monocyte * cos(angles[1]) + M2 * cos(angles[2]) + M1 * cos(angles[3]))
data$y <- with(data, Monocyte * sin(angles[1]) + M2 * sin(angles[2]) + M1 * sin(angles[3]))

# 创建坐标轴数据
axis_df <- data.frame(
  x_start = rep(0, 3),
  y_start = rep(0, 3),
  x_end = c(1, cos(angles[2]), cos(angles[3])),
  y_end = c(0, sin(angles[2]), sin(angles[3]))
)

# 创建坐标轴标签数据
axis_labels <- data.frame(
  x = c(1.1, 1.1 * cos(angles[2]), 1.1 * cos(angles[3])),
  y = c(0, 1.1 * sin(angles[2]), 1.1 * sin(angles[3])),
  label = c("Naive", "Memory", "Treg")
)

# 假定 data 中含 Cytotoxic、Exhausted、Regulatory、T_celltype 等列

# 1. 对这三列做scale标准化
data_scaled <- data
cols <- c("Monocyte", "M2", "M1")
data_scaled[cols] <- scale(data[cols])

# 2. 定义坐标轴角度
angles <- c(0, 90, 210) * pi / 180

# 3. 根据标准化后的数据计算x、y坐标
data_scaled$x <- with(data_scaled, Monocyte * cos(angles[1]) +
                        M2 * cos(angles[2]) +
                        M1 * cos(angles[3]))
data_scaled$y <- with(data_scaled, Monocyte * sin(angles[1]) +
                        M2 * sin(angles[2]) +
                        M1 * sin(angles[3]))

# 4. 重新生成坐标轴和标签 (和原样即可)
axis_df <- data.frame(
  x_start = rep(0, 3),
  y_start = rep(0, 3),
  x_end   = c(1, cos(angles[2]), cos(angles[3])),
  y_end   = c(0, sin(angles[2]), sin(angles[3]))
)
axis_labels <- data.frame(
  x     = c(1.1, 1.1 * cos(angles[2]), 1.1 * cos(angles[3])),
  y     = c(0,   1.1 * sin(angles[2]), 1.1 * sin(angles[3])),
  label = c("Monocyte", "M2", "M1")
)

# 5. 绘图
library(ggrepel)
library(ggplot2)
ggplot() +
  geom_segment(
    data = axis_df,
    aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
    arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
    linewidth = 0.8,
    color = "gray30"
  ) +
  geom_text(
    data = axis_labels,
    aes(x, y, label = label),
    size = 5,
    color = "darkred"
  ) +
  geom_point(
    data = data_scaled,
    aes(x, y, color = celltype),
    size = 3
  ) +
  geom_text_repel(
    data = data_scaled,
    aes(x, y, label = celltype),
    size = 4,
    box.padding = 0.5,
    max.overlaps = Inf
  ) +
  coord_fixed(xlim = c(-2, 2), ylim = c(-2, 2)) + # 范围可根据数据分布微调
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(1, 1, 1, 1, "cm")
  ) +
  annotate(
    "path",
    x = cos(seq(0, 2*pi, length.out = 100)),
    y = sin(seq(0, 2*pi, length.out = 100)),
    color = "gray80",
    linewidth = 0.3
  ) +
  scale_color_brewer(palette = "Set3")



MarkerNameVector <- c("Monocytenew1", "M2new2", "M1new3")

GroupHeatmap(
  srt = scRNA_mye,
  features = MarkerNameVector, group.by = c("celltype"),
  heatmap_palette = "RdYlBu",add_dot = T, add_bg = F, nlabel = 0, show_row_names = TRUE,show_column_names = T,
  cell_annotation_palette = c("Dark2"),row_names_side = "left"
  #add_dot = TRUE, add_reticle = TRUE
)


####
library(slingshot)
sce <- RunSlingshot(srt = sce, group.by = "celltype", reduction = "UMAP",start = c('TIMP1_Mph'))
FeatureDimPlot(sce, features = paste0("Lineage", 1:3), reduction = "UMAP", theme_use = "theme_blank")

CellDimPlot(sce, group.by = "celltype", reduction = "UMAP", lineages = paste0("Lineage", 1:3), lineages_span = 0.1)

############################################################################
##老俊俊富集分析
###########################################################################
Idents(sce) <- sce@meta.data$celltype
table(Idents(sce))


library(ClusterGVis)
library(scRNAtoolVis)
library(org.Hs.eg.db)
library(ggplot2)
library(Seurat)
library(dplyr)
#Idents(sce) <- sce$celltype
# Find markers for each cluster vs. all others (only positive markers)
sce.markers.all <- FindAllMarkers(sce,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25)

sce.markers.all <- sce.markers.all
# Get top 4 marker genes per cluster
sce.markers <- sce.markers.all %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)


# Prepare data for ClusterGVis visualization
st.data <- prepareDataFromscRNA(object = sce,
                                diffData = sce.markers,
                                showAverage = TRUE)




# 富集分析：对每个 cluster 的基因做 GO BP 富集
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)

# 查看前几条结果
head(enrich, 3)

############################################################################
#绘制条形图:
# 调色板（每个 cluster 一种配色）
palette = c("Grays","Light Grays","Blues2","Blues3","Purples2","Purples3","Reds2","Reds3","Greens2","Greens3")

# 绘图循环：每个 cluster 一个条形图
gglist <- lapply(seq_along(unique(enrich$group)), function(x) {
  tmp <- enrich %>%
    dplyr::filter(group == unique(enrich$group)[x]) %>%
    dplyr::arrange(desc(pvalue))
  
  tmp$Description <- factor(tmp$Description, levels = tmp$Description)
  
  p <- ggplot(tmp) +
    geom_col(aes(x = -log10(pvalue), y = Description, fill = -log10(pvalue)),
             width = 0.75) +
    geom_line(aes(x = log10(ratio), y = as.numeric(Description)), color = "grey50") +
    geom_point(aes(x = log10(ratio), y = Description), size = 3, color = "orange") +
    theme_bw() +
    scale_y_discrete(position = "right",
                     labels = function(x) stringr::str_wrap(x, width = 40)) +
    scale_x_continuous(sec.axis = sec_axis(~., name = "log10(ratio)")) +
    colorspace::scale_fill_binned_sequential(palette = palette[x]) +
    ylab("")
  
  return(p)
})

# 给每个图命名
names(gglist) <- paste("C", 1:6, sep = "")  #注意修改参数

#插入到主图:
# 输出包含嵌入 GO 富集 barplot 的 PDF 图
pdf("sc_ggplot_go.pdf", height = 12, width = 14.5, onefile = FALSE)
visCluster(object = st.data,
                plot.type = "both",
                line.side = "left",
                column_names_rot = 45,
                markGenes = sce.markers$gene,
                cluster.order = 1:6,#注意修改参数
                ggplot.panel.arg = c(5, 0.5, 16, "grey90", NA),
                gglist = gglist)

dev.off()




