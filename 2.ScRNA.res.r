


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


setwd("G:\\1.Sprouts\\CRC.Pathology\\2.ScRNA.res\\")


sce  <- qread("sce.qs")

sce@assays[["RNA"]]@data <- sce@assays[["RNA"]]@counts

Idents(sce) <- sce@meta.data$celltype
Bcs <- subset(sce,idents=c("Bcs"))
Mast <-subset(sce,idents=c("Mast"))
Mye <- subset(sce,idents=c("Mye"))
Tcs <- subset(sce,idents=c("Tcs"))

qsave(Bcs,file = "Bcs.qs")
qsave(Mast,file = "Mast.qs")
qsave(Mye,file = "Mye.qs")
qsave(Tcs,file = "Tcs.qs")

sce2 =sce
sce = sce[,!(sce@meta.data$RNA_snn_res.3 %in% c(44,46,37))] 
############################################################################
###             样本读取
###########################################################################
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")# 计算线粒体占比
sce[["percent.Ribo"]] <- PercentageFeatureSet(sce, pattern = "^RP[SL]")# 计算rRNA占比

sce<-subset(sce,subset =  
              nFeature_RNA < 5000 & 
              percent.mt < 15)

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
            npcs = 50,pc.genes=VariableFeatures)     #PCA分析
##整合方法1，单个样本间进行整合
##max.iter.cluster默认是10，决定运行RunHarmony结果再迭代多少次后完成了收敛
#lamba参数，默认是1，决定整合力度，值越小，整合力度越大
ElbowPlot(sce , ndims=50, reduction="pca") 
pc.num=1:15

############################################################################
###                               1.3数据   harmony                      ### 
###########################################################################
sce   <- sce  %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(sce, 'harmony')

############################################################################
###                               1.4数据   选定分辨率                  ### 
###########################################################################
options(future.globals.maxSize = 256*1000 * 1024^2)

sce=FindNeighbors(sce, dims = 1:15, reduction = "harmony")
for (res in c(0.01, 0.02,0.03,0.05, 0.1, 0.2,0.3, 0.5,1.5,2.0,2.5,3.0
              )) {
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
  RunUMAP(reduction = "harmony", dims = pc.num,min.dist = 0.8) %>%
  #RunTSNE(reduction = "harmony", dims = pc.num) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution = 0.01) %>%
  identity()

# sce <- sce  %>%
#   RunTSNE(reduction = "harmony", dims = pc.num,min.dist = 1.2) %>%
#   #RunTSNE(reduction = "harmony", dims = pc.num) %>%
#   FindNeighbors(reduction = "harmony", dims = pc.num) %>%
#   FindClusters(resolution = 0.1) %>%
#   identity()

#min_dist -- 是低维空间中点之间的最小距离。这个参数控制了UMAP将点紧密地聚集在一起的程度，数值越低，嵌入的点就越紧密。
#更大的min_dist值会使UMAP将点更松散地包装在一起，而专注于保存广泛的拓扑结构。
pc.num=1:15
sce <- sce  %>%
  RunUMAP(reduction = "harmony", n.components = 3,
          dims = pc.num,min.dist = 0.8) %>%
  #RunTSNE(reduction = "harmony", dims = pc.num) %>%
  FindNeighbors(reduction = "harmony", dims = pc.num) %>%
  FindClusters(resolution = 0.01) %>%
  identity()

# sce <- sce  %>%
#   RunTSNE(reduction = "harmony", dims = pc.num,min.dist = 1.2) %>%
#   #RunTSNE(reduction = "harmony", dims = pc.num) %>%
#   FindNeighbors(reduction = "harmony", dims = pc.num) %>%
#   FindClusters(resolution = 0.1) %>%
#   identity()


CellDimPlot3D(srt = sce, group.by = "LUADtype",pt.size = 0.5)

CellDimPlot3D(srt = sce, group.by = "celltype",pt.size = 0.5)

CellDimPlot3D(srt = sce, group.by = "RNA_snn_res.0.2",pt.size = 0.5)

############################################################################
###                               1.5Topgene                           ### 
##########################################################################
Idents(sce) <- sce@meta.data$celltype2
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
  top_n(n = 20, wt = avg_log2FC)


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



############################################################################
###                               COSG计算特征                          ### 
###########################################################################

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
  n_genes_user = 3000
)

# 记录结束时间
endtime <- Sys.time()

# 计算耗时
timecost <- endtime - starttime
print(timecost)

COSG_markers$names
markers =  as.data.frame(COSG_markers$names)


write.csv(markers,file = "markers.csv")
#######################################写出GMT文件
# 只保留前30个marker
markers_top30 <- markers[1:30, ]

# 写出GMT文件
# GMT格式：每行为 基因集名称 \t 描述 \t 基因1 \t 基因2 \t ...

gmt_file <- "markers_top30.gmt"

# 打开文件连接
con <- file(gmt_file, "w")

# 遍历每一列（每个细胞类型）
for (col_name in colnames(markers_top30)) {
  # 获取该列的基因列表（去除NA值）
  genes <- markers_top30[[col_name]]
  genes <- genes[!is.na(genes)]
  
  # 构建GMT格式的行
  # 格式：基因集名称\t描述\t基因1\t基因2\t...
  gmt_line <- paste(c(col_name, col_name, genes), collapse = "\t")
  
  # 写入文件
  writeLines(gmt_line, con)
}

# 关闭文件连接
close(con)

print(paste("GMT文件已保存至:", gmt_file))

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


DimPlot(sce, reduction = "umap",label = F,
        raster=FALSE, group.by = "celltype",pt.size = 0.01,
        cols = mycolors)

DimPlot(sce, reduction = "umap",label = T,
        raster=FALSE, group.by = "RNA_snn_res.0.2",pt.size = 0.1,
        cols = mycolors)+ scale_colour_manual(values = mycolors)


DimPlot(sce, reduction = "umap",label = F,
        raster=FALSE, group.by = "degree",pt.size = 0.1,
        cols = mycolors)+ scale_colour_manual(values = mycolors)



qsave(sce,"sce.qs")
sce <- qread("sce.qs")
############################################################################
###                                                       ### 
###########################################################################
FeaturePlot(object = sce, 
            features = "PLP1",
            cols = c("#dcdcdc", "#ea9d6d", "#de825e"),  # 三色渐变
            pt.size = 0.05,                         # 点的大小
            order = TRUE)                          # 高表达的细胞显示在上层
FeatureDimPlot3D(srt = sce.obj, features = c("CD3D","KRT19","MZB1","CD79A","CD68","DCN","CPA3","VWF"),pt.size = 0.8,force = TRUE)
############################################################################
###                                                    ### 
###########################################################################
library(data.table)

# 文件路径
file1 <- "G:/1.Sprouts/CRC.Pathology/0.RawData/SC.data/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt"
file2 <- "G:/1.Sprouts/CRC.Pathology/0.RawData/SC.data/GSE144735_processed_KUL3_CRC_10X_annotation.txt"

# 读取（推荐 fread，快且稳）
anno1 <- fread(file1)
anno2 <- fread(file2)

rownames(anno1)  = anno1$Index
anno1$prefix <- "CRC2"
rownames(anno1) <- paste0("CRC2_", rownames(anno1))
anno1$cellid = rownames(anno1)

rownames(anno2)  = anno2$Index
anno2$prefix <- "CRC1"
rownames(anno2) <- paste0("CRC1_", rownames(anno2))
anno2$cellid = rownames(anno2)
ann = rbind(anno1,anno2)
# 查看
ann = as.data.frame(ann)
rownames(ann)  = ann$cellid

# 是否完全一致（顺序也一致）
identical(rownames(ann), colnames(sce))
common_cells <- intersect(colnames(sce), rownames(ann))

ann <- ann[common_cells,]
sce <- AddMetaData(sce, metadata = ann)
############################################################################
###                         1.6 看一下marker检验                         ### 
###########################################################################



genes_to_check = c(
  "CCL5","NKG7", "CD3D", "CD3E", "CD8A","CD4",#0：T cell:
  "EPCAM","KRT19","CLDN4","GKN1","PGC",#1:Epi:
  "IGHG1","CD38",'MZB1', 'SDC1',"JCHAIN","SLAMF7",#2:plasma:
  "C1QC","CSF1R","CD68","LYZ","S100A8",#3:Macro:
  'LUM' , 'FGF7',"DCN","COL1A1","CXCL14",#4:Fib:
  "CD79A",'MS4A1','CD19',"BANK1","BLK",#5:B/plasma:
  "CDH5","VWF","FABP4","PLVAP",'PECAM1',#6:Endo:
  "RGS5","TAGLN","ACTA2","MYL9","IGFBP7",#7:Pericyte:
  'CPA3' , 'KIT', 'TPSAB1','TPSB2'#8:Mast cell
)
#pdf("Mph_marker.pdf",width = 4.08,height = 8.16)
DotPlot(sce, group.by = 'RNA_snn_res.0.2',
        features = unique(genes_to_check))+
  coord_flip()+theme_bw()+#去除背景，旋转图片
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=45,hjust = 1,vjust=0.5))+#文字90度呈现  
  scale_color_gradientn(values = seq(0,1,0.2),
                        colors = c("#5FB4C4",'#A89490','#B74638'))#颜色渐变设置  
#dev.off()
sce.obj =sce
sce.obj= sce.obj[, sample(1:ncol(sce.obj),round(ncol(sce.obj)/2)) ]
FeatureDimPlot3D(srt = sce.obj, features = c("CD3D","KRT19","MZB1","CD79A","CD68","DCN","CPA3","VWF"),pt.size = 0.8,force = TRUE)


############################################################################
###                               注释                                  ### 
###########################################################################
Idents(sce) <- sce@meta.data$RNA_snn_res.0.2
table(sce$RNA_snn_res.0.2)

sce <- RenameIdents(sce, 
                    "0" = "Tcs",
                    "1" = "Epi",
                    "2" = "Mye",
                    "3" = "plasma",
                    "4" = "Fib",
                    "5" = "Bcs",
                    "6" = "Epi",
                    "7" = "EC",
                    "8" = "Pericyte",
                    "9" = "Epi",
                    "10"= "Mast"
)

sce$celltype <- Idents(sce)
Idents(sce) <- sce$celltype

qsave(sce,"sce.qs")

############################################################################
###                               细胞分布偏好                           ### 
###########################################################################
###ROE
# >>RO/E----
ROIE <- function(crosstab, filter = NULL){
  ## Calculate the Ro/e value from the given crosstab
  ##
  ## Args:
  #' @crosstab: the contingency table of given distribution
  ##
  ## Return:
  ## The Ro/e matrix 
  if(is.null(filter)){filter = 10}
  rowsum.matrix <- matrix(0, nrow = nrow(crosstab), ncol = ncol(crosstab))
  rowsum.matrix[,1] <- rowSums(crosstab)
  rowsum.matrix[rowsum.matrix <= filter] <- 0
  colsum.matrix <- matrix(0, nrow = ncol(crosstab), ncol = ncol(crosstab))
  colsum.matrix[1,] <- colSums(crosstab)
  colsum.matrix[colsum.matrix <= filter] <- 0
  allsum <- sum(crosstab)
  roie <- divMatrix(crosstab, rowsum.matrix %*% colsum.matrix / allsum)
  row.names(roie) <- row.names(crosstab)
  colnames(roie) <- colnames(crosstab)
  roie <- roie[rowSums(roie)!=0,]
  return(roie)
}

divMatrix <- function(m1, m2){
  ## Divide each element in turn in two same dimension matrixes
  ##
  ## Args:
  #' @m1: the first matrix
  #' @m2: the second matrix
  ##
  ## Returns:
  ## a matrix with the same dimension, row names and column names as m1. 
  ## result[i,j] = m1[i,j] / m2[i,j]
  dim_m1 <- dim(m1)
  dim_m2 <- dim(m2)
  if( sum(dim_m1 == dim_m2) == 2 ){
    div.result <- matrix( rep(0,dim_m1[1] * dim_m1[2]) , nrow = dim_m1[1] )
    row.names(div.result) <- row.names(m1)
    colnames(div.result) <- colnames(m1)
    for(i in 1:dim_m1[1]){
      for(j in 1:dim_m1[2]){
        if(m2[i,j] == 0){
          div.result[i,j] <- 0
        }else{
          div.result[i,j] <- m1[i,j] / m2[i,j]
        }
      }
    }   
    return(div.result)
  }
  else{
    warning("The dimensions of m1 and m2 are different")
  }
}

# >>ROIE_plot----
ROIE_plot <- function(ROIE_table, max = 3, color = 1, type = 1) {
  if (type == 1) {
    # Melt the ROIE table
    ROIE_table <- melt(ROIE_table)
    ROIE_table$text <- round(ROIE_table$value, 2)
    ROIE_table$value[ROIE_table$value > max] <- max
    
    # Define color palettes based on the color parameter
    if (color == 1) {
      colors <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(100)
    } else if (color == 2) {
      colors <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(100)
    } else if (color == 3) {
      colors <- colorRampPalette(c("#3783c4", "#8daddb", "#fefdfd", "#ed8171", "#c7252c"))(100)
    } else {
      stop("Invalid color option. Choose 1, 2, or 3.")
    }
    
    # Create the plot
    p <- ggplot(data = ROIE_table, aes(Var2, Var1, fill = value)) +
      geom_tile() + 
      scale_fill_gradientn(name = "Ro/e", colours = colors, na.value = "lightgrey") +
      geom_text(aes(label = text)) +
      labs(x = "", y = "") +
      theme_cowplot(font_size = 12) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), 
            axis.line = element_blank(), 
            axis.ticks = element_blank()) +
      guides(fill = guide_colourbar(barwidth = 1, barheight = 10))
    
    return(p)
  } else if (type == 2) {
    # Melt the ROIE table
    ROIE_table <- melt(ROIE_table)
    ROIE_table$value[ROIE_table$value > 9] <- 9
    
    # Convert to matrix
    plot.data <- acast(ROIE_table, Var1 ~ Var2, value.var = "value")
    
    # Define colors for cell_fun
    cols <- setNames(object = c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
                     nm = c("1", "5", "7",">7"))
    
    # Custom cell function
    cell_fun <- function(j, i, x, y, width, height, fill) {
      value <- plot.data[i, j]
      if (value < as.numeric(names(cols)[1])) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[1], fill = cols[1]))
        grid.text("\u00B1", x = x, y = y)
      } else if (value < as.numeric(names(cols)[2])) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[2], fill = cols[2]))
        grid.text("+", x = x, y = y)
      } else if (value < as.numeric(names(cols)[3])) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[3], fill = cols[3]))
        grid.text("++", x = x, y = y)
      } else {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[4], fill = cols[4]))
        grid.text("+++", x = x, y = y)
      }
    } 
    # Generate heatmap
    hm <- Heatmap(plot.data, cell_fun = cell_fun,
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  width = unit(6, "cm"), height = unit(3, "cm"), 
                  show_heatmap_legend = FALSE)
    
    # Generate legend
    lgd <- Legend(labels = names(cols), title = expression(R[o/e]), legend_gp = gpar(fill = cols))
    
    # Draw heatmap and legend
    draw(hm, annotation_legend_list = list(lgd))
    
  } else if (type == 3) {
    # Melt the ROIE table
    ROIE_table <- melt(ROIE_table)
    ROIE_table$value[ROIE_table$value > max] <- max
    
    # Convert to matrix
    plot.data <- acast(ROIE_table, Var1 ~ Var2, value.var = "value")
    
    # Define colors for cell_fun
    cols <- setNames(object = c("#FEE6CE", "#FDC08C", "#F5904B", "#E6550D"),
                     nm = c("1", "2", "3", ">3"))
    
    # Custom cell function
    cell_fun <- function(j, i, x, y, width, height, fill) {
      value <- plot.data[i, j]
      if (value < as.numeric(names(cols)[1])) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[1], fill = cols[1]))
        grid.text("\u00B1", x = x, y = y)
      } else if (value < as.numeric(names(cols)[2])) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[2], fill = cols[2]))
        grid.text("+", x = x, y = y)
      } else if (value < as.numeric(names(cols)[3])) {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[3], fill = cols[3]))
        grid.text("++", x = x, y = y)
      } else {
        grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = cols[4], fill = cols[4]))
        grid.text("+++", x = x, y = y)
      }
    } 
    
    # Generate heatmap
    hm <- Heatmap(plot.data, cell_fun = cell_fun,
                  cluster_rows = FALSE, cluster_columns = FALSE,
                  width = unit(10, "cm"), height = unit(3, "cm"), 
                  show_heatmap_legend = FALSE)
    
    # Generate legend
    lgd <- Legend(labels = names(cols), title = expression(R[o/e]), legend_gp = gpar(fill = cols))
    
    # Draw heatmap and legend
    draw(hm, annotation_legend_list = list(lgd))
  } else {
    stop("Invalid type option. Choose 1 or 2.")
  }
}

library(reshape2)
library(cowplot )
library(RColorBrewer)
plot.data<- ROIE(t(table(sce$celltype, sce$Class)))
#plot.data<- ROIE((table(scRNAdata$seurat_clusters, scRNAdata$Group2)))
p <- ROIE_plot(plot.data, max = 3, color = 1, type = 1)
print(p)


############################################################################
###                             换一个风格  绘制细胞数量               ### 
###########################################################################
###注意重启R，应为有的功能互相串扰导致报错
###注意重启R，应为有的功能互相串扰导致报错
###注意重启R，应为有的功能互相串扰导致报错
sce$sumgroup = sce$Class
metadata <- FetchData(sce, c("sumgroup", "celltype"))
cell_counts <- metadata %>% dplyr::count(sumgroup, celltype)
total_cells_by_celltype <- cell_counts %>% group_by(celltype) %>% summarise(total = sum(n))
total_cells <- cell_counts %>% group_by(sumgroup) %>% summarise(total = sum(n))
cell_counts <- cell_counts %>% left_join(total_cells, by = "sumgroup")
cell_counts <- cell_counts %>% mutate(percentage = n / total * 100)
write.csv(cell_counts, 'cell_percentage.csv')



# 绘制细胞数量柱状图
p1 <- ggplot(total_cells_by_celltype, aes(x = total, y = celltype, fill = celltype)) + 
  geom_bar(stat = "identity") + 
  theme_classic() + 
  theme(
    text = element_text(size = 15, colour = "black"), 
    axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"), 
    axis.text.y = element_text(colour = "black"), 
    legend.position = "none"
  ) + 
  labs(x = '# cells', y = NULL) + 
  scale_fill_manual(values = mycolors) +  # 使用 scale_fill_manual 指定颜色
  geom_text(aes(label = total),           # 添加数字标签
            hjust = 0.3,                  # 调整标签位置（右对齐）
            size = 4,                     
            color = "black")

p1



# 绘制分组比例柱状图
# 修改细胞比例计算方式（每个细胞类型内各组的比例）
cell_counts <- metadata %>% 
  dplyr::count(sumgroup, celltype) %>%                  # 按组和细胞类型计数
  group_by(celltype) %>%                      # 按细胞类型分组
  mutate(
    total_celltype = sum(n),                  # 计算细胞类型总细胞数
    percent_in_celltype = n / total_celltype * 100  # 计算组内占比
  ) %>% 
  ungroup()

p2 <- ggplot(cell_counts, aes(x = celltype, y = n, fill = sumgroup)) + 
  geom_bar(
    stat = "identity", 
    position = "fill",                         # 自动将 y 轴标准化为 0-1
    width = 0.8, 
    size = 0, 
    colour = "black"
  ) + 
  geom_text(
    aes(label = sprintf("%.1f", percent_in_celltype)), 
    position = position_fill(vjust = 0.5), 
    size = 4, 
    color = "black"
  ) + 
  scale_fill_manual(values = mycolors) + 
  scale_y_continuous(  
    labels = scales::percent  # 使用 scales 包的百分比格式化函数
  ) + 
  theme_classic() + 
  labs(x = NULL, y = "# percentage") + 
  coord_flip() + 
  theme(
    axis.ticks.length = unit(0.1, "cm"), 
    text = element_text(size = 15, colour = "black"), 
    axis.text.x = element_text(
      angle = 45, 
      hjust = 1, 
      colour = "black"
    ), 
    axis.text.y = element_text(colour = "black")
  )

p2


# 将两个图合并
library(patchwork)

p1 + p2 + plot_layout(nrow = 1)

ggsave('barplot_counts_percentage.pdf', width = 11.5, height = 5.5)



qsave(sce,"sce.f.qs")
