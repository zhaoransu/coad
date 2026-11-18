
setwd("G:\\1.Sprouts\\CRC.Pathology\\4.Sc.ML\\")
library(Seurat)
library(tidyverse)
library(Matrix)
library(stringr)
library(dplyr)
library(patchwork)
library(ggplot2)
library(pROC)
library(ggplot2)
library(keras)
library(dplyr)
library(caret)
library(MASS)
library(class)
library(rpart)
library(randomForest)
library(xgboost)
library(glmnet)
library(e1071)
library(ComplexHeatmap)
library(kernlab)
library(gbm)
library(lightgbm)
library(ada)
library(qs)
source("F1S.R")
set.seed(123)

# 评分前缀
af_name <- "EMT"

#0.8作为训练集
train_ration <- 0.8 # 训练集比例
heatmapcolor <- c("#BD3C29", "#AE7000","#925E9FFF","#00468BFF") # 设置颜色
selectResult = c(    # 指定ROC/SVM输出的结果
  "RF",
  "SVM"
)
final_select <- "RF" # 最终预测使用的模型


.libPaths(c("D:\\R_tools_install\\R-4.4.2\\Seuratobject_V5", .libPaths()))
library(SeuratObject)
packageVersion("SeuratObject")
library(SCP)
CellDimPlot3D(srt = af, group.by = "Predict",pt.size = 1)
CellDimPlot(srt = af, group.by = "celltype2",pt.size = 1)






sce  <- qread("G:\\1.Sprouts\\CRC.Pathology\\2.ScRNA.res\\sce.qs")


Idents(sce) <- sce@meta.data$celltype
Epi <- subset(sce,idents=c("Epi"))


#########################################################
#########################################################
sce = Epi
af = sce
Idents(sce) <- sce$celltype
Idents(af) <- af$Predict
CAFC3 <- subset(af,ident="positive")
cell_names <- colnames(CAFC3)
sce=SetIdent(object = sce,cells = cell_names,value = 'EMT.p')
sce@meta.data$Predict<- sce@active.ident#
table(sce$Predict)

EMT.p <- subset(sce,idents=c("EMT.p"))

qsave(EMT.p,file = 'G:\\1.Sprouts\\CRC.Pathology\\12.Epi.sub\\EMT.p.qs')

qsave(sce,file = 'G:\\1.Sprouts\\CRC.Pathology\\12.Epi.sub\\sce.EMT.p.qs')
#########################################################
#########################################################
af <- qread("G:\\1.Sprouts\\CRC.Pathology\\3.Mph.sub\\Mph.qs")

# 1. 定义基因向量
ecm <- c("ACTA2", "TAGLN", "BGN", "COL8A1", "COL15A1", "IGFBP7", "TPM1", "TPM2", 
         "COL10A1", "POSTN", "MYL9", "COL13A1", "COL14A1", "MYH11", "MYLK", 
         "ACTG2", "FN1", "LUM", "DCN", "VCAN", "COL5A1", "COL5A2", 
         "COL63A", "PDGFA", "MMP2", "COL1A2", "RGS5")


emt <- c(
  "COL3A1","COL5A2","COL5A1","FBN1","COL1A1","FN1","COL6A3","SERPINE1",
  "COL1A2","COL4A1","COL4A2","VCAN","IGFBP3","TGFBI","SPARC","LUM",
  "LAMC1","LOX","LAMC2","CCN2","TAGLN","COL7A1","LOXL2","COL6A2",
  "ITGAV","THBS2","COL16A1","NNMT","TPM1","CDH2","MMP2","COL11A1",
  "THBS1","FAP","BGN","SERPINH1","FSTL1","POSTN","THY1","SPP1",
  "TNC","TFPI2","NID2","ITGB5","MMP3","VIM","LOXL1","FBLN5",
  "COL12A1","ELN","CDH11","COMP","SPOCK1","BMP1","IL32","LAMA3",
  "TIMP1","QSOX1","TIMP3","VCAM1","CCN1","EDIL3","CALD1","MAGEE1",
  "FBLN1","SGCB","ECM1","LAMA2","FSTL3","TPM2","INHBA","DAB2",
  "EMP3","BASP1","ITGA5","MGP","VEGFA","CXCL1","WNT5A","SDC1",
  "PLOD2","PCOLCE","GREM1","ITGB1","COL5A3","RHOB","HTRA1","FGF2",
  "SNTB1","GADD45A","MEST","LRRC15","TNFRSF11B","CD59","ACTA2",
  "EFEMP2","MATN2","PCOLCE2","SERPINE2","GPC1","ABI3BP","FUCA1",
  "SLIT3","LAMA1","PMEPA1","COL8A2","FBN2","IGFBP2","PFN2","SDC4",
  "CD44","GADD45B","CXCL8","GLIPR1","ANPEP","P3H1","VEGFC","MMP14",
  "SGCD","PLOD1","MATN3","MYL9","SLC6A8","CALU","PRRX1","TNFRSF12A",
  "FMOD","ID2","GEM","PLAUR","MYLK","TGFB1","SFRP1","PLOD3","IL6",
  "APLP1","FBLN2","MSX1","PTX3","FZD8","JUN","FERMT2","DKK1",
  "SNAI2","DST","TPM4","DCN","GJA1","PMP22","IGFBP4","COPA",
  "LRP1","ITGA2","FLNA","MFAP5","PTHLH","TGFBR3","SFRP4","LGALS1",
  "RGS4","CDH6","SAT1","NT5E","DPYSL3","PPIB","TGM2","SGCG",
  "ITGB3","PDLIM4","CTHRC1","ECM2","CRLF1","AREG","IL15","MCM7",
  "GAS1","PRSS2","CADM1","OXTR","SCG2","CXCL6","MMP1","TNFAIP3",
  "CAPG","CAP2","MXRA5","FOXC2","NTM","ENO2","FAS","BDNF",
  "ADAM12","PVR","CXCL12","PDGFRB","SLIT2","NOTCH2","COLGALT1",
  "GPX7","WIPF1"
)

# 2. 关键修改：将向量转化为“有名字的列表”。你可以把 "ECM_Signature" 换成你想要的通路名字
gs <- list(ECM_Signature = emt) 

file <- "EMT.gmt"

# 3. 定义写出函数
write.gmt <- function(gs, file){
  sink(file)
  lapply(names(gs), function(i){
    # i 是基因集的名字 (这里是 "ECM_Signature")
    # 'NA' 是GMT文件的第二列(描述列)的占位符
    # gs[[i]] 是具体的基因向量
    cat( paste(c(i, 'NA', gs[[i]]), collapse='\t') )
    cat('\n')
  })
  sink()
}

# 4. 执行写出
write.gmt(gs, file)



signature_list_cd4 <- list(
  ecm = ecm
)

# 2. 打分
DefaultAssay(af) <- 'RNA'
table(af$celltype)
library(Seurat)
af <- AddModuleScore(
  object = af, 
  features = signature_list_cd4, 
  name = c("ecm")
)

af <- qread("G:\\1.Sprouts\\CRC.Pathology\\3.Mph.sub\\Mph.qs")

FeaturePlot(af, features = 'ecm1')
#########################################################
#########################################################
#af <- qread("E:\\2.WorkStudio\\2.RWorkstudio\\Merge.qs")
af <- qread("G:\\1.Sprouts\\CRC.Xenium\\4.ScRNAseq.res\\sce.f.qs")
af = af[, sample(1:ncol(af),round(ncol(af)/40)) ]

gmt <- clusterProfiler::read.gmt("EMT.gmt")
gene_sets <- split(gmt$gene, gmt$term) 

# 保留 Seurat 对象中存在的基因
seurat_genes <- rownames(af)
gene_sets_filtered <- lapply(gene_sets, function(genes) {
  intersect(genes, seurat_genes)
})

# 计算基因集得分
af <- AddModuleScore(
  object = af,
  features = gene_sets_filtered,
  name = af_name
)

score_n <- paste0(af_name, 1:length(gene_sets_filtered))
# 简单看看评分的分布情况
VlnPlot(af, features = score_n, ncol = 4)  # ncol调整列数
FeaturePlot(af, features = score_n, ncol = 4)

##################################################### 01.高斯分布
######################################
#########################
library(mclust)

df <- af@meta.data %>% as.data.frame()
df <- df[,grep(af_name, colnames(df))]
#
cell_classification <- matrix(NA, nrow = nrow(df), ncol = ncol(df))
rownames(cell_classification) <- rownames(df)
colnames(cell_classification) <- colnames(df)

# 对每个基因集评分拟合高斯分布并分类
for (gene_set in colnames(df)) {
  scores <- df[[gene_set]]
  model <- Mclust(scores, G = 2) # 拟合高斯模型
  high_class <- which.max(model$parameters$mean) # 高置信度细胞
  cell_classification[, gene_set] <- ifelse(model$classification == high_class, 1, 0)
}

# 构建样本集合
high_confidence_cells <- rownames(cell_classification)[rowSums(cell_classification) == ncol(cell_classification)]
low_confidence_cells <- rownames(cell_classification)[rowSums(cell_classification) == 0]
if (length(high_confidence_cells)>length(low_confidence_cells)) {
  high_confidence_cells <- sample(high_confidence_cells, length(low_confidence_cells))
}else{
  low_confidence_cells <- sample(low_confidence_cells, length(high_confidence_cells))
}

##################################################### 02.构建特征矩阵
######################################
#########################
raw_af <- af
af$confidence_group <- ifelse(Cells(af) %in% high_confidence_cells, "High", 
                              ifelse(Cells(af) %in% low_confidence_cells, "Low", NA))
af <- subset(af, subset = confidence_group %in% c("High", "Low"))  # 保留标记分组的细胞
# 鉴定前500个差异基因
Idents(af) <- af$confidence_group
markers <- FindMarkers(af, ident.1 = "High", ident.2 = "Low", only.pos = FALSE)
markers <- markers[!grepl("-", rownames(markers)),]
top_markers <- head(rownames(markers[order(markers$avg_log2FC, decreasing = TRUE), ]), 100)
#
feature_matrix <- as.data.frame(t(as.matrix(af@assays$RNA@data[top_markers, ])))
feature_matrix$confidence_label <- ifelse(rownames(feature_matrix) %in% high_confidence_cells, 1, 0)

# 按标签分组细胞
high_cells <- rownames(feature_matrix[feature_matrix$confidence_label == 1, ])
low_cells <- rownames(feature_matrix[feature_matrix$confidence_label == 0, ])

# 分别对高置信和低置信细胞进行抽样
train_high <- sample(high_cells, size = floor(train_ration * length(high_cells)))
train_low <- sample(low_cells, size = floor(train_ration * length(low_cells)))
train_cells <- c(train_high, train_low)
#
remaining_high <- setdiff(high_cells, train_high)
remaining_low <- setdiff(low_cells, train_low)

# 分别抽样验证集 1 和验证集 2，确保平衡
val1_high <- sample(remaining_high, size = floor(0.5 * length(remaining_high)))
val1_low <- sample(remaining_low, size = floor(0.5 * length(remaining_low)))
val1_cells <- c(val1_high, val1_low)
val2_cells <- setdiff(c(remaining_high, remaining_low), val1_cells)
#
train_exp <- feature_matrix[train_cells, top_markers]
train_labels <- feature_matrix[train_cells, "confidence_label"]
#
val1_matrix <- feature_matrix[val1_cells, top_markers]
val2_matrix <- feature_matrix[val2_cells, top_markers]
test_explist <- list(val1_matrix, val2_matrix)
#
all_labels <- data.frame(V1=rownames(feature_matrix), V2=feature_matrix$confidence_label)
rownames(all_labels) <- all_labels$V1


# ##################################################### 字符串纠错
# ######################################
# # 假设 train_exp 是一个数据框
# # 获取所有以 "Rik" 结尾的列名
# rik_columns <- grep("Rik$", colnames(train_exp), value = TRUE)
# 
# # 删除这些列
# train_exp <- train_exp[, !(colnames(train_exp) %in% rik_columns)]
# val1_matrix <- val1_matrix[, !(colnames(val1_matrix) %in% rik_columns)]
# val2_matrix <- val2_matrix[, !(colnames(val2_matrix) %in% rik_columns)]
# 
# # train_exp <- subset(train_exp, select = -c(`5330417C22Rik`))
# # val1_matrix <- subset(val1_matrix, select = -c(`5330417C22Rik`))
# # val2_matrix <- subset(val2_matrix, select = -c(`5330417C22Rik`))
# 
# 
# 
# # 假设 train_exp1 和 train_exp2 是两个矩阵或数据框
# # 检查两个矩阵的列名是否一致
# identical(colnames(val2_matrix), colnames(val1_matrix))
# identical(colnames(train_exp), colnames(val1_matrix))
##################################################### 02.多种机器学习模型
######################################
#########################
af <- afmulti_machine(train_exp, train_labels, test_explist, all_labels, raw_af, top_markers, selectResult, final_select)
# 看一下最终模型的预测结果降维图
DimPlot(af, reduction = "umap", group.by = "Predict")

mycolors = c("#9ab882","#b02c3e","#584f8a","#d9792b","#C8E3C1","#A4C9DE","#F6CBDE","#A05528","#9bbaca","#b6acc2","#387091")
DimPlot(af, reduction = "umap",label = F,
        raster=FALSE, group.by = "Predict",pt.size = 0.05,
        cols = mycolors)+ scale_colour_manual(values = mycolors)


# 预测结果储存在af$Predict中，positive意味着是
table(af$Predict,af$celltype)

library(Scissor)
gplots::balloonplot(
  table(af$celltype,af$Predict)
)

save.image("af.RData")





















# 提示词：
# 这四张图的数值，绘制了四张图，但是我希望绘制一个汇总的大热图来展示，
# 你能提取其中的数值，
# 并整理成可以绘制大热图的数值形式吗？并给出分裂热图的绘制代码？
# 从左到右依次为：all_Accuracy_statical，all_AUC_statical，all_F-score_statical，all_Recall_statical



# ==========================================
# 1. 加载包
# ==========================================
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# ==========================================
# 2. 准备数据 (使用之前的 df 数据)
# ==========================================
df <- data.frame(
  Model = rep(c("SVM", "AdaBoost", "LDA", "XGBoost", "Logistic", "RF", 
                "StepWise-Logistic", "ENR", "Ridge", "Lasso", "NB", "DT"), 4),
  Metric = rep(c("Accuracy", "AUC", "F-score", "Recall"), each = 12),
  Train_set = c(0.971, 0.995, 0.927, 1.000, 0.931, 1.000, 0.917, 0.907, 0.896, 0.912, 0.869, 0.880,
                0.971, 0.995, 0.927, 1.000, 0.931, 1.000, 0.917, 0.907, 0.896, 0.912, 0.869, 0.880,
                0.971, 0.995, 0.927, 1.000, 0.934, 1.000, 0.921, 0.914, 0.904, 0.917, 0.865, 0.883,
                0.975, 0.994, 0.926, 1.000, 0.976, 1.000, 0.971, 0.979, 0.979, 0.979, 0.840, 0.906),
  Test_set_1 = c(0.933, 0.916, 0.899, 0.910, 0.893, 0.899, 0.876, 0.876, 0.876, 0.871, 0.854, 0.831,
                 0.933, 0.916, 0.899, 0.910, 0.893, 0.899, 0.876, 0.876, 0.876, 0.871, 0.854, 0.831,
                 0.933, 0.917, 0.898, 0.913, 0.896, 0.902, 0.880, 0.885, 0.885, 0.878, 0.856, 0.830,
                 0.944, 0.933, 0.888, 0.944, 0.921, 0.933, 0.910, 0.955, 0.955, 0.933, 0.865, 0.820),
  Test_set_2 = c(0.939, 0.917, 0.933, 0.917, 0.922, 0.906, 0.911, 0.906, 0.894, 0.900, 0.850, 0.794,
                 0.939, 0.917, 0.933, 0.917, 0.922, 0.906, 0.911, 0.906, 0.894, 0.900, 0.850, 0.794,
                 0.939, 0.919, 0.933, 0.921, 0.925, 0.910, 0.916, 0.911, 0.902, 0.905, 0.846, 0.800,
                 0.944, 0.944, 0.922, 0.967, 0.956, 0.956, 0.967, 0.967, 0.967, 0.956, 0.822, 0.822)
)

model_order <- c("SVM", "AdaBoost", "LDA", "XGBoost", "Logistic", "RF", 
                 "StepWise-Logistic", "ENR", "Ridge", "Lasso", "NB", "DT")
df$Model <- factor(df$Model, levels = model_order)

get_matrix <- function(metric_name) {
  mat <- df %>% 
    dplyr::filter(Metric == metric_name) %>% 
    dplyr::select(Model, Train_set, Test_set_1, Test_set_2) %>% 
    tibble::column_to_rownames("Model") %>% 
    as.matrix()
  return(mat)
}

mat_acc    <- get_matrix("Accuracy")
mat_auc    <- get_matrix("AUC")
mat_fscore <- get_matrix("F-score")
mat_recall <- get_matrix("Recall")

# ==========================================
# 3. 设置丰富的高级注释 (Annotations)
# ==========================================

# 3.1 设定模型的类别 (用于行注释)
model_class <- c("SVM" = "Other", "AdaBoost" = "Tree/Ensemble", "LDA" = "Linear",
                 "XGBoost" = "Tree/Ensemble", "Logistic" = "Linear", "RF" = "Tree/Ensemble",
                 "StepWise-Logistic" = "Linear", "ENR" = "Linear", "Ridge" = "Linear",
                 "Lasso" = "Linear", "NB" = "Other", "DT" = "Tree/Ensemble")

# 提取排好序的模型类别
class_ordered <- model_class[model_order]

# 3.2 创建左侧的【行注释】(区分算法类型)
row_anno <- rowAnnotation(
  Algorithm = class_ordered,
  col = list(Algorithm = c("Tree/Ensemble" = "#4DAC26", # 绿色
                           "Linear" = "#D01C8B",        # 紫色
                           "Other" = "#FDAE61")),       # 橙色
  annotation_name_gp = gpar(fontface = "bold")
)

# 3.3 创建顶部的【列注释】(区分训练集/测试集)
col_anno <- HeatmapAnnotation(
  Dataset_Type = c("Train", "Test", "Test"),
  col = list(Dataset_Type = c("Train" = "#2C7BB6", "Test" = "#D7191C")),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontface = "bold", fontsize = 10)
)

# 热图主体颜色
col_fun <- colorRamp2(c(0.79, 0.90, 1.00), c("#3caea3", "white", "#f48c86"))

# ==========================================
# 4. 创建带列注释的单个热图函数
# ==========================================
create_ht <- function(mat, title) {
  Heatmap(mat, 
          name = title, 
          col = col_fun,
          cluster_rows = FALSE, 
          cluster_columns = FALSE, 
          column_title = title,
          column_title_gp = gpar(fontsize = 12, fontface = "bold"),
          column_names_rot = 45,
          top_annotation = col_anno, # ⭐在这里加入列注释
          rect_gp = gpar(col = "gray80", lwd = 1),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.3f", mat[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
          },
          show_heatmap_legend = (title == "Accuracy") # 只显示一个数值图例
  )
}

# ==========================================
# 5. 拼接热图并全局输出
# ==========================================
# 注意这里：用 "+" 号把左侧的行注释 (row_anno) 加在最前面！
ht_list <- row_anno + 
  create_ht(mat_acc, "Accuracy") + 
  create_ht(mat_auc, "AUC") + 
  create_ht(mat_fscore, "F-score") + 
  create_ht(mat_recall, "Recall")

# 绘图，并添加全局主标题
draw(ht_list, ht_gap = unit(5, "mm"), 
     row_title = "Machine Learning Models", 
     row_title_gp = gpar(fontsize = 14, fontface = "bold"),
     column_title = "Comprehensive Performance Evaluation of 12 Machine Learning Models", # 全局大标题
     column_title_gp = gpar(fontsize = 16, fontface = "bold", col = "darkblue"))


































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
library(ggplot2)
plot.data <- ROIE(t(table(af$celltype, af$Predict)))
#plot.data<- ROIE((table(scRNAdata$seurat_clusters, scRNAdata$Group2)))
p <- ROIE_plot(plot.data, max = 3, color = 1, type = 1)
print(p)



############################################################################
###                             换一个风格  绘制细胞数量               ### 
###########################################################################
###注意重启R，应为有的功能互相串扰导致报错
###注意重启R，应为有的功能互相串扰导致报错
###注意重启R，应为有的功能互相串扰导致报错
sce$sumgroup = sce$Tissue
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
