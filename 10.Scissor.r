
setwd("G:\\1.Sprouts\\CRC.Pathology\\10.Scissor\\")
#devtools::install_github("jinworks/scAB")
#devtools::install_github('sunduanchen/Scissor')
library(Scissor)
library(Seurat)
library(preprocessCore)
library(qs)
#library(scAB)


source("source_scissor.R")


sce  <- qread("G:\\1.Sprouts\\CRC.Pathology\\3.Caf.sub\\Mph.qs")

sc_dataset = sce
set.seed(123)

cells.use <- sample(
  colnames(sce),
  size = round(ncol(sce) / 2)
)

sc_dataset <- sce[, cells.use]

# 重新计算邻居图
sc_dataset <- FindNeighbors(
  sc_dataset,
  assay = "RNA",
  reduction = "pca"
)


matrix <- read.csv("tumor_matrix.csv",row.names = 1,check.names = FALSE)

clinic <- read.csv("Clinic.csv",row.names = 1)
# # ========== 区分TCGA肿瘤和正常样本 ==========
# # 提取样本类型代码（第14-15位）
# sample_type_code <- substr(colnames(matrix), 14, 15)
# 
# # 区分肿瘤和正常样本
# is_tumor <- sample_type_code >= "01" & sample_type_code <= "09"
# is_normal <- sample_type_code >= "10" & sample_type_code <= "19"
# 
# # 查看样本分布
# cat("样本分布:\n")
# cat("  肿瘤样本数:", sum(is_tumor), "\n")
# cat("  正常样本数:", sum(is_normal), "\n")
# 
# # 提取肿瘤样本矩阵
# tumor_data <- matrix[, is_tumor]
# 
# # 简化列名为 TCGA-AA-3531 格式（前12位）
# colnames(tumor_data) <- substr(colnames(tumor_data), 1, 12)
# 
# 
# 
# matrix = tumor_data
clinic <- clinic[rownames(clinic) %in% colnames(matrix),]

# matrix <- matrix[, colnames(matrix) %in% rownames(clinic)]
# write.csv(matrix,file = "tumor_matrix.csv")
#####################################################################################################
###################################################对矩阵去重#######################################
#####################################################################################################
bulk_dataset = matrix


phenotype = clinic
colnames(phenotype)=c("time","status")
phenotype$time <- as.numeric(phenotype$time)  
phenotype$status <- as.numeric(phenotype$status)  
bulk_dataset <- bulk_dataset[, colSums(is.na(bulk_dataset)) == 0]  

bulk_dataset[] <- lapply(bulk_dataset, function(x) as.numeric(as.character(x)))  
if(T){
  start_time1 <- Sys.time()
  
  # devtools::install_github('sunduanchen/Scissor')
  library(Seurat)
  library(Scissor)
  infos1 <- Scissor(as.matrix(bulk_dataset), 
                    sc_dataset, phenotype, alpha = 0.001, 
                    family = "cox", 
                    Save_file = './Scissor_survival.RData')
  
  end_time1 <- Sys.time()
  execution_time1 <- end_time1 - start_time1
  print(execution_time1)
}




names(infos1)
length(infos1$Scissor_pos)
infos1$Scissor_pos[1:4]
length(infos1$Scissor_neg)

Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
Scissor_select[infos1$Scissor_pos] <- "Scissor+"
Scissor_select[infos1$Scissor_neg] <- "Scissor-"
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
tsne_scissor <- DimPlot(sc_dataset, reduction = 'umap', 
                        group.by = 'scissor',
                        cols = c('grey','royalblue','indianred1'), 
                        pt.size = 0.1, order = c("Scissor-","Scissor+"))
tsne_scissor
table(sc_dataset$scissor,sc_dataset$celltype)

gplots::balloonplot(
  table(sc_dataset$celltype,sc_dataset$scissor)
)
patchwork::wrap_plots(plots = list(celltype,tsne_scissor), ncol = 2)

