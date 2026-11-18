# 定义函数
perform_nichenet_analysis <- function(B4_mat, lr_network, ligand_target_matrix, geneset_oi, top_n = 30, output_pdf = "SCLCepiEEpiC2_NicheNetr.pdf") {
  # 计算表达的基因
  expressed_genes_receiver <- B4_mat %>% apply(1, function(x) { sum(x > 0) / length(x) }) %>% .[. >= 0.1] %>% names()
  background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
  
  # 定义潜在的配体和受体
  ligands <- lr_network %>% pull(from) %>% unique()
  expressed_ligands <- ligands
  receptors <- lr_network %>% pull(to) %>% unique()
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  lr_network_expressed <- lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors)
  potential_ligands <- lr_network_expressed %>% pull(from) %>% unique()
  
  # 进行NicheNet的配体活性分析
  ligand_activities <- predict_ligand_activities(geneset = geneset_oi, 
                                                 background_expressed_genes = background_expressed_genes, 
                                                 ligand_target_matrix = ligand_target_matrix, 
                                                 potential_ligands = potential_ligands)
  best_upstream_ligands <- ligand_activities %>% top_n(top_n, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
  
  # 推断顶级配体的靶基因并在热图中可视化
  active_ligand_target_links_df <- best_upstream_ligands %>% lapply(get_weighted_ligand_target_links, geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()
  active_ligand_target_links_df <- na.omit(active_ligand_target_links_df)
  
  active_ligand_target_links <- prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.25)
  
  order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
  order_targets <- active_ligand_target_links_df$target %>% unique()
  order_targets <- intersect(order_targets, rownames(active_ligand_target_links))  
  vis_ligand_target <- active_ligand_target_links[order_targets, order_ligands] %>% t()
  color_palette <- colorRampPalette(c("#4b82bc", "#fefdfd","#bb3030", "#951d26"))(100)
  
  # 创建热图并应用新的颜色渐变
  B4_specific <- vis_ligand_target %>%
    make_heatmap_ggplot("APOE_Mph maturation ligands", "PLCG2_Epi-specific maturation signature", color = "#bb3030", legend_position = "top", x_axis_position = "top", legend_title = "Regulatory potential") + 
    scale_fill_gradientn(colors = color_palette, breaks = c(0, 0.005, 0.01)) + 
    theme(axis.text.x = element_text(face = "italic"))
  
  # 保存结果到PDF文件
  pdf(output_pdf, width = 7, height = 8)
  print(B4_specific)
  dev.off()
}


library(Seurat)
library(ggplot2)
library(future)
library(tidyverse)
library(nichenetr)
library(qs)
options(future.globals.maxSize = 10*1024^3)

setwd("G:\\1.Sprouts\\CRC.Pathology\\17.Nichesnet\\")

sce <- qread("G:\\1.Sprouts\\CRC.Pathology\\16.Cellchat\\sce.f.qs")


Idents(sce) <- sce$celltype2
allCells=names(Idents(sce))
allType = levels(Idents(sce))

# choose_Cells <- unlist(lapply(allType, function(x){
#   cgCells <- allCells[Idents(sce) == x]   #注意修改参数
#   cgSize <- length(cgCells)
#   if (cgSize > 500) {
#     cgCells <- sample(cgCells, 500)
#   }
#   else {
#     cgCells <- cgCells
#   }
#   cgCells
# }))
# 
# cg_sce = sce[, allCells %in% choose_Cells] #注意修改参数
# cg_sce
# table(Idents(cg_sce))
# Bcell <- cg_sce

Bcell = sce
Idents(Bcell) <- Bcell$celltype2
B4 <- subset(Bcell,idents='APOE_Mph')

## take our B4
B4_mat <- as.matrix(B4@assays$RNA@data)

### Read in NicheNet’s ligand-target prior model, ligand-receptor network and weighted integrated networks
ref_dir <- "G:\\1.Sprouts\\CRC.Pathology\\17.Nichesnet\\"
ligand_target_matrix <- readRDS(paste0(ref_dir,"ligand_target_matrix.rds"))
lr_network <- readRDS(paste0(ref_dir,"lr_network.rds"))
weighted_networks <- readRDS(paste0(ref_dir,"weighted_networks.rds"))

pbmc.markers.all <- read.csv("G:\\1.Sprouts\\CRC.Pathology\\17.Nichesnet\\sce.ALL_markers.csv")
### Define the gene set of interest and a background of genes

# pbmc.markers.all$gene <- pbmc.markers.all$X

B_marker <- pbmc.markers.all %>%
  dplyr::filter(cluster=="PLCG2_Epi")%>%
  dplyr::top_n(n = 100, wt = avg_log2FC)

B_marker$gene <- rownames(B_marker)
B_marker$gene <- B_marker$X

geneset_oi <- B_marker$gene %>% .[. %in% rownames(ligand_target_matrix)]
geneset_oi <- grep("-", geneset_oi, value = TRUE, invert = TRUE)
### Define a background of genes
perform_nichenet_analysis(B4_mat, lr_network, ligand_target_matrix, geneset_oi, top_n = 30, 
                          output_pdf = "G:\\1.Sprouts\\CRC.Pathology\\17.Nichesnet\\NICHNET.pdf")





FeaturePlot(sce,features = 'DLL4')








.libPaths(c("D:\\R_tools_install\\R-4.4.2\\Seuratobject_V5", .libPaths()))
library(SeuratObject)
packageVersion("SeuratObject")
library(SCP)


c("CITED2", "ELF3", "HES1", "HEXIM1", "HIST1H2BG", "ID2", "IER2", "IER5", "JUNB", "MALAT1", "PLEC", "STAT3", "UBC", "WDR74", "PLCG2", "LDLR", "SQSTM1", "NFKBIZ", "AHI1", "POLR2A")

all_nichenet_genes <- c("RTN4", "TTR", "DLL1", "DLL4", "CCL5", "APOB", "ADAM15", "CNTN1", "DUSP18", "NID1", "LIF", "CLDN11", "TSHB", "CTSG", "PVR", "WNT3A", 
                        "PSEN1", "KLK3", "EFNA5", "FGF8", "CCL26", "NTN3", "CLSTN1", "NOV", "GMFB", "FGF21", "GZMB", "RARRES2", "ICAM5")

scRNA_mye = qread('G:\\1.Sprouts\\CRC.Pathology\\3.Mph.sub\\Mph.qs')



pdf('GroupHeatmap.pdf', width = 6.79, height = 8)
GroupHeatmap(
  srt = scRNA_mye,
  features = all_nichenet_genes, group.by = c("celltype"),
  heatmap_palette = "RdYlBu",add_dot = F, add_bg = F, nlabel = 0, show_row_names = TRUE,show_column_names = T,cluster_rows = T,
  cell_annotation_palette = c("Dark2"),row_names_side = "left"
  #add_dot = TRUE, add_reticle = TRUE
)
dev.off()
