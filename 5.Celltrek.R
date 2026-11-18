setwd("G:\\1.Sprouts\\CRC.Pathology\\5.Celltrek\\")
#library(devtools)
#install_github("navinlabcode/CellTrek")
options(stringsAsFactors = F)
#library("CellTrek")
library("dplyr")
library("Seurat")
library("viridis")
library("ConsensusClusterPlus")
#BiocManager::install("ConsensusClusterPlus")
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(qs)
library(Seurat)
library(patchwork)
library(qs)
library(qs)
library(ggsci)
mycolors = c(pal_d3("category20")(20),
             pal_d3("category20b")(20),
             pal_d3("category20c")(10),
             pal_d3("category10")(10))
############################################################################
###                            1.Celltrek开始                           ###
###########################################################################
# ==========================================
# 1. 全局变量设置 (只需修改这里)
# ==========================================
inputpath   <- "G:\\1.Sprouts\\CRC.Pathology\\5.Celltrek\\SP.Visium\\10X.2\\"  # 你的 Spatial 数据所在文件夹路径
inputfile   <- "filtered_feature_bc_matrix.h5"            # h5 表达矩阵文件名
sample_name <- "10X.2"                                   # 样本名称 (用于切片命名、ID替换和最终保存)

# ==========================================
# 2. 读取 10x 空间转录组数据
# ==========================================
# 自动使用上方定义的 inputpath, inputfile 和 sample_name
sce_st_cortex <- Load10X_Spatial(
  data.dir = inputpath,
  filename = inputfile,
  assay = "Spatial",
  slice = sample_name
)

# 查看元数据和初始形态
head(sce_st_cortex@meta.data)
SpatialDimPlot(sce_st_cortex, alpha = 0)

# 和单细胞数据一样，人为更改这些命名，使其与 sample_name 统一
sce_st_cortex@project.name <- sample_name
Idents(sce_st_cortex)      <- sample_name
sce_st_cortex$orig.ident   <- sample_name


plot1 <- VlnPlot(sce_st_cortex, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(sce_st_cortex, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
# ggsave(filename = "../Output/test.pdf")

p3 <- VlnPlot(sce_st_cortex, features = "nFeature_Spatial", pt.size = 0.1, cols = "tomato", group.by = "orig.ident") + NoLegend()
p4 <- FeatureScatter(sce_st_cortex, feature1 = "nCount_Spatial", feature2= "nFeature_Spatial", group.by = "orig.ident") + NoLegend()
wrap_plots(p3, p4)

# 调整 future.globals.maxSize 选项来增加允许的全局变量大小，防止 SCTransform 报错
options(future.globals.maxSize = 1000 * 1024^3)

# 运行 SCTransform
sce_st_cortex <- SCTransform(sce_st_cortex, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

# 提取表达量变化最高的10个基因
top10 <- head(VariableFeatures(sce_st_cortex), 10)
print("Top 10 Variable Features:")
print(top10)

# 自动根据 sample_name 生成保存的文件名 (例如: "P24.T2.qs")
# 如果你想保存在之前定义的 inputpath 路径下，可以将 file 改为: file = file.path(inputpath, paste0(sample_name, ".qs"))
save_filename <- paste0(sample_name, ".qs")

qsave(sce_st_cortex, file = save_filename)

############################################################################
###                            空间数据处理
###########################################################################

sce_st_cortex <- qread("colon3.qs")
SpatialDimPlot(sce_st_cortex)

############################################################################
###                           单细胞数据处理
###########################################################################
sce_sc <- qread("sce.f.qs")
table(sce_sc$celltype2)

sce_sc = sce_sc[,(sce_sc@meta.data$celltype2 %in% c("SEPP1_Mph","SPP1_Mph","FCER1A_Mph",
                                                    "TIMP1_Mph","APOE_Mph","TIMP13_Mph",
                                                    "KRT16_Epi","PLCG2_Epi","ITLN1_Epi",
                                                    "FABP1_Epi","LDHB_Epi","SLPI_Epi"
                                                    ))] 

Idents(sce_sc) = sce_sc$celltype2 
table(Idents(sce_sc))
########备份
sce_sc = sce_sc2
########
sce_sc= sce_sc[, sample(1:ncol(sce_sc),round(ncol(sce_sc)/40)) ]


#sce_sc <- sce_sc[, !((Idents(sce_sc)) %in% c("Fib_ACTG1","Fib_CCL11","Fib_CCL5","Fib_CFD"))]
## Rename the cells/spots with syntactically valid names
sce_st_cortex <- RenameCells(sce_st_cortex, new.names=make.names(Cells(sce_st_cortex)))
sce_sc <- RenameCells(sce_sc, new.names=make.names(Cells(sce_sc)))



sce_traint <- traint(st_data=sce_st_cortex, sc_data=sce_sc, sc_assay='RNA', cell_names='celltype2')

## We can check the co-embedding result to see if there is overlap between these two data modalities
#DimPlot(sce_traint, group.by = "type") 

sce_celltrek <- celltrek(st_sc_int=sce_traint, int_assay='traint', sc_data=sce_sc, sc_assay = 'RNA', 
                         reduction='pca', intp=T, 
                         intp_pnt=10000, 
                         intp_lin=F, 
                         nPCs=30, 
                         ntree=100, 
                         dist_thresh=0.4, 
                         top_spot=5, 
                         spot_n=5, 
                         repel_r=20, 
                         repel_iter=20, 
                         keep_model=T)$celltrek



sce_celltrek$celltype2 <- factor(sce_celltrek$celltype2, levels=sort(unique(sce_celltrek$celltype2)))

# CellTrek::celltrek_vis(sce_celltrek@meta.data %>% dplyr::select(coord_x, coord_y, celltype:id_new),
#                       sce_celltrek@images$anterior1@image, sce_celltrek@images$anterior1@scale.factors$lowres)




sce_sgraph_KL <- CellTrek::scoloc(sce_celltrek, col_cell='celltype2', use_method='KL', eps=1e-50)

sce_sgraph_KL_mst_cons <- sce_sgraph_KL$mst_cons


pheatmap::pheatmap(sce_sgraph_KL_mst_cons,color=colorRampPalette(c("#F0F2EE","#2E9298","#156486"))(100),border=FALSE)

# sce_sgraph_KL_mst_cons = read.csv("sce_sgraph_KL_mst_cons.csv",row.names = 1)
# write.csv(sce_sgraph_KL_mst_cons,file = "sce_sgraph_KL_mst_cons.csv")

#mycolors = c("#9ab882","#b02c3e")
#pdf("celltrek_PM6.pdf",width = 6.71,height = 6.36)
mycolors = c("#9ab882","#584f8a","#b02c3e","#d9792b","#A05528","#ce75ac",
             "#C8E3C1","#A4C9DE","#F6CBDE","#A05528","#9bbaca","#b6acc2","#387091")
SpatialDimPlot(sce_celltrek, label.size = 10,stroke = 0.1,
               group.by = "celltype2",
               image.alpha = 0.7,
               pt.size.factor = 1.2)+ 
  scale_fill_manual(values = mycolors)





# 定义 12-gene TLS signature
tls_12_genes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", 
                  "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")

# 在 Seurat 对象中打分
sce_celltrek <- AddModuleScore(sce_celltrek, 
                                features = list(tls_12_genes), 
                                name = "TLS_Score")










SpatialPlot(sce_celltrek, features = 'TLS_Score1')






SpatialPlot(sce_celltrek, features = 'TGFB1')