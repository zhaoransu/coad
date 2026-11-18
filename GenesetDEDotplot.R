
setwd("G:\\1.Sprouts\\CRC.Pathology\\11.GenesetDEDotplot\\")
library(qs)

library(Seurat)
library(magrittr)
library(plyr)
library(openxlsx)
library(GSVA)
library(AUCell)
library(clusterProfiler)
library(limma)
library(ggplot2)
library(dplyr)
library(cowplot)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor



#BiocManager::install("AUCell")
# 单细胞数据的准备
# 目标：生成一个注释好的Seurat对象，以备作图需求
# Seurat对象应包含以下信息
# 样本分组（group）和细胞类型（celltype）

# 设置存放输入数据的路径

Mphsub  <- qread("G:\\1.Sprouts\\CRC.Pathology\\3.Mph.sub\\Mph.qs")
# 读取通路信息，计算各细胞类型的通路得分
# 设置存放输入数据的路径
data.path <- file.path("InputData") 
geneset <- read.gmt("h.all.v7.5.1.symbols.gmt")
geneset <- split(geneset$gene, geneset$term)
genesetInfo <- read.delim("GenesetInfo.txt", sep = ",")
genesetInfo <- subset(genesetInfo, Classification != "")
levels = c("Immune", "Metabolism", "Signaling", "Proliferation")
genesetInfo$Classification <- factor(genesetInfo$Classification, levels)
genesetInfo <- arrange(genesetInfo, genesetInfo$Classification, genesetInfo$geneset)
genesetInfo$geneset <- factor(genesetInfo$geneset, levels = genesetInfo$geneset)
score <- AUCell_run(Mphsub@assays$RNA@counts, geneset)
score <- score@assays@data$AUC



#Mphsub@meta.data$sumgroup = Mphsub@meta.data[["sample_type"]]

# 在不同细胞类型中，对两组进行差异分析
compare="Tumor-Normal" # AC组与PA组进行比较
adjust.method = "bonferroni" # p值矫正方法，设定为bonferroni，用于作图的显著性指标为FDR
DP.list <- lapply(levels(Mphsub$celltype), function(celltype){
  group = setNames(object = Mphsub$sumgroup[Mphsub$celltype == celltype],
                   nm = colnames(Mphsub)[Mphsub$celltype == celltype])
  design = model.matrix(~ 0 + factor(group))
  colnames(design) = levels(factor(group))
  rownames(design) = names(group)
  
  contrast.matrix = makeContrasts(compare, levels = design) # should be Test-Control
  fit = lmFit(score[, names(group)], design)
  fit2 = contrasts.fit(fit, contrast.matrix)
  fit2 = eBayes(fit2) 
  DPs = topTable(fit2, coef=1, n=Inf, adjust.method = adjust.method)
  DPs$celltype = celltype
  DPs$Pathway = rownames(DPs)
  return(DPs)
})
DPs <- do.call(rbind, DP.list)

# 输出差异分析表格
write.table(DPs, file = "output_hallmark.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)









############################################################################################
###################################作图######################################################
# 准备画图数据
DPs <- read.table("output_hallmark.txt", header = T,sep = '\t')
plot.data <- DPs
plot.data <- subset(plot.data, Pathway %in% genesetInfo$geneset)
plot.data$celltype <- factor(plot.data$celltype)
plot.data$Pathway <- factor(plot.data$Pathway, levels = genesetInfo$geneset) # 通路顺序与genesetInfo一致
plot.data$FDR<- cut(plot.data$adj.P.Val, breaks = c(0, 1e-125, 1e-75, 1e-25, 1),
                    include.lowest = T)
plot.data$FDR <- factor(as.character(plot.data$FDR),
                        levels = rev(levels(plot.data$FDR)))
levels(plot.data$Pathway) <- tolower(gsub("HALLMARK_", "", levels(plot.data$Pathway)))

# 准备绘图
color = c("#4682B4", "#FFFFFF", "#CD2626") # logFC对应的颜色
class.color = c("Immune" = "#D58986", "Metabolism" = "#80554C",
                "Signaling" = "#71AC7A", "Proliferation" = "#E8D4B4") # 通路分类对应的颜色
p1 <- ggplot(plot.data, aes(x = Pathway, y = celltype, color = logFC, size = FDR)) +
  geom_point() +
  scale_color_gradient2(low = color[1], mid = color[2], high = color[3]) +
  geom_hline(yintercept = seq(min(as.numeric(plot.data$celltype))-0.5,
                              max(as.numeric(plot.data$celltype))+0.5),
             color = "grey80") +
  geom_vline(xintercept = seq(min(as.numeric(plot.data$Pathway))-0.5,
                              max(as.numeric(plot.data$Pathway))+0.5),
             color = "grey80") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.line = element_blank())

p2 <- ggplot(genesetInfo, aes(x = geneset, y = 1, fill = Classification)) +
  geom_tile() +
  theme_classic() +
  scale_fill_manual(values = class.color) +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        axis.ticks = element_blank(), axis.line = element_blank())

pdf("GenesetDEDotplot.pdf",width = 14,height = 4.25)
plot_grid(p1, p2, ncol = 1, align = 'v', rel_heights = c(10, 1)) # 拼接图像
dev.off()



write.table(DPs, file = "output_hallmark.txt",
            sep = "\t", row.names = F, col.names = T, quote = F)

write.csv(plot.data,file = "plot.data.csv")
plot.data2 <- read.csv("plot.data2.csv")
