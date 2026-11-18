
setwd("G:\\1.Sprouts\\CRC.Pathology\\16.Cellchat\\")
library(CellChat)
library(patchwork)
library(Seurat)
library(qs)
options(stringsAsFactors = FALSE)

#af <- qread("G:\\1.生信豆芽菜\\0.脂肪瘤\\Sarcoma2.qs")
sce <- qread("sce.f.qs")

table(sce@meta.data[["celltype2"]])
Idents(sce) = sce$celltype2


cellchat = createCellChat(object = sce, group.by = "celltype2")
levels(cellchat@idents) #展示以下现在的细胞分组
groupSize <- as.numeric(table(cellchat@idents)) #细胞亚群各组数量
#设置配体受体交互数据库 
CellChatDB = CellChatDB.human #如果是老鼠的话使用内置“CellChatDB.mouse”数据
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")#取出相应分类用作分析数据库
cellchat@DB <- CellChatDB.use#将数据库内容载入cellchat对象中
cellchat <- subsetData(cellchat,features = NULL)#取出表达数据
cellchat <- identifyOverExpressedGenes(cellchat)#寻找高表达的基因#
#future::plan("multiprocess", workers = 4)#可以跳过
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI.human)
#计算通信概率并推断cellchat网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat<- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)#将细胞通讯预测结果以数据框的形式取出
DT::datatable(df.net)
write.csv(df.net,'01.df.net.csv')
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(2,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
#细胞通信网络的可视化
cellchat@netP$pathways
pathways.show <- c("PLAU") 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
# Circle plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Chord diagram,将细胞分为不同的类型，弦图会出现裂开的弧边
#group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
#names(group.cellType) <- levels(cellchat@idents)
#netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

#计算每个配体受体对整体信号通路的贡献，并可视化由单个配体受体对调节的细胞通信
netAnalysis_contribution(cellchat, signaling = pathways.show)
#我们还可以可视化由单个配体受体对调节的细胞-细胞通信。
#我们提供一个函数extractEnrichedLR来提取给定信号通路的所有重要相互作用（L-R对）和相关信号基因。
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = c(2,5)
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")


#自动保存所有推断网络的模块以进行快速探索
pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
for (i in 1:length(pathways.show.all)) {
  # Visualize communication network associated with both signaling pathway and individual L-R pairs
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 4, height = 2, units = 'in', dpi = 300)
}
#我们还可以使用netVisual_bubble显示从某些细胞组到其他细胞组的所有重要相互作用（L-R对）。


table(sce@meta.data[["celltype2"]])

pdf(file = "netVisual_bubble_mash.pdf",width=4.65,height = 6.30)
netVisual_bubble(cellchat, sources.use = 5, 
                 targets.use = c(7,8,9,10,11,12), remove.isolate = FALSE)
dev.off()
#sources.use = 4表示第四个4细胞类型，targets.use = c(6:10)表示从剩下的6到10个细胞的通讯关系，重点在于关注哪个细胞
netVisual_bubble(cellchat, sources.use = 5, 
                 targets.use = c(1,2,3,4,6,7,8,9,10,11,12,14,15,16,17,18,19,20), 
                 signaling = c("PLAU"), remove.isolate = FALSE)
#指定配受体对在细胞中查看
#类似于气泡图，cellchat提供了绘制和弦图的功能，可视化效果不是很好，在配受体数量较多的时候
netVisual_chord_gene(cellchat, sources.use = 5, 
                     targets.use = c(1,2,3,4,6,7,8,9,10,11,12,14,15,16,17,18,19,20), lab.cex = 0.5,legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = 5, 
                     targets.use = c(1,2,3,4,6,7,8,9,10,11,12,14,15,16,17,18,19,20), 
                     legend.pos.x = 15,small.gap = 4)
#sources.use = c(1,2,3,4)表示被选取的细胞，targets.use表示被针对的细胞
netVisual_chord_gene(cellchat,sources.use = 5, 
                     targets.use = c(1,2,3,4,6,7,8,9,10,11,12,14,15,16,17,18,19,20), 
                     #signaling = c("PTN"),
                     legend.pos.x = 8)
#同理，我们可以指定配受体，看看细胞之间的联系，如果基因名称重叠，您可以通过降低small.gap值来调整参数。
#small.gap数字越大，基因名称标签的离散度就越大


#细胞通信网络系统分析
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
pathways.show <- c("PLAU") 
# pathways.show <- c("GALECTIN") 
# pathways.show <- c("PARs") 
# pathways.show <- c("SPP1") 
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 14, height = 3.5, font.size = 10)
#在 2D 空间中可视化占主导地位的发送器（源）和接收器（目标）
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("PLAU"))
gg1 + gg2
#识别对某些细胞组的传出或传入信号贡献最大的信号，可以回答以下问题：哪些信号对某些细胞组的传出或传入信号贡献最大。
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",font.size = 6.5,
                                         width = 10,height = 15,font.size.title = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 6.5,
                                         width = 10,height = 15,font.size.title = 10)
ht1+ht2

save.image('cellchat.RData')
