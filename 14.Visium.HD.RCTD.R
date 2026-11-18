# 设置工作目录
setwd("G:\\1.Sprouts\\CRC.Pathology\\14.VisiumHD\\")

# 加载所需 R 包
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
#devtools::install_github("dmcable/spacexr", build_vignettes = FALSE)
library(spacexr)
library(Matrix)
library(stringr)
library(qs)
library(arrow)
# devtools::install_local("D:\\Rpackage_github_install_location\\spacexr.zip")
# 在此处开始你的分析代码 ---------------------------------
# 例如：
# data <- Read10X(data.dir = "path_to_10x_data")
# seurat_obj <- CreateSeuratObject(counts = data)
sce  <- qread("G:\\1.Sprouts\\CRC.Pathology\\5.Celltrek\\sce.f.qs")

# 读取 Seurat 对象
ImmuRef <- sce

# 提取 counts 矩阵
CountRef <- GetAssayData(ImmuRef, slot = "counts")

# 提取细胞类型注释
CTRef <- ImmuRef@meta.data$celltype2

# 转换为因子
CTRef <- as.factor(CTRef)

# 设置细胞名作为名称
names(CTRef) <- rownames(ImmuRef@meta.data)

# 构建 reference 对象
# reference <- Reference(
#   counts = CountRef[, names(CTRef)],
#   labels = CTRef,
#   libsize = colSums(CountRef)
# )

reference <- Reference(
  CountRef[, names(CTRef)],
  CTRef,
  colSums(CountRef)
)

# ===============================
# Deconvolution with RCTD (HD)
# ===============================
counts <- Read10X_h5('G:\\1.Sprouts\\CRC.Pathology\\0.RawData\\HD_data\\P2CRC\\filtered_feature_bc_matrix.h5')


# 2. 读取空间坐标
coords <- read_parquet(
  "G:\\1.Sprouts\\CRC.Pathology\\0.RawData\\HD_data\\P2CRC\\Metadata.parquet",
  as_data_frame = TRUE
)

coords <- as.data.frame(coords)
rownames(coords) <- coords$barcode

# 3. 对齐 barcode 顺序
coords <- coords[colnames(counts), ]

# 4. 仅保留 x / y 坐标
coords <- coords[, 3:4]

# 5. 计算每个 spot 的 UMI 数
nUMI <- colSums(counts)

# 6. 构建 SpatialRNA 对象
puck <- SpatialRNA(
  coords = coords,
  counts = counts,
  nUMI = nUMI
)

barcodes <- colnames(puck@counts)

# 7. 创建并运行 RCTD
myRCTD <- create.RCTD(
  puck,
  reference = reference,
  max_cores = 20
)

myRCTD <- run.RCTD(
  myRCTD,
  doublet_mode = "doublet"
)

# 8. 导出结果
RCTD_result <- data.frame(
  ID = rownames(myRCTD@results$results_df),
  myRCTD@results$results_df
)

sample_dir = "G:\\1.Sprouts\\CRC.Pathology\\14.VisiumHD\\"
write.table(
  RCTD_result,
  file = file.path(sample_dir, "RCTD_result.P1.txt"),
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)





saveRDS(myRCTD,file = 'myRCTD.p1crc.rds')
