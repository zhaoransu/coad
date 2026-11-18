
setwd("G:\\1.Sprouts\\CRC.Pathology\\13.Epi.KM\\")
############################################################################
###                                                                   ### 
###########################################################################
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
library(GSEABase)
library(randomcoloR)
library(AnnoProbe)
if (!require('R.utils')) install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto' )
library(BiocManager)
#BiocManager::install("impute")



rt = read.csv("tumor_matrix.csv",row.names = 1,check.names = F)
os = read.csv("Clinic.csv",row.names = 1)
os <- os[rownames(os) %in% colnames(rt),]
colnames(os) = c("futime","fustat")

library(GSVA)
library(GSVA)
library(limma)
library(GSEABase)
library(limma)
#可视化
exp_data1 = rt
exp_data1 = as.matrix(exp_data1)

geneSet=getGmt("D1.afgmt.gmt", geneIdType=SymbolIdentifier())#读取我们所提供的免疫细胞基因文件
##免疫浸润分析
gsva <- gsva(exp_data1, geneSet, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}#定义ssGSEA_Score矫正函数
norm_ssGSEA=normalize(gsva)#对ssGSEA_Score进行矫正
#norm_ssGSEA=rbind(id=colnames(norm_ssGSEA),norm_ssGSEA)
write.csv(gsva,file="KM.norm_ssGSEA_score.csv")




df = read.csv("KM.norm_ssGSEA_score.csv",row.names = 1)
df = t(df)

kmdf = df

 
############################################################################
### 准备数据
###########################################################################
library(survminer)
library(survival)
kmdf = as.data.frame(kmdf)


col_names <- colnames(kmdf)
# 遍历所有列
for (target_col in col_names) {
  # 确保只处理数值型列
  if (is.numeric(kmdf[[target_col]])) {
    # 计算上下四分位数
    q1 <- quantile(kmdf[[target_col]], 0.25, na.rm = TRUE)
    q3 <- quantile(kmdf[[target_col]], 0.75, na.rm = TRUE)
    
    # 构建新的风险列名
    risk_col_name <- paste0(target_col, ".risk")
    
    # 生成新列，根据四分位数分配 highrisk、lowrisk 或 MID
    kmdf[[risk_col_name]] <- ifelse(kmdf[[target_col]] >= q3, "highrisk",
                                    ifelse(kmdf[[target_col]] <= q1, "lowrisk", "MID"))
  }
}

# 删除所有包含 "MID" 标签的行
# 找出所有以 ".risk" 结尾的列
risk_cols <- grep("\\.risk$", colnames(kmdf), value = TRUE)

# 对于每个风险列，筛选出不含 "MID" 的行
# for (risk_col in risk_cols) {
#   kmdf <- kmdf[kmdf[[risk_col]] != "MID", ]
# }


kmdf = cbind(os,kmdf)

write.csv(kmdf,file="kmdf.csv")


kmdf = read.csv("kmdf.csv",row.names = 1,check.names = F)

# ###########################################################################

# ###########################################################################

# ==================== 完整修复脚本：生成文件 + KM分析 + P值&HR汇总 ====================
library(survival)
library(survminer)
library(dplyr)

# ==================== 步骤1：读取原始数据并动态获取所有 .risk 列 ====================
cat("正在读取原始数据并动态识别所有 .risk 列...\n")

original_df <- read.csv("kmdf.csv", row.names = 1, check.names = FALSE)

# 动态获取所有以 ".risk" 结尾的列名（更稳健，不会因数字范围变化而出错）
risk_cols <- colnames(original_df)[grepl("\\.risk$", colnames(original_df))]

cat("检测到以下亚群列：", paste(risk_cols, collapse = ", "), "\n\n")

# ==================== 步骤2：生成并清洗每个 .risk 文件 ====================
cat("正在生成并清洗所有亚群文件...\n")

for (col in risk_cols) {
  # 提取当前亚群所需的三列
  sub_df <- original_df[, c("futime", "fustat", col)]
  
  # 只保留 highrisk 和 lowrisk，删除包含 MID 的行
  sub_df <- sub_df[sub_df[[col]] %in% c("highrisk", "lowrisk"), ]
  sub_df <- na.omit(sub_df)
  
  # 生成文件名（去掉 ".risk" 后缀）
  subgroup_id <- gsub("\\.risk$", "", col)
  file_name <- paste0(subgroup_id, ".risk.csv")
  
  # 保存为独立 CSV 文件
  write.csv(sub_df, file_name, row.names = TRUE, quote = FALSE)
  
  cat("已生成：", file_name, "（样本数：", nrow(sub_df), "）\n")
}

cat("\n✅ 所有 .risk.csv 文件生成并清洗完成！\n\n")

# ==================== 步骤3：批量 KM 分析 + 收集 P值和HR ====================
cat("正在进行 KM 生存分析并计算 HR...\n")

results <- data.frame(
  Subgroup = character(),
  pValue = character(),
  HR = numeric(),
  Lower95 = numeric(),
  Upper95 = numeric(),
  stringsAsFactors = FALSE
)

for (col in risk_cols) {
  # 构造对应文件名
  subgroup_id <- gsub("\\.risk$", "", col)
  file_name <- paste0(subgroup_id, ".risk.csv")
  
  # 读取清洗后的文件
  kmdf <- read.csv(file_name, row.names = 1, check.names = FALSE)
  
  # 数据类型转换
  kmdf$futime <- as.numeric(kmdf$futime)
  kmdf$fustat <- as.numeric(kmdf$fustat)
  
  # 生存差异检验（p值）
  diff <- survdiff(Surv(futime, fustat) ~ kmdf[[col]], data = kmdf)
  pValue <- 1 - pchisq(diff$chisq, df = 1)
  pValue_fmt <- format(signif(pValue, 4), scientific = TRUE)
  
  # 计算 HR 和 95% CI（使用 coxph）
  cox_fit <- coxph(Surv(futime, fustat) ~ kmdf[[col]], data = kmdf)
  hr_summary <- summary(cox_fit)
  
  HR <- hr_summary$conf.int[1, "exp(coef)"]
  Lower95 <- hr_summary$conf.int[1, "lower .95"]
  Upper95 <- hr_summary$conf.int[1, "upper .95"]
  
  # 保存结果
  results <- rbind(results, data.frame(
    Subgroup = col,
    pValue = pValue_fmt,
    HR = round(HR, 3),
    Lower95 = round(Lower95, 3),
    Upper95 = round(Upper95, 3),
    stringsAsFactors = FALSE
  ))
  
  # 控制台输出单个亚群结果
  cat(paste0("亚群 ", col, "  p值 = ", pValue_fmt, 
             "   HR = ", round(HR, 3), 
             " (95%CI: ", round(Lower95, 3), " - ", round(Upper95, 3), ")\n"))
  
  # # 生成并保存 KM 生存曲线
  # fit <- survfit(Surv(futime, fustat) ~ kmdf[[col]], data = kmdf)
  # 
  # pdf(paste0("KM_plot_", subgroup_id, ".risk.pdf"), width = 8, height = 6)
  # print(ggsurvplot(fit, 
  #                  data = kmdf,
  #                  pval = TRUE,
  #                  risk.table = TRUE,
  #                  xlab = "Time (days)",
  #                  ylab = "Survival Probability",
  #                  legend.title = paste0("Risk Group (", col, ")"),
  #                  palette = c("#E41A1C", "#377EB8")))
  # dev.off()
}

# ==================== 步骤4：生成 P值和HR汇总表格 ====================
cat("\n正在生成 P值 + HR 汇总表格...\n")

write.csv(results, "KM_Summary_Pvalue_HR.csv", row.names = FALSE, quote = FALSE)

# 在控制台以表格形式打印
print(results)

cat("\n🎉 所有处理已完成！\n")
cat("📁 生成的文件包括：\n")
cat("   - 所有 X.risk.csv（已清洗，仅保留 highrisk/lowrisk）\n")
cat("   - KM_plot_X.risk.pdf（各亚群生存曲线）\n")
cat("   - KM_Summary_Pvalue_HR.csv（P值和HR汇总表格）\n")

# ###########################################################################

# ###########################################################################
kmdf = read.csv('F:\\Sprouts\\LUAD\\11.Caf.sub\\gpc3.csv')

kmdf$Caf_GPC3.risk

kmdf = read.csv("4.risk.csv",row.names = 1,check.names = F)
kmdf$futime <- as.numeric(kmdf$futime)  
kmdf$fustat <- as.numeric(kmdf$fustat) 

diff=survdiff(Surv(futime, fustat) ~ Caf_GPC3.risk,data = kmdf)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ Caf_GPC3.risk, data = kmdf)

#pdf(file="survival.pdf",onefile = FALSE,width = 8,height =8)
ggsurvplot(fit,               # 创建的拟合对象
           data = kmdf,         # 指定变量数据来源
           conf.int = F,     # 显示置信区间
           pval = TRUE,         # 添加P值
           surv.median.line = "hv",  # 添加中位生存时间线
           risk.table = TRUE,   # 添加风险表
           risk.table.height = 0.25,
           xlab = "Time in Days",
           ylab = "Overall Survival Probability",
           legend = c(0.8, 0.75),  # 指定图例位置
           legend.title = "",   # 设置图例标题，这里设置不显示标题，用空格替代
           legend.labs = c("Low risk","High risk"),  # 指定图例分组标签
           break.x.by = 365,    # 设置x轴刻度间距
           palette = c("blue", "red")
           )     # 自定义调色板
#dev.off()
#summary(fit)


# ###########################################################################
#绘制玫瑰图
# ###########################################################################

results = read.csv("results.all.celltype.km.csv")
results$negLog10PValue <- -log10(results$pValue)

ggplot(results, aes(x = celltype, y = negLog10PValue, fill = celltype)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = mycolors) + 
  geom_hline(yintercept = 1.34, # 根据自己的数据调整
             color = "grey", linetype="dashed") + 
  coord_polar() + 
  theme_classic() + 
  theme(axis.text = element_blank(), axis.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank())



# 首先按照negLog10PValue对celltype进行排序
results$celltype <- factor(results$celltype,
                           levels = results$celltype[order(results$negLog10PValue)])

# 然后使用你的绘图代码
ggplot(results, aes(x = celltype, y = negLog10PValue, fill = celltype)) + 
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = mycolors) + 
  geom_hline(yintercept = 1.34, 
             color = "grey", linetype = "dashed") + 
  coord_polar() + 
  theme_classic() + 
  theme(axis.text = element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(), 
        axis.ticks = element_blank())


# ###########################################################################
#森林图
# ###########################################################################

library(forestploter)
library(dplyr)

library(dplyr) # 确保 dplyr 已经加载

data_fp <- results %>%
  mutate(
    # 1. 强制转换为数值型
    pValue = as.numeric(pValue),
    HR = as.numeric(HR),
    Lower95 = as.numeric(Lower95),
    Upper95 = as.numeric(Upper95),
    
    # 2. 格式化 P 值 
    `P Value` = ifelse(pValue < 0.001, "< 0.001", sprintf("%.3f", pValue)),
    
    # 3. 格式化 HR (95% CI) 字符串
    `HR (95% CI)` = sprintf("%.2f (%.2f - %.2f)", HR, Lower95, Upper95),
    
    # 4. 创建一个空列，用来给画图留出位置
    ` ` = paste(rep(" ", 30), collapse = " ")
  ) %>%
  # 关键修改：在这里加上 dplyr:: 强制指定使用 dplyr 包的 select
  dplyr::select(Subgroup, `P Value`, ` `, `HR (95% CI)`, HR, Lower95, Upper95)

# 检查一下
head(data_fp)

# 检查一下数据处理是否成功（不出错就可以继续跑后面的画图代码了）
head(data_fp)

# 2. 设置森林图的美化主题
tm <- forest_theme(
  base_size = 10,          # 基础字体大小
  base_family = "sans",    # 字体
  ci_pch = 15,             # 形状：方形
  ci_col = "#2166AC",      # 图形主色调 (经典蓝色)
  ci_fill = "#ec6538",       
  ci_alpha = 0.8,
  ci_lty = 1,
  ci_lwd = 2,              # 线条粗细
  refline_lwd = 1,         # 参考线粗细
  refline_lty = "dashed",  # 参考线虚线
  refline_col = "red"      # 参考线颜色红色
)


# 3. 绘制带有数据表的森林图
forest(
  data = data_fp[, c(1, 2, 3, 4)], # 选择左边要展示为文本的列 (包括空列)
  est = data_fp$HR,                # HR
  lower = data_fp$Lower95,         # 下界
  upper = data_fp$Upper95,         # 上界
  ci_column = 3,                   # 告诉函数把图形画在第几列 (我们刚才创建的空列)
  ref_line = 1,                    # 参考线位置
  xlim = c(0.5, 2.5),              # X轴范围 (根据你的数据调整)
  ticks_at = c(0.5, 1, 1.5, 2, 2.5),# X轴刻度
  theme = tm                       # 应用主题
)

# 显示图形
print(p)
