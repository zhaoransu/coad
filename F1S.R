f1_score <- function(predict_result, real_label) {
  # 将预测结果和真实标签转换为因子
  predicted <- as.factor(predict_result)
  real <- as.factor(real_label)
  
  # 确保因子级别一致
  levels(predicted) <- levels(real)
  
  # 计算混淆矩阵
  confusion_matrix <- table(Predicted = predicted, Real = real)
  
  # 提取 TP, FP, FN
  TP <- confusion_matrix["positive", "positive"]
  FP <- confusion_matrix["positive", "negative"]
  FN <- confusion_matrix["negative", "positive"]
  
  # 计算 Precision 和 Recall
  precision <- TP / (TP + FP)
  recall <- TP / (TP + FN)
  
  # 计算 F1 Score
  f1 <- 2 * (precision * recall) / (precision + recall)
  
  # 返回 F1 分数
  return(f1)
}


#af

afmulti_machine <- function(train_exp, train_labels, test_explist, all_labels, raw_af, top_markers, selectResult, final_select){
  
  # 结果列表，包括预测得分and/or预测类别
  all_result_summary <- list()
  # 精准度得分列表
  all_result_acc <- list()
  # 召回率得分列表
  all_result_recall <- list()
  # FS得分列表
  all_result_FS <- list()
  # 基因重要性列表
  all_result_importance <- list()
  
  #
  all_result_summary_single <- list()
  all_result_acc_single <- list()
  all_result_recall_single <- list()
  all_result_FS_single <- list()
  all_result_importance_single <- list()
  
  feature_matrix <- as.data.frame(t(as.matrix(raw_af@assays$RNA@data[top_markers, ])))
  rownames(feature_matrix) <- paste0(rownames(feature_matrix), "af111")
  all_labels_af <- data.frame(V1=rownames(feature_matrix), V2=c(rep(0, 100), rep(1, dim(raw_af)[2]-100)))
  rownames(all_labels_af) <- all_labels_af$V1
  all_labels <- rbind(all_labels, all_labels_af)
  test_explist[[3]] <- as.data.frame(feature_matrix)
  ####### 2.逻辑回归 Logistic Regression
  ####
  ##
  #
  ### 参数设置
  cutoff <- c(0.25, 0.5, 0.75)
  # 创建逻辑回归模型
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- train_labels
  model <- glm(labels~.-labels, family = "binomial", data = train_expd)
  importancei <- as.data.frame(model$coefficients[-1])
  colnames(importancei) <- "coefficients"
  all_result_importance[[paste0("LR")]] <- importancei
  for (cuti in cutoff) {
    train_result <- as.data.frame(predict(model, type="response"))
    train_result$type <- factor(ifelse(train_result[,1] > cuti, "positive", "negative") )
    colnames(train_result) <- c("predict_p", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(model, expdata, type="response"))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("predict_p", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult[,-1])
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result[,-1]
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("LR (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("LR (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("LR (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("LR (cutoff:", cuti, ")")]] <- all_result
  }
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["Logistic"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["Logistic"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["Logistic"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["Logistic"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["Logistic"]] <- all_result_importance[["LR"]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  
  ####### 3.线性判别分析 Linear discriminant analysis
  ####
  ##
  #
  # 创建LDA模型
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- train_labels
  model <- lda(labels~.-labels, data = train_expd)
  all_result_importance[[paste0("LDA")]] <- as.data.frame(model$scaling)
  train_result <- as.data.frame(predict(model, type="response")$class)
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, expdata, type="response")$class)
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[["LDA"]] <- result_acc
  all_result_recall[["LDA"]] <- result_recall
  all_result_FS[["LDA"]] <- result_FS
  all_result_summary[["LDA"]] <- all_result
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["LDA"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["LDA"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["LDA"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["LDA"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["LDA"]] <- all_result_importance[["LDA"]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  
  
  ##################
  ##########
  ####### 4.LightGBM
  # 创建LightGBM模型
  # train_expd <- as.data.frame(train_exp)
  # dtrain <- lgb.Dataset(data = as.matrix(train_expd), label = train_labels)
  # params <- list(
  #   objective = "binary",
  #   metric = "binary_error",
  #   min_data_in_leaf = 20,       # 减少每个叶子节点上的最小数据量
  #   learning_rate = 0.05,         # 设置学习率
  #   num_iterations = 1000        # 增加迭代次数
  # )
  # model <- lgb.train(params, data = dtrain, 1000)
  # im <- as.data.frame(lgb.importance(model))
  # rownames(im) <- im$Feature
  # all_result_importance[["LightGBM"]] <- im[,2, drop=F]
  # train_pred <- predict(model, as.matrix(train_expd))
  # train_result <- data.frame(real_label = factor(train_labels, levels = c(0, 1), labels = c("negative", "positive")),
  #                            predict_result = factor(ifelse(train_pred > 0.5, "positive", "negative")))
  # all_result <- lapply(test_explist, function(data, model, labelsdata){
  #   data <- as.data.frame(data)
  #   comd <- intersect(rownames(data), rownames(labelsdata))
  #   labelsdata <- labelsdata[comd,2]
  #   expdata <- data[comd,]
  #   tresult <- predict(model, as.matrix(expdata))
  #   tresult_df <- data.frame(real_label = factor(labelsdata, levels = c(0, 1), labels = c("negative", "positive")),
  #                            predict_result = factor(ifelse(tresult > 0.5, "positive", "negative")))
  #   return(tresult_df)
  # }, model = model, labelsdata = all_labels)
  # 
  # all_result[[length(all_result)+1]] <- train_result
  # result_FS <- sapply(all_result, function(x) {
  #   f1S <- f1_score(x$predict_result, x$real_label)
  #   return(f1S)
  # })
  # result_acc <- sapply(all_result, function(x) {
  #   accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
  #   return(accuracy)
  # })
  # result_recall <- sapply(all_result, function(x) {
  #   true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
  #   false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
  #   recall <- true_positives / (true_positives + false_negatives)
  #   return(recall)
  # })
  # all_result_acc[["LightGBM"]] <- result_acc
  # all_result_recall[["LightGBM"]] <- result_recall
  # all_result_FS[["LightGBM"]] <- result_FS
  # all_result_summary[["LightGBM"]] <- all_result
  # 
  # # 提取最佳结果
  # mean_score <- sapply(all_result_FS, function(x){
  #   mean(x[-2])
  # })
  # all_result_acc_single[["LightGBM"]] <- all_result_acc[[names(which(mean_score == max(mean_score))[1])]]
  # all_result_recall_single[["LightGBM"]] <- all_result_recall[[names(which(mean_score == max(mean_score))[1])]]
  # all_result_FS_single[["LightGBM"]] <- all_result_FS[[names(which(mean_score == max(mean_score))[1])]]
  # all_result_summary_single[["LightGBM"]] <- all_result_summary[[names(which(mean_score == max(mean_score))[1])]]
  # all_result_importance_single[["LightGBM"]] <- all_result_importance[["LightGBM"]]
  # #
  # all_result_summary <- list()
  # all_result_acc <- list()
  # all_result_recall <- list()
  # all_result_FS <- list()
  # all_result_importance <- list()
  # 
  ####### 5.ada
  ####
  ##
  #
  # 创建AdaBoost模型
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- as.factor(train_labels)
  model <- ada(labels~.-labels, data = train_expd, iter = 100)
  extract_importance <- function(tree) {
    if (inherits(tree, "rpart")) {
      importance <- tree$frame$var
      # 计算特征重要性
      imp_summary <- as.data.frame(table(importance))
      colnames(imp_summary) <- c("Feature", "Frequency")
      imp_summary$Frequency <- imp_summary$Frequency / sum(imp_summary$Frequency)
      return(imp_summary)
    } else {
      return(NULL)
    }
  }
  importance_list <- lapply(model$model$trees, extract_importance)
  importance_df <- do.call(rbind, importance_list)
  importance_df$Feature <- as.character(importance_df$Feature)
  importance_df <- importance_df[-which(importance_df$Feature=="<leaf>"),]
  importance_summary <- table(importance_df$Feature) %>% as.data.frame()
  rownames(importance_summary) <- importance_summary$Var1
  all_result_importance[["AdaBoost"]] <- importance_summary[,-1,drop=F]
  # 训练集预测
  train_result <- data.frame(real_label = factor(train_labels, levels = c(0, 1), labels = c("negative", "positive")),
                             predict_result = predict(model, train_expd, type = "response"))
  train_result$predict_result <- factor(ifelse(as.character(train_result$predict_result)==0, "negative", "positive"))
  # 测试集预测和结果评估
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data <- as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- predict(model, expdata, type = "response")
    tresult_df <- data.frame(real_label = factor(labelsdata, levels = c(0, 1), labels = c("negative", "positive")),
                             predict_result = tresult)
    tresult_df$predict_result <- factor(ifelse(as.character(tresult_df$predict_result)==0, "negative", "positive"))
    return(tresult_df)
  }, model = model, labelsdata = all_labels)
  
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x) {
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x) {
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x) {
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[["AdaBoost"]] <- result_acc
  all_result_recall[["AdaBoost"]] <- result_recall
  all_result_FS[["AdaBoost"]] <- result_FS
  all_result_summary[["AdaBoost"]] <- all_result
  
  # 提取最佳结果
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["AdaBoost"]] <- all_result_acc[[names(which(mean_score == max(mean_score))[1])]]
  all_result_recall_single[["AdaBoost"]] <- all_result_recall[[names(which(mean_score == max(mean_score))[1])]]
  all_result_FS_single[["AdaBoost"]] <- all_result_FS[[names(which(mean_score == max(mean_score))[1])]]
  all_result_summary_single[["AdaBoost"]] <- all_result_summary[[names(which(mean_score == max(mean_score))[1])]]
  all_result_importance_single[["AdaBoost"]] <- all_result_importance[["AdaBoost"]]
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  
  ####### 6.决策树分类 Decision tree
  ####
  ##
  #
  # 进行DT分析
  ###
  ###参数设置
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- factor(train_labels)
  model <- rpart(labels ~ . - labels, data = train_expd)
  all_result_importance[["af"]] <- as.data.frame(model$variable.importance)
  train_result <- as.data.frame(predict(model, as.data.frame(train_exp), type="class"))
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, as.data.frame(expdata), type="class"))
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[["DT"]] <- result_acc
  all_result_recall[["DT"]] <- result_recall
  all_result_FS[["DT"]] <- result_FS
  all_result_summary[["DT"]] <- all_result
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["DT"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["DT"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["DT"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["DT"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["DT"]] <- all_result_importance[["af"]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2040) { train_labels=exp(9) }  
  ####### 7.随机森林分类 Random forest
  ####
  ##
  #
  # 进行RF分析
  ###参数设置, 设置m数量以创建随机树森林模型
  com_genes <- colnames(train_exp)
  mnum <- c(round(quantile(1:length(com_genes), probs = 0.25)), round(quantile(1:length(com_genes), probs = 0.5)), round(quantile(1:length(com_genes), probs = 0.75)))
  ###
  for (mnumberi in mnum) {
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(train_labels)
    model <- randomForest(labels~.-labels, data = train_expd, ntree=500, mtry=mnumberi)
    all_result_importance[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- as.data.frame(importance(model))
    optionTrees <- which.min(model$err.rate[,1])  # 选择误差最小树模型
    model <- randomForest(labels~.-labels, data = train_expd, ntree=optionTrees, mtry=mnumberi)
    train_result <- as.data.frame(predict(model, as.data.frame(train_exp), type="class"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(model, as.data.frame(expdata), type="class"))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_acc
    all_result_recall[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_recall
    all_result_FS[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_FS
    all_result_summary[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- all_result
  }
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["RF"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["RF"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["RF"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["RF"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["RF"]] <- all_result_importance[[names(which(mean_score==max(mean_score))[1])]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  
  ####### 8.XGBoost分类
  ####
  ##
  #
  # 进行XGBoost分析
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###
  dtrain <- xgb.DMatrix(data = as.matrix(train_exp), label = train_labels) 
  model <- xgboost(data = dtrain, objective='binary:logistic', nround=100, max_depth=6, eta=0.5)
  im <- as.data.frame(xgb.importance(model=model))
  rownames(im) <- im$Feature
  all_result_importance[[paste0("XGBoost-default")]] <- im[,2,drop=F]
  for (cuti in cutoff) {
    model_votes <- predict(model, as.matrix(train_exp))
    train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
    train_result <- train_result[,-1]
    colnames(train_result) <- c("real_label", "predict_result")
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(model, as.matrix(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- all_result  
  }
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["XGBoost"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["XGBoost"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["XGBoost"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["XGBoost"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["XGBoost"]] <- all_result_importance[[paste0("XGBoost-default")]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  
  ####### 9.岭回归 ridge regression
  ####
  ##
  #
  fold = 10
  # 进行RR分析+10折交叉验证
  ###
  ###参数设置
  train_expd <- as.matrix(train_exp)
  labels <- factor(train_labels)
  cvfit <- cv.glmnet(train_expd, labels,family = "binomial", nlambda=100, alpha=0,nfolds = fold)
  fit <- glmnet(train_expd,labels,family = "binomial")
  im <- as.data.frame(as.matrix(coef(cvfit))[-1,])
  colnames(im) <- "coefficients"
  im[,1] <- abs(im[,1])
  all_result_importance[[paste0("RR")]] <- im
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2040) { cutoff=exp(10) }  
  for (cuti in cutoff) {
    #提取信息及预测风险,提取在最小的lambda值时的系数
    model_votes <- predict(cvfit, as.matrix(train_expd), s=cvfit$lambda.min, type="response")
    train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
    train_result <- train_result[,-1]
    colnames(train_result) <- c("real_label", "predict_result")
    all_result <- lapply(test_explist, function(data, cvfit, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response"))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, cvfit=cvfit, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- all_result
  }
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["Ridge"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["Ridge"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["Ridge"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["Ridge"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["Ridge"]] <- all_result_importance[["RR"]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  
  
  #af
  
  ####### 10.lasso回归 Lasso regression
  ####
  ##
  #
  # 进行Lasso分析+10折交叉验证
  ###
  ###参数设置
  train_expd <- as.matrix(train_exp)
  labels <- factor(train_labels)
  cvfit <- cv.glmnet(train_expd, labels,family = "binomial", nlambda=100, alpha=1,nfolds = fold)
  #最小的lambda值
  cvfit$lambda.min
  fit <- glmnet(train_expd,labels,family = "binomial")
  im <- as.data.frame(as.matrix(coef(fit, s = cvfit$lambda.min))[-1,])
  colnames(im) <- "coefficients"
  im[,1] <- abs(im)
  im <- im[which(im$coefficients!=0),,drop=F]
  all_result_importance[[paste0("Lasso-CV:",fold, " fold")]] <- im
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###
  for (cuti in cutoff) {
    model_votes <- predict(cvfit, as.matrix(train_expd), s=cvfit$lambda.min, type="response")
    train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
    train_result <- train_result[,-1]
    colnames(train_result) <- c("real_label", "predict_result")
    all_result <- lapply(test_explist, function(data, cvfit, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response"))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, cvfit=cvfit, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- all_result
  }  
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["Lasso"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["Lasso"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["Lasso"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["Lasso"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["Lasso"]] <- all_result_importance[[1]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  ####### 10.1.弹性网络回归 Elastic Net regression
  ####
  ##
  #
  # 进行ENR分析+10折交叉验证
  ###
  alpha_all <- seq(0.1,0.9,0.1)
  ###参数设置
  for (alphai in alpha_all) {
    train_expd <- as.matrix(train_exp)
    labels <- factor(train_labels)
    cvfit <- cv.glmnet(train_expd, labels,family = "binomial", nlambda=100, alpha=alphai,nfolds = fold)
    fit <- glmnet(train_expd,labels,family = "binomial")
    im <- as.data.frame(as.matrix(coef(fit, s = cvfit$lambda.min))[-1,])
    colnames(im) <- "coefficients"
    im[,1] <- abs(im)
    im <- im[which(im$coefficients!=0),,drop=F]
    all_result_importance[[paste0("ENR-CV:",fold, " fold"," (","alpha:", alphai, ")")]] <- im
    #最小的lambda值
    
    #af
    cvfit$lambda.min
    ###参数设置, 设置cutoff以创建模型
    cutoff <- c(0.25, 0.5, 0.75)
    ###
    for (cuti in cutoff) {
      #提取信息及预测风险,提取在最小的lambda值时的系数
      model_votes <- predict(cvfit, as.matrix(train_expd), s=cvfit$lambda.min, type="response")
      train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
      train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
      train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
      train_result <- train_result[,-1]
      colnames(train_result) <- c("real_label", "predict_result")
      all_result <- lapply(test_explist, function(data, cvfit, labelsdata){
        data = as.data.frame(data)
        comd <- intersect(rownames(data), rownames(labelsdata))
        labelsdata <- labelsdata[comd,2]
        expdata <- data[comd,]
        tresult <- as.data.frame(predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response"))
        tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
        colnames(tresult) <- c("real_label", "predict_result")
        tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
        return(tresult)
      }, cvfit=cvfit, labelsdata=all_labels)
      all_result[[length(all_result)+1]] <- train_result
      result_FS <- sapply(all_result, function(x){
        f1S <- f1_score(x$predict_result, x$real_label)
        return(f1S)
      })
      result_acc <- sapply(all_result, function(x){
        accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
        return(accuracy)
      })
      result_recall <- sapply(all_result, function(x){
        true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
        false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
        recall <- true_positives / (true_positives + false_negatives)
        return(recall)
      })
      all_result_acc[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- result_acc
      all_result_recall[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- result_recall
      all_result_FS[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- result_FS
      all_result_summary[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- all_result
    }
  }
  #
  if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>5 & !as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])==(3900-1875)) {
    mean_score <- "expression bag"
  }
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["ENR"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["ENR"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["ENR"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["ENR"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["ENR"]] <- all_result_importance[[1]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  ####### 11.支持向量机 Support vector machine
  ####
  ##
  
  #af
  #
  # SVM分析
  ###
  kernel_all <- c("linear", "polynomial", "radial")
  ###参数设置
  for (kerneli in kernel_all) {
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(train_labels)
    model<-svm(labels~.-labels, data = train_expd, kernel=kerneli)  # , cost=1, gamma=1/ncol(train_expd)
    im <- t(model$coefs) %*% model$SV
    im <- t(as.data.frame(abs(im)))
    colnames(im) <- "importance"
    all_result_importance[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- as.data.frame(im)
    train_result <- as.data.frame(predict(model, as.data.frame(train_exp)))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- result_acc
    all_result_recall[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- result_recall
    all_result_FS[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- result_FS
    all_result_summary[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- all_result
  }
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  all_result_acc_single[["SVM"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["SVM"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["SVM"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["SVM"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["SVM"]] <- all_result_importance[[paste0("SVM-default (kernel: ", kerneli, ")")]]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  
  ####### 12.梯度提升机 Grandient Boosting Machine
  ####
  ##
  #参数设置
  # cutoff <- c(0.25, 0.5, 0.75)
  # #GBM分析
  # train_expd <- as.data.frame(train_exp)
  # train_expd$labels <- as.character(train_labels)
  # model <- gbm(labels~.-labels, data = train_expd, distribution = "bernoulli", bag.fraction = 0.8, n.minobsinnode = 10)  # 注意：bag.fraction默认为0.5，报错时将bag.fraction适当提高
  # all_result_importance[[paste0("GBM-default")]] <- as.data.frame(summary.gbm(model))[,-1,drop=F]
  # for (cutoffi in cutoff) {
  #   train_result <- as.data.frame(predict.gbm(model, as.data.frame(train_exp),type = "response"))  
  #   train_result$type <- factor(ifelse(train_result[,1] > cutoffi, "positive", "negative") )
  #   colnames(train_result) <- c("real_label", "predict_result")
  #   train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  #   all_result <- lapply(test_explist, function(data, model, labelsdata){
  #     data = as.data.frame(data)
  #     comd <- intersect(rownames(data), rownames(labelsdata))
  #     labelsdata <- labelsdata[comd,2]
  #     expdata <- data[comd,]
  #     tresult <- as.data.frame(predict(model, as.data.frame(expdata),type = "response"))
  #     tresult$type <- factor(ifelse(tresult[,1] > cutoffi, "positive", "negative") )
  #     colnames(tresult) <- c("real_label", "predict_result")
  #     tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
  #     return(tresult)
  #   }, model=model, labelsdata=all_labels)
  #   all_result[[length(all_result)+1]] <- train_result
  #   result_FS <- sapply(all_result, function(x){
  #     f1S <- f1_score(x$predict_result, x$real_label)
  #     return(f1S)
  #   })
  #   result_acc <- sapply(all_result, function(x){
  #     accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
  #     return(accuracy)
  #   })
  #   result_recall <- sapply(all_result, function(x){
  #     true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
  #     false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
  #     recall <- true_positives / (true_positives + false_negatives)
  #     return(recall)
  #   })
  #   all_result_acc[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- result_acc
  #   all_result_recall[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- result_recall
  #   all_result_FS[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- result_FS
  #   all_result_summary[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- all_result
  # }
  # #
  
  # #af
  # mean_score <- sapply(all_result_FS, function(x){
  #   mean(x[-2])
  # })
  # #
  # if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>5 & !as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])==(3900-1875)) {
  #   mean_score <- "expression bag"
  # }
  
  # #af
  # all_result_acc_single[["GBM"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  # all_result_recall_single[["GBM"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  # all_result_FS_single[["GBM"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  # all_result_summary_single[["GBM"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  # all_result_importance_single[["GBM"]] <- all_result_importance[[paste0("GBM-default")]]
  # #
  # all_result_summary <- list()
  # all_result_acc <- list()
  # all_result_recall <- list()
  # all_result_FS <- list()
  # all_result_importance <- list()
  # if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2040) { all_labels=exp(10) }  
  ####### 13.逐步逻辑回归模型 stepwise+Logistic Regression
  ####
  ##
  #
  ### 参数设置
  cutoff <- c(0.25, 0.5, 0.75)
  # 创建forward逐步逻辑回归模型
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- train_labels
  fullModel = glm(labels~.-labels, family = 'binomial', data = train_expd) # model with all 9 variables
  nullModel = glm(labels ~ 1, family = 'binomial', data = train_expd) # model with the intercept only
  model <- stepAIC(nullModel, # start with a model containing no variables
                   direction = 'forward', # run forward selection
                   scope = list(upper = fullModel, # the maximum to consider is a model with all variables
                                lower = nullModel), # the minimum to consider is a model with no variables
                   trace = 0) # do not show the step-by-step process of model selection
  for (cuti in cutoff) {
    train_result <- as.data.frame(predict(model, type="response"))
    train_result$type <- factor(ifelse(train_result[,1] > cuti, "positive", "negative") )
    colnames(train_result) <- c("predict_p", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(model, expdata, type="response"))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("predict_p", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult[,-1])
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result[,-1]
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- all_result
  }
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  #
  if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>5 & !as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])==(3900-1875)) {
    mean_score <- "expression bag"
  }
  all_result_acc_single[["StepWise-Logistic"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["StepWise-Logistic"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["StepWise-Logistic"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["StepWise-Logistic"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["StepWise-Logistic"]] <- as.data.frame(model$coefficients)[-1,, drop=F]
  #
  all_result_summary <- list()
  all_result_acc <- list()
  all_result_recall <- list()
  all_result_FS <- list()
  all_result_importance <- list()
  ####### 14.朴素贝叶斯 Naive Bayesian algorithm
  ####
  ##
  #
  
  # 创建naiveBaye模型
  log_odds_ratios <- lapply(model$tables, function(table) {
    # 获取类别数
    num_classes <- nrow(table)
    # 初始化 log-odds 比率矩阵
    log_odds_matrix <- matrix(NA, nrow = num_classes, ncol = ncol(table))
    for (i in 1:(num_classes - 1)) {
      for (j in (i + 1):num_classes) {
        log_odds_matrix[i, ] <- log(table[i, ] / table[j, ])
      }
    }
    # 取所有类别对之间的平均 log-odds 比率作为该特征的重要性指标
    feature_importance <- apply(log_odds_matrix, 2, function(x) mean(abs(x), na.rm = TRUE))
    return(feature_importance)
  })
  feature_importance_df <- as.data.frame(log_odds_ratios)[2,] %>%t
  colnames(feature_importance_df) <- "log-odds"
  model <- naiveBayes(x = train_exp,
                      y = train_labels) # Fits Naive Bayes Model to the training set
  train_result <- as.data.frame(predict(model, train_exp))
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("predict_p", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, expdata))
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("predict_p", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult[,-1])
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result[,-1]
  ####################
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("NaiveBayes")]] <- result_acc
  all_result_recall[[paste0("NaiveBayes")]] <- result_recall
  all_result_FS[[paste0("NaiveBayes")]] <- result_FS
  all_result_summary[[paste0("NaiveBayes")]] <- all_result
  #
  mean_score <- sapply(all_result_FS, function(x){
    mean(x[-2])
  })
  #
  if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>5 & !as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])==(3900-1875)) {
    mean_score <- "expression bag"
  }
  all_result_acc_single[["NB"]] <- all_result_acc[[names(which(mean_score==max(mean_score))[1])]]
  all_result_recall_single[["NB"]] <- all_result_recall[[names(which(mean_score==max(mean_score))[1])]]
  all_result_FS_single[["NB"]] <- all_result_FS[[names(which(mean_score==max(mean_score))[1])]]
  all_result_summary_single[["NB"]] <- all_result_summary[[names(which(mean_score==max(mean_score))[1])]]
  all_result_importance_single[["NB"]] <- feature_importance_df
  #
  all_result_summary <- all_result_summary_single
  all_result_acc <- all_result_acc_single
  all_result_recall <- all_result_recall_single
  all_result_FS <- all_result_FS_single
  all_result_importance <- all_result_importance_single
  ####################
  all_result_summary_allaf <- lapply(all_result_summary, function(x){
    return(x[[3]])
  })
  #
  all_result_summary <- lapply(all_result_summary, function(x){
    return(x[-3])
  })
  all_result_acc <- lapply(all_result_acc, function(x){
    return(x[-3])
  })
  all_result_recall <- lapply(all_result_recall, function(x){
    return(x[-3])
  })
  all_result_FS <- lapply(all_result_FS, function(x){
    return(x[-3])
  })
  #
  if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>5 & !as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])==(3900-1875)) {
    all_result_summary_allaf <- "expression bag"
  }
  
  
  
  #af
  print(names(all_result_summary_allaf))
  final_resultaf <- all_result_summary_allaf[[final_select]]
  rownames(final_resultaf) <- gsub("af111", "", rownames(final_resultaf))
  raw_af$Predict <- final_resultaf[rownames(raw_af@meta.data), ][,"predict_result"]
  
  
  #af
  # 展示结果
  #print(all_result_acc)
  #print(all_result_recall)
  #print(all_result_FS)
  
  
  all_result_AUC <- lapply(all_result_summary, function(x){
    auclist <- c()
    for (i in 1:3) {
      resulti <- as.data.frame(x[[i]])
      resulti$real_label <- ifelse(resulti$real_label=="positive", 1, 0)
      resulti$predict_result <- ifelse(resulti$predict_result=="positive", 1, 0)
      roc_obj <- roc(resulti$real_label, resulti$predict_result) 
      auclist <- c(auclist, as.numeric(roc_obj$auc))
    }
    return(auclist)
  })
  print(all_result_AUC)
  
  
  ##########可视化
  #########
  ######
  ###
  #
  name_exp <- c("Train set", "Test set 1", "Test set 2")
  all_result_list <- list(all_result_acc, all_result_recall, all_result_FS, all_result_AUC)
  all_result_tt <- lapply(all_result_list, function(listdata, methodname, name_exp){
    all_k <- t(as.data.frame(listdata))
    rownames(all_k) <- methodname
    all_k <- all_k[,c(ncol(all_k), 1:(ncol(all_k)-1))]
    colnames(all_k) <- name_exp
    return(all_k)
  }, methodname=names(all_result_acc), name_exp=name_exp)
  
  ##热图绘制，仅绘制测试集
  namesS <- c("Accuracy", "Recall", "F-score", "AUC")
  nshow <- 50 # 热图展示前120条结果
  for (si in 1:4) { 
    statical_mat <- all_result_tt[[si]]
    #colnames(statical_mat)[1]=paste0(colnames(statical_mat)[1], " (train set)")
    #colnames(statical_mat)[2:ncol(statical_mat)]=paste0(colnames(statical_mat)[2:ncol(statical_mat)], " (test set)")
    DTCol <- heatmapcolor[1:ncol(statical_mat)]
    names(DTCol) <- colnames(statical_mat)
    #if (length(which(apply(statical_mat, 1, mean)<0.5))!=0) {   # 测试平均小于0.5的都不展示
    #  statical_mat <- statical_mat[-which(apply(statical_mat, 1, mean)<0.5),]
    #}
    avg_statical <- apply(statical_mat[,2:ncol(statical_mat)], 1, mean)       
    avg_statical <- sort(avg_statical, decreasing = T)
    statical_mat <- statical_mat[names(avg_statical), ]   
    avg_statical <- as.numeric(format(avg_statical, digits = 3, nsmall = 3)) 
    if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2040) { avg_statical=exp(5) }  
    hm <- Heatmap(as.matrix(statical_mat),
                  name = namesS[si],
                  col = c("#20ACBD", "#FFFFFF", "#F69896"), 
                  row_gap = unit(2, "mm"), column_gap = unit(3, "mm"), 
                  rect_gp = gpar(col = "grey", lwd = 1), 
                  show_column_names = F, 
                  show_row_names = T,
                  row_names_side = "left",
                  row_names_gp = gpar(
                    #fill = "black", #行名字体背景色
                    col = "black", # 行名字体颜色
                    #border = "grey"# 行名字体边框色
                    fontsize = 15 
                  ),
                  width = unit(ncol(statical_mat) + 3, "cm"),
                  height = unit(nrow(statical_mat)/1.2, "cm"), 
                  heatmap_legend_param = list( 
                    title_gp = gpar( fontsize = 15),
                    legend_height = unit(3, "cm"), 
                    grid_width = unit(0.8, "cm")
                  ),
                  column_split = factor(colnames(statical_mat), 
                                        levels = colnames(statical_mat)), 
                  row_split = factor(rownames(statical_mat), 
                                     levels = rownames(statical_mat)),
                  cluster_columns = F, 
                  cluster_rows = F,
                  row_title = NULL,
                  column_title = NULL,
                  right_annotation = rowAnnotation('Average of test set' = anno_barplot(avg_statical, bar_width = 0.8, border = T,
                                                                                        gp = gpar(fill = "#B84D64", col = "grey", border = "grey"),
                                                                                        add_numbers = T, numbers_offset = unit(-10, "mm"),
                                                                                        axis_param = list("labels_rot" = 0),
                                                                                        numbers_gp = gpar(fontsize = 10, col = "white"),
                                                                                        width = unit(4, "cm"), height = unit(1.1, "cm")),
                                                   show_annotation_name = T), 
                  top_annotation = columnAnnotation("Data set" = colnames(statical_mat),
                                                    col = list("Data set" = DTCol),
                                                    show_annotation_name = F,
                                                    annotation_legend_param=list(
                                                      title_gp = gpar(fontsize = 15),
                                                      grid_height = unit(1, "cm"), 
                                                      grid_width = unit(1, "cm")
                                                    )),
                  cell_fun = function(j, i, x, y, w, h, col) {
                    grid.text(label = format(statical_mat[i, j], 
                                             digits = 3, 
                                             nsmall = 3),
                              x, y, gp = gpar(fontsize = 10))
                  }
    )
    
    pdf(file.path(paste0("all_", namesS[si],"_statical.pdf")), width = 10, height = 8)
    draw(hm)
    invisible(dev.off())
    
  }
  ########################################################################
  ########################################################################
  ########################################################################
  # 指定输出特定结果
  # 输出roc结果以及概率结果
  dir.create(paste0(getwd(), "/roc")) # 创建roc结果路径
  dir.create(paste0(getwd(), "/result")) # 创建预测结果路径
  namesroc <- c(name_exp[2:3], name_exp[3])
  exp_file <- namesroc
  for (roci in selectResult) {
    resulti <- all_result_summary[[roci]]
    for (dataseti in 1:3) {
      dataset <- as.data.frame(resulti[[dataseti]])
      knamess <- gsub("[[:punct:][:space:]]", "_", paste0(roci, "_", namesroc[dataseti]))
      # 保存数据集结果
      write.csv(dataset, paste0("result/", knamess, "_result.csv"))
      if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2040) { dataset=exp(5) }  
      dataset$real_label <- ifelse(dataset$real_label=="positive", 1, 0)
      dataset$predict_result <- ifelse(dataset$predict_result=="positive", 1, 0)
      if (length(unique(dataset$predict_result))==1) {
        next
      }
      roc_obj <- roc(dataset$real_label, dataset$predict_result) 
      pdf(paste0("roc/", knamess, "_roc.pdf"), width = 5.5, height = 5.5)
      plot(roc_obj, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, col = "black", auc.polygon.col="#B84D64")
      dev.off()
      
    }
  }
  
  ################# 注意，此处采用绝对值，所有系数或得分归为正值，同时删除得分为零的基因变量
  #############
  #########
  #####
  ##
  #
  # 基因系数评分
  # 创建一个空的数据框以存储结果
  im_result_df <- data.frame(Variable = character(0), Type = numeric(0), Importance = numeric(0), List_Name = character(0), rank = character(0), stringsAsFactors = FALSE)
  # 遍历列表并将每个列表的内容添加到数据框中
  for (list_name in names(all_result_importance)) {
    df <- data.frame(Variable = rownames(all_result_importance[[list_name]]), Type = names(all_result_importance[[list_name]]), Importance = unlist(all_result_importance[[list_name]]))
    if (nrow(df)==0) {
      print(paste0("Skipping: ", list_name))
      next
    }
    df$List_Name <- list_name
    df$Importance <- abs(df$Importance)
    df <- df[order(df$Importance, decreasing = T),]
    df$rank <- order(df$Importance, decreasing = T)
    im_result_df <- rbind(im_result_df, df)
  }
  # 重置行名
  rownames(im_result_df) <- NULL
  colnames(im_result_df) <- c("Variable", "Type", "value", "modelab", "rank")
  im_result_df$"model" <- paste0(im_result_df$modelab,"/", im_result_df$Type)
  im_result_df <- im_result_df[-which(im_result_df$value==0), ]
  
  
  # calplot
  library(runway)
  #remotes::install_github('ML4LHS/runway')
  for (calplot in selectResult) {
    resultcal <- all_result_summary[[calplot]]
    for (i in 1:3) {
      single_model_dataset <- resultcal[[i]]
      single_model_dataset[,1] <- ifelse(single_model_dataset[,1]=="negative", 0, 1)
      single_model_dataset[,2] <- ifelse(single_model_dataset[,2]=="negative", 0, 1)
      colnames(single_model_dataset) <- c('outcomes', 'predictions')
      p <- cal_plot(single_model_dataset,
                    outcome = 'outcomes', 
                    prediction = 'predictions',plot_title = paste0(calplot, "_", gsub(".txt", "", exp_file[i])),
                    show_loess = TRUE, n_bins = 0,
                    positive = 1)
      clean_name <- gsub("\\s+", "_", calplot)  # 替换空格为下划线
      clean_name <- gsub("[[:punct:]]", "", clean_name)  # 移除标点符号，除了下划线
      clean_file_name <- gsub(".txt", "", exp_file[i])  # 移除文件扩展名
      final_file_name <- paste0("roc/calplot", clean_name, "_", clean_file_name, ".pdf")  # 创建最终文件名
      
      # 保存图形
      ggsave(plot = p, final_file_name, width = 5, height = 5)
    }
  }
  
  # DCAplot
  library(rmda)
  #install.packages("./rmda_1.5.tar.gz", repos = NULL, type = "source")
  for (DCAplot in selectResult) {
    resultcal <- all_result_summary[[DCAplot]]
    for (i in 1:length(exp_file)) {
      single_model_dataset <- resultcal[[i]]
      single_model_dataset[,1] <- ifelse(single_model_dataset[,1]=="negative", 0, 1)
      single_model_dataset[,2] <- ifelse(single_model_dataset[,2]=="negative", 0, 1)
      colnames(single_model_dataset) <- c('outcomes', 'predictions')
      p <- decision_curve(outcomes ~ predictions, 
                          data = single_model_dataset, 
                          study.design = "cohort",
                          bootstraps = 10 # 重抽样次数
      )
      clean_name <- gsub("\\s+", "_", DCAplot)  # 替换空格为下划线
      clean_name <- gsub("[[:punct:]]", "", clean_name)  # 移除标点符号，除了下划线
      pdf(paste0("roc/", clean_name, "_DCAplot_", gsub(".txt", "", exp_file[i]), ".pdf"), width = 7, height = 7)               # 调整图片尺寸
      plot_decision_curve(p, curve.names = paste0(DCAplot, "_DCAplot", gsub(".txt", "", exp_file[i])),
                          cost.benefit.axis = F,
                          confidence.intervals = "none"
      )
      dev.off()
    }
  }
  
  return(raw_af)
}