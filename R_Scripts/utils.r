# Given a data frame (df) returns the same df without the class columns
factors.remove <- function (df)
{
  factors <- names(Filter(is.factor, df))
  return(df[,which(names(df) != factors)])
}

factors.get <- function (df)
{
  factors <- names(Filter(is.factor, df))
  return(df[,which(names(df) == factors)])
}

auc_all <- function (df)
{
  predictor = factors.get(df)
  
  level_count = length(levels(predictor))
  
  df_names <- c()
  df_aucs <- c()
  
  for(name in names(factors.remove(df))) {
    if (level_count > 2) {
      r <- pROC::multiclass.roc(predictor, df[,name])
      a <- as.numeric(r$auc)  
    }else {
      a <- pROC::auc(predictor, df[,name])
    }
    df_names <- c(df_names, name)
    df_aucs <- c(df_aucs, a)
  }
  
  return(data.frame(names=df_names, aucs=df_aucs))
}

svm_rfe__single_step <- function(X, y, ...) {
  
  library(kernlab)
  
  # keep track of the iteration during which
  # a feature was eliminated
  ii <- rep(NA, ncol(X))
  i <- 0
  while ( any(is.na(ii)) ) {
    # indices of remaining features
    not_elim_yet <- which(is.na(ii))
    # train the classifier on the remaining features
    fit <- ksvm(X[,not_elim_yet], y, scaled = FALSE, ...)
    # compute the feature weights
    sv_i <- alphaindex(fit)[[1]]
    w <- t( coef(fit)[[1]] ) %*% X[ sv_i, not_elim_yet ]
    # eliminate the feature with the smallest squared weight
    to_elim <- not_elim_yet[ head(order( w * w ), 1) ]
    ii[to_elim] <- i
    i <- i + 1
  }
  # convert iterations into ranks
  print (i)
  print (ii)
  i - ii
}

calculate.confusion <- function(states, clusters)
{
  # generate a confusion matrix of cols C versus states S
  d <- data.frame(state = states, cluster = clusters)
  td <- as.data.frame(table(d))
  # convert from raw counts to percentage of each label
  pc <- matrix(ncol=max(clusters),nrow=0) # k cols
  for (i in 1:9) # 9 labels
  {
    total <- sum(td[td$state==td$state[i],3])
    pc <- rbind(pc, td[td$state==td$state[i],3]/total)
  }
  rownames(pc) <- td[1:9,1]
  return(pc)
}

assign.cluster.labels <- function(cm, k)
{
  # take the cluster label from the highest percentage in that column
  cluster.labels <- list()
  for (i in 1:k)
  {
    cluster.labels <- rbind(cluster.labels, row.names(cm)[match(max(cm[,i]), cm[,i])])
  }
  
  # this may still miss some labels, so make sure all labels are included
  for (l in rownames(cm)) 
  { 
    if (!(l %in% cluster.labels)) 
    { 
      cluster.number <- match(max(cm[l,]), cm[l,])
      cluster.labels[[cluster.number]] <- c(cluster.labels[[cluster.number]], l)
    } 
  }
  return(cluster.labels)
}

calculate.accuracy <- function(states, clabels)
{
  matching <- Map(function(state, labels) { state %in% labels }, states, clabels)
  tf <- unlist(matching, use.names=FALSE)
  return (sum(tf)/length(tf))
}

# Remove Outliers from Dataframe using Grubbs test
grubbs.RemoveOutliersOnDataframe <- function(x, factorCol) {
    for (name in colnames(x)) {
        test <- x[,c(name)]
        grubbs.result <- grubbs.test(test)
        pv <- grubbs.result$p.value
        if (is.na(pv)) {
            pv <- 1
        }
        while(pv < 0.05) {
            outlierValue <- as.numeric(strsplit(grubbs.result$alternative," ")[[1]][3])
            roundedtest <- round(test, 6)
            indexOfOutlier <- which(roundedtest ==  round(outlierValue, 6))[1]
            x <- x[-c(indexOfOutlier),]
            if(!missing(factorCol)) {
              factorCol <- factorCol[-c(indexOfOutlier)]
            } 
            # else {
            #   x + y
            # }
            test <- x[,c(name)]
            grubbs.result <- grubbs.test(test)
            pv <- grubbs.result$p.value
            # p-value goes to 0 when all the values are equal
            if (is.na(pv) || (sum(abs(diff(test))) == 0)) {
                pv <- 1
            }
        }
    }
    return (cbind(factorCol,x))
}

calc_ROC_2Classes <- function(probabilities, known_truth, model.name=NULL)
{
  library(dplyr)
  outcome <- as.numeric(factor(known_truth))-1
  pos <- sum(outcome) # total known positives
  neg <- sum(1-outcome) # total known negatives
  pos_probs <- outcome*probabilities # probabilities for known positives
  neg_probs <- (1-outcome)*probabilities # probabilities for known negatives
  true_pos <- sapply(probabilities,
                     function(x) sum(pos_probs>=x)/pos) # true pos. rate
  false_pos <- sapply(probabilities,
                      function(x) sum(neg_probs>=x)/neg)
  if (is.null(model.name))
    result <- data.frame(true_pos, false_pos)
  else
    result <- data.frame(true_pos, false_pos, model.name)
  result %>% arrange(false_pos, true_pos)
}

####################################
## Hack to return multiple values ##
####################################
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}

# Function that recieves dataframe and select the specified number of features using a specific feature selection method
selectFeaturesWithFSMethods <- function(data_train, data_train_whc, number_features, fs) { #whc means with highly correlated as this is multivariate and should deal with those relations
  if (fs == 'FSCR') { # Fisher score (FSCR)
    library(PredPsych)
    fisher.score <- fscore(data_train, 1, seq(2,ncol(data_train)), silent = FALSE)
    df_fisher_score <- data.frame(fisher.score)
    colnames(df_fisher_score) <- c("fscore")
    df_fisher_score <- df_fisher_score[order(df_fisher_score$fscore, decreasing = TRUE), , drop = FALSE]
    namesOfFeature <- rownames(df_fisher_score)
    selectedFeatures <- namesOfFeature[seq(1, number_features)]
    
    label_selectedFeatures <- c("Label", selectedFeatures)
    joint_df <- data_train[ , label_selectedFeatures]
    
  }
  else if (fs == 'RELF') { # Relief (RELF)
    library(CORElearn)
    relief_score <- CORElearn::attrEval(LABEL_COL, data_train, "Relief")
    df_relief_score <- data.frame(relief_score)
    colnames(df_relief_score) <- c("reliefscore")
    df_relief_score <- df_relief_score[order(df_relief_score$reliefscore, decreasing = TRUE), , drop = FALSE]
    namesOfFeature <- rownames(df_relief_score)
    selectedFeatures <- namesOfFeature[seq(1, number_features)]
    
    label_selectedFeatures <- c("Label", selectedFeatures)
    joint_df <- data_train[ , label_selectedFeatures]
    
  }
  else if (fs == 'TSCR') { # T-test (t-score) (TSCR)
    library(stats)
    group <- data_train$Label
    t_test_score <- lapply(data_train[,-1], function(x) { stats::t.test(x ~ group)$statistic })
    t_test_score_unlist <- unlist(t_test_score,use.names=TRUE)
    df_t_test_score <- data.frame(matrix(t_test_score_unlist))
    rownames(df_t_test_score) <- c(colnames(data_train[,-1]))
    colnames(df_t_test_score) <- c("ttestscore")
    df_t_test_score <- df_t_test_score[order(df_t_test_score$ttestscore, decreasing = TRUE), , drop = FALSE]
    namesOfFeature <- rownames(df_t_test_score)
    selectedFeatures <- namesOfFeature[seq(1, number_features)]
    
    label_selectedFeatures <- c("Label", selectedFeatures)
    joint_df <- data_train[ , label_selectedFeatures]
    
  }
  else if (fs == 'CHSQ') { # Chi-square (CHSQ)
    group <- data_train$Label
    chi_square_score <- lapply(data_train[,-1], function(x) { chisq.test(cbind(-1 * min(x) + x, group))$statistic })
    chi_square_score_unlist <- unlist(chi_square_score,use.names=TRUE)
    df_chi_square_score <- data.frame(matrix(chi_square_score_unlist))
    rownames(df_chi_square_score) <- c(colnames(data_train[,-1]))
    colnames(df_chi_square_score) <- c("chisquaretestscore")
    df_chi_square_score <- df_chi_square_score[order(df_chi_square_score$chisquaretestscore, decreasing = TRUE), , drop = FALSE]
    namesOfFeature <- rownames(df_chi_square_score)
    selectedFeatures <- namesOfFeature[seq(1, number_features)]
    
    label_selectedFeatures <- c("Label", selectedFeatures)
    joint_df <- data_train[ , label_selectedFeatures]
    
  }
  else if (fs == 'WLCX') { # Wilcoxon (WLCX)
    group <- data_train$Label
    wilcox_score <- lapply(data_train[,-1], function(x) { wilcox.test( x, as.numeric(group))$statistic })
    wilcox_score_unlist <- unlist(wilcox_score,use.names=TRUE)
    df_wilcox_score <- data.frame(matrix(wilcox_score_unlist))
    rownames(df_wilcox_score) <- c(colnames(data_train[,-1]))
    colnames(df_wilcox_score) <- c("wilcoxtestscore")
    df_wilcox_score <- df_wilcox_score[order(df_wilcox_score$wilcoxtestscore, decreasing = FALSE), , drop = FALSE]
    namesOfFeature <- rownames(df_wilcox_score)
    selectedFeatures <- namesOfFeature[seq(1, number_features)]
    
    label_selectedFeatures <- c("Label", selectedFeatures)
    joint_df <- data_train[ , label_selectedFeatures]
    
  }
  else if (fs == 'GINI') { # Gini index (GINI)
    library(CORElearn)
    gini_index_score <- CORElearn::attrEval(LABEL_COL, data_train, "Gini")
    df_gini_index_score <- data.frame(gini_index_score)
    colnames(df_gini_index_score) <- c("giniscore")
    df_gini_index_score <- df_gini_index_score[order(df_gini_index_score$giniscore, decreasing = TRUE), , drop = FALSE]
    namesOfFeature <- rownames(df_gini_index_score)
    selectedFeatures <- namesOfFeature[seq(1, number_features)]
    
    label_selectedFeatures <- c("Label", selectedFeatures)
    joint_df <- data_train[ , label_selectedFeatures]
    
  }
  else if (fs == 'MIM') { # Mutual information maximization (MIM)
    library(entropy)
    group <- data_train$Label
    mim_score <- lapply(data_train[,-1], function(x) { mi.empirical(discretize2d(x, as.numeric(group), numBins1=16, numBins2=16)) })
    mim_score_unlist <- unlist(mim_score, use.names=TRUE)
    df_mim_score <- data.frame(matrix(mim_score_unlist))
    rownames(df_mim_score) <- c(colnames(data_train[,-1]))
    colnames(df_mim_score) <- c("mimscore")
    df_mim_score <- df_mim_score[order(df_mim_score$mimscore, decreasing = TRUE), , drop = FALSE]
    namesOfFeature <- rownames(df_mim_score)
    selectedFeatures <- namesOfFeature[seq(1, number_features)]
    
    label_selectedFeatures <- c("Label", selectedFeatures)
    joint_df <- data_train[ , label_selectedFeatures]
  }
  else if (fs == 'mRMR') { # Minimum redundancy maximum relevance (MRMR)
    # MINIMUM REDUNDANCY MAXIMUM RELEVANCE
    library(mRMRe)
    df_asNumeric <- data_train_whc
    df_asNumeric$Label <- as.numeric(df_asNumeric$Label)
    dd <- mRMR.data(data = df_asNumeric)#[ , -c(LABEL_COL)]
    # For mRMR.classic-like results
    # Number of features MRMR
    num_features <-  number_features
    ensemble <- mRMR.ensemble(data = dd, target_indices = c(1), solution_count = 1, feature_count = num_features)
    causality(ensemble)
    cols_mrmr <- unlist(solutions(ensemble))
    df_filtered <- df_asNumeric[, c(cols_mrmr)]
    
    factorCol <- as.factor(data_train_whc[, LABEL_COL])
    joint_df <- cbind(data_train_whc[, LABEL_COL], df_filtered)
    colnames(joint_df) <- c('Label', colnames(df_filtered))
    #joint_df$Label <- as.factor(data_train_whc[, LABEL_COL])
    
  }
  else if (fs == 'JMI') { # Joint mutual information (JMI)
    library(praznik)
    group <- data_train_whc$Label
    df_jmi <- praznik::JMIM(data_train_whc[,-1], group, k =  number_features) #, threads = 0
    
    selectedFeaturesJMI <- names(df_jmi$score)
    
    label_selectedFeaturesJMI <- c("Label", selectedFeaturesJMI)
    joint_df <- data_train_whc[ , label_selectedFeaturesJMI]
    
  }
  return (list('df'=joint_df, 'fs_method'=fs, 'n_features'=number_features))
}

############################################################################
## https://johanndejong.wordpress.com/2018/04/04/nested-cross-validation/ ##
############################################################################

glm_train_and_validate <- function(data, data_withHighCorr, fold, n_feat, fs) {
  library(stats)
  library(PRROC)
  list_data_fs_fold_fs_nfeat <- selectFeaturesWithFSMethods(data_train = data[-fold,], data_train_whc = data_withHighCorr[-fold,], number_features = n_feat, fs=fs)
  data_fs_fold <- list_data_fs_fold_fs_nfeat$df
  fs_method <- list_data_fs_fold_fs_nfeat$fs_method
  n_features <- list_data_fs_fold_fs_nfeat$n_features
  #c(data_fs_fold, fs_method, n_features) := selectFeaturesWithFSMethods(data_train = data[-fold,], data_train_whc = data_withHighCorr[-fold,], number_features = n_feat, fs=fs)
  cols_fs <- colnames(data_fs_fold)
  # Train an SVM, excluding the fold
  fit <- stats::glm(Label ~ ., family = binomial(), data = data_fs_fold#,#data[-fold,],
              )
  # Predict the fold
  if (fs != 'mRMR' && fs != 'JMI') {
    yh <- predict(fit, newdata = data[fold,cols_fs], type = "response")#probabilities
    # Compare the predictions to the labels
    posneg <- split(yh, data$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  } else {
    yh <- predict(fit, newdata = data_withHighCorr[fold,cols_fs], type = "response")#probabilities
    # Compare the predictions to the labels
    posneg <- split(yh, data_withHighCorr$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  }
  # Return the AUC under the ROC
  PRROC::roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv_glm <- function( data, dataWHC, k, NumFeaturesList, FeatureSelectionList, seed = NULL) {
  # Set the seed, if given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  list_NF_FS <- list(nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_NF_FS <- do.call(expand.grid, list_NF_FS)
  
  # For each value of the hyperparameter C ...
  auc <- mapply(function(n_feat,fs) {
    folds <- createFolds(data$Label, k = k)
    # For each fold ...
    
    sapply(folds, function(fold) {
      # Train an SVM, and validate on the fold
      glm_train_and_validate(data, dataWHC, fold, n_feat, fs)
    })
  }, n_feat=grid_search_NF_FS$nFS, fs=grid_search_NF_FS$FS)
  auc
}

ncv_glm <- function(data, dataWHC, k, NumFeaturesList, FeatureSelectionList, n_reps = 1, seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  list_NF_FS <- list(nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_NF_FS <- do.call(expand.grid, list_NF_FS)
  
  folds <- createFolds(data$Label, k = k)
  # For each fold ...
  auc <- replicate(n_reps, {
    sapply(folds, function(fold) {
      # Do a cross-validation for each C
      auc <- cv_glm(data[-fold,], dataWHC[-fold,], k, NumFeaturesList, FeatureSelectionList, seed = seed)
      # Select the C with the highest AUC
      max_nFS <- grid_search_NF_FS$nFS[which.max(apply(auc, 2, mean))]
      max_FS <- grid_search_NF_FS$FS[which.max(apply(auc, 2, mean))]
      
      glm_train_and_validate(data, dataWHC, fold, max_nFS, max_FS)
    })
  })
  auc
}

lda_train_and_validate <- function(data, data_withHighCorr, fold, n_feat, fs) {
  library(MASS)
  library(PRROC)
  data_fs_fold <- selectFeaturesWithFSMethods(data_train = data[-fold,], data_train_whc = data_withHighCorr[-fold,], number_features = n_feat, fs=fs)
  cols_fs <- colnames(data_fs_fold$df)
  # Train an SVM, excluding the fold
  data_df <- data_fs_fold$df
  fit <- MASS::lda(Label ~ ., data = data_df)
  # Predict the fold
  if (fs != 'mRMR' && fs != 'JMI') {
    yh <- predict(fit, newdata = data[fold,cols_fs], type = "probabilities")
    # Compare the predictions to the labels
    posneg <- split(yh$posterior[,1], data$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  } else {
    yh <- predict(fit, newdata = data_withHighCorr[fold,cols_fs], type = "probabilities")#probabilities
    # Compare the predictions to the labels
    posneg <- split(yh$posterior[,1], data_withHighCorr$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  }
  # yh <- predict(fit, newdata = data[fold,cols_fs], type = "probabilities")# predict(fit, newdata = data[fold,cols_fs], type = "probabilities")
  # # Compare the predictions to the labels
  # posneg <- split(yh$posterior[,1], data$Label[fold])
  # Return the AUC under the ROC
  PRROC::roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv_lda <- function( data, dataWHC, k, NumFeaturesList, FeatureSelectionList, seed = NULL) {
  # Set the seed, if given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  list_NF_FS <- list(nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_NF_FS <- do.call(expand.grid, list_NF_FS)
  
  # For each value of the hyperparameter C ...
  auc <- mapply(function(n_feat,fs) {
    folds <- createFolds(data$Label, k = k)
    # For each fold ...
    sapply(folds, function(fold) {
      # Train an SVM, and validate on the fold
      lda_train_and_validate(data, dataWHC, fold, n_feat, fs)
    })
  }, n_feat=grid_search_NF_FS$nFS, fs=grid_search_NF_FS$FS)
  auc
}

ncv_lda <- function(data, dataWHC, k, NumFeaturesList, FeatureSelectionList, n_reps = 1, seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  list_NF_FS <- list(nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_NF_FS <- do.call(expand.grid, list_NF_FS)
  
  folds <- createFolds(data$Label, k = k)
  # For each fold ...
  auc <- replicate(n_reps, {
    sapply(folds, function(fold) {
    # Do a cross-validation for each C
    auc <- cv_lda(data[-fold,], dataWHC[-fold,], k, NumFeaturesList, FeatureSelectionList, seed = seed)
    # Select the C with the highest AUC
    max_nFS <- grid_search_NF_FS$nFS[which.max(apply(auc, 2, mean))]
    max_FS <- grid_search_NF_FS$FS[which.max(apply(auc, 2, mean))]
    
    lda_train_and_validate(data, dataWHC, fold, max_nFS, max_FS)
  })
  })
  auc
}

glmnet_train_and_validate <- function(data, data_withHighCorr, alpha, lambda, fold, n_feat, fs) {
  library(glmnet)
  library(PRROC)
  data_train_k_1_fold <- data[-fold,]
  data_train_whc_k_1_fold <- data_withHighCorr[-fold,]
  
  preprocessParams_k_1_fold <- caret::preProcess(data_train_k_1_fold, method=c("center", "scale"), na.remove=TRUE, verbose=FALSE);
  data_train_k_1_fold <- stats::predict(preprocessParams_k_1_fold, data_train_k_1_fold)
  
  preprocessParams_whc_k_1_fold <- caret::preProcess(data_train_whc_k_1_fold, method=c("center", "scale"), na.remove=TRUE, verbose=FALSE);
  data_train_whc_k_1_fold <- stats::predict(preprocessParams_whc_k_1_fold, data_train_whc_k_1_fold)
  
  data_fs_fold <- selectFeaturesWithFSMethods(data_train = data_train_k_1_fold, data_train_whc = data_train_whc_k_1_fold, number_features = n_feat, fs=fs)
  cols_fs <- colnames(data_fs_fold$df)
  data_df <- data_fs_fold$df
  # Train an SVM, excluding the fold
  fit <- glmnet::glmnet(as.matrix(data_df[,-c(1)]), as.matrix(data_df[,c(1)]), family=c("binomial"), alpha = alpha, lambda = lambda)
  # Predict the fold
  if (fs != 'mRMR' && fs != 'JMI') {
    data_nhc <- data[fold,]
    data_nhc <- stats::predict(preprocessParams_k_1_fold, data_nhc)
    data_pred <- data_nhc[,cols_fs]
    # data_pred <- stats::predict(preprocessParams_k_1_fold, data_pred)
    
    yh <- predict(fit, newx = as.matrix(data_pred[,-c(1)]), type = "response")
    # Compare the predictions to the labels
    posneg <- split(yh[,1], data$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  } else {
    data_whc <- data_withHighCorr[fold,]
    data_whc <- stats::predict(preprocessParams_whc_k_1_fold, data_whc)
    data_whc_pred <- data_whc[,cols_fs]
    # 
    # data_whc_pred <- stats::predict(preprocessParams_whc_k_1_fold, data_whc_pred)
    # 
    yh <- predict(fit, newx = as.matrix(data_whc_pred[,-c(1)]), type = "response")#probabilities
    # Compare the predictions to the labels
    posneg <- split(yh[,1], data_withHighCorr$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  }
  # Return the AUC under the ROC
  PRROC::roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv_glmnet <- function( data, dataWHC, k, AlphaList, LambdaList, NumFeaturesList, FeatureSelectionList, seed = NULL) {
  # Set the seed, if given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  list_Alpha_Lambda_NF_FS <- list(alpha = AlphaList, lambda = LambdaList, nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_Alpha_Lambda_NF_FS <- do.call(expand.grid, list_Alpha_Lambda_NF_FS)
  
  # For each value of the hyperparameter C ...
  auc <- mapply(function(alpha, lambda, n_feat, fs) {
    folds <- createFolds(data$Label, k = k)
    print(paste0("alpha: ",alpha))
    print(paste0("lambda: ",lambda))
    print(paste0("n_feat: ",n_feat))
    print(paste0("fs: ",fs))
    # For each fold ...
    sapply(folds, function(fold) {
      # Train an GLMNET, and validate on the fold
      glmnet_train_and_validate(data, dataWHC, alpha, lambda, fold, n_feat, fs)
    })
  }, alpha=grid_search_Alpha_Lambda_NF_FS$alpha, lambda=grid_search_Alpha_Lambda_NF_FS$lambda, n_feat=grid_search_Alpha_Lambda_NF_FS$nFS, fs=grid_search_Alpha_Lambda_NF_FS$FS)
  auc
}

ncv_glmnet <- function(data, dataWHC, k, AlphaList, LambdaList, NumFeaturesList, FeatureSelectionList, n_reps = 1, seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  list_HP_NF_FS <- list(alpha = AlphaList, lambda = LambdaList, nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_HP_NF_FS <- do.call(expand.grid, list_HP_NF_FS)
  
  folds <- createFolds(data$Label, k = k)
  # For each fold ...
  auc <- replicate(n_reps, {
    sapply(folds, function(fold) {
      # Do a cross-validation for each C
      auc <- cv_glmnet(data[-fold,], dataWHC[-fold,], k, AlphaList, LambdaList, NumFeaturesList, FeatureSelectionList, seed = seed)
      # Select the C with the highest AUC
      max_alpha <- grid_search_HP_NF_FS$alpha[which.max(apply(auc, 2, mean))]
      max_lambda <- grid_search_HP_NF_FS$lambda[which.max(apply(auc, 2, mean))]
      max_nFS <- grid_search_HP_NF_FS$nFS[which.max(apply(auc, 2, mean))]
      max_FS <- grid_search_HP_NF_FS$FS[which.max(apply(auc, 2, mean))]
      
      glmnet_train_and_validate(data, dataWHC, max_alpha, max_lambda, fold, max_nFS, max_FS)
    })
  })
  auc
}

# KNN
# fit.knn <- caret::train(frm, data=joint_df, method="knn", metric=metric, trControl=control)
# knn(train, test, cl, k = 1, l = 0, prob = tru, use.all = TRUE)
knn_train_and_validate <- function(data, data_withHighCorr, k_knn, fold, n_feat, fs) {
  library(glmnet)
  library(PRROC)
  data_fs_fold <- selectFeaturesWithFSMethods(data_train = data[-fold,], data_train_whc = data_withHighCorr[-fold,], number_features = n_feat, fs=fs)
  cols_fs <- colnames(data_fs_fold$df)
  data_df <- data_fs_fold$df
  # Train an SVM, excluding the fold
  # fit <- class::knn(data_df[,-c(1)], data[fold,cols_fs], cl=data[-fold,]$Label, k_knn, prob = TRUE)
  # data_test <- data[fold,cols_fs]
  if (fs != 'mRMR' && fs != 'JMI') {
    data_test <- data[fold,cols_fs]
  } else {
    data_test <- data_withHighCorr[fold,cols_fs]
  }
  fit <- class::knn(data_df[,-c(1)], data_test[,-c(1)], cl=data_df[,c(1)], k_knn, prob = TRUE)
  # Predict the fold
  # yh <- predict(fit, newdata = data[fold,cols_fs], type = "probabilities")
  #yh <- predict(fit, newdata = data_test[,-c(1)], type = "probabilities")
  prob <- attr(fit, "prob")
  posneg <- split(prob, data$Label[fold])
  # Compare the predictions to the labels
  # posneg <- split(yh[,1], data$Label[fold])
  # Return the AUC under the ROC
  PRROC::roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv_knn <- function( data, dataWHC, k, KList, NumFeaturesList, FeatureSelectionList, seed = NULL) {
  # Set the seed, if given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  list_K_NF_FS <- list(k_knn = KList, nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_K_NF_FS <- do.call(expand.grid, list_K_NF_FS)
  
  # For each value of the hyperparameter C ...
  auc <- mapply(function(k_knn, n_feat,fs) {
    folds <- createFolds(data$Label, k = k)
    # For each fold ...
    sapply(folds, function(fold) {
      # Train an GLMNET, and validate on the fold
      knn_train_and_validate(data, dataWHC, k_knn, fold, n_feat, fs)
    })
  }, k_knn=grid_search_K_NF_FS$k_knn, n_feat=grid_search_K_NF_FS$nFS, fs=grid_search_K_NF_FS$FS)
  auc
}

ncv_knn <- function(data, dataWHC, k, KList, NumFeaturesList, FeatureSelectionList, n_reps = 1, seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  list_HP_NF_FS <- list(k_knn = KList, nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_HP_NF_FS <- do.call(expand.grid, list_HP_NF_FS)
  
  folds <- createFolds(data$Label, k = k)
  # For each fold ...
  auc <- replicate(n_reps, {
    sapply(folds, function(fold) {
    # Do a cross-validation for each C
    auc <- cv_knn(data[-fold,], dataWHC[-fold,], k, KList, NumFeaturesList, FeatureSelectionList, seed = seed)
    # Select the C with the highest AUC
    max_k <- grid_search_HP_NF_FS$k_knn[which.max(apply(auc, 2, mean))]
    max_nFS <- grid_search_HP_NF_FS$nFS[which.max(apply(auc, 2, mean))]
    max_FS <- grid_search_HP_NF_FS$FS[which.max(apply(auc, 2, mean))]
    
    knn_train_and_validate(data, dataWHC, max_k, fold, max_nFS, max_FS)
  })
  })
  auc
}

# Naive Bayes
# fit.nb <- caret::train(frm, data=joint_df, method="nb", metric=metric, trControl=control)
nb_train_and_validate <- function(data, data_withHighCorr, fold, n_feat, fs) {
  library(e1071)
  library(PRROC)
  data_fs_fold <- selectFeaturesWithFSMethods(data_train = data[-fold,], data_train_whc = data_withHighCorr[-fold,], number_features = n_feat, fs=fs)
  cols_fs <- colnames(data_fs_fold$df)
  data_df <- data_fs_fold$df
  # Train an SVM, excluding the fold
  data_df$Label <- factor(data_df$Label)
  fit <- e1071::naiveBayes(Label~., data=data_df) #knn(data_fs_fold, data[fold,cols_fs], cl=data_fs_fold$Label, k_knn, prob = TRUE)
  # Predict the fold
  if (fs != 'mRMR' && fs != 'JMI') {
    data$Label 
    yh <- predict(fit, newdata = data[fold,cols_fs], type = "raw")
    # Compare the predictions to the labels
    posneg <- split(yh[,1], data$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  } else {
    yh <- predict(fit, newdata = data_withHighCorr[fold,cols_fs], type = "raw")#probabilities
    # Compare the predictions to the labels
    posneg <- split(yh[,1], data_withHighCorr$Label[fold])#posneg <- split(yh[1,], data$Label[fold])
  }
  # yh <- predict(fit, newdata = data[fold,cols_fs], type = "probabilities")
  # # Compare the predictions to the labels
  # posneg <- split(yh[,1], data$Label[fold])
  # Return the AUC under the ROC
  PRROC::roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv_nb <- function( data, dataWHC, k, NumFeaturesList, FeatureSelectionList, seed = NULL) {
  # Set the seed, if given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  list_NF_FS <- list(nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_NF_FS <- do.call(expand.grid, list_NF_FS)
  
  # For each value of the hyperparameter C ...
  auc <- mapply(function(n_feat,fs) {
    folds <- createFolds(data$Label, k = k)
    # For each fold ...
    sapply(folds, function(fold) {
      # Train an GLMNET, and validate on the fold
      nb_train_and_validate(data, dataWHC, fold, n_feat, fs)
    })
  }, n_feat=grid_search_NF_FS$nFS, fs=grid_search_NF_FS$FS)
  auc
}

ncv_nb <- function(data, dataWHC, k, NumFeaturesList, FeatureSelectionList, n_reps = 1, seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  list_NF_FS <- list(nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_NF_FS <- do.call(expand.grid, list_NF_FS)
  
  folds <- createFolds(data$Label, k = k)
  # For each fold ...
  auc <- replicate(n_reps, {
    sapply(folds, function(fold) {
      # Do a cross-validation for each C
      auc <- cv_nb(data[-fold,], dataWHC[-fold,], k, NumFeaturesList, FeatureSelectionList, seed = seed)
      # Select the C with the highest AUC
      max_nFS <- grid_search_NF_FS$nFS[which.max(apply(auc, 2, mean))]
      max_FS <- grid_search_NF_FS$FS[which.max(apply(auc, 2, mean))]
      
      nb_train_and_validate(data, dataWHC, fold, max_nFS, max_FS)
    })
  })
  auc
}

# Function for one round of training and validating an SVM
ksvm_train_and_validate <- function(data, data_withHighCorr, fold, C, sigma, n_feat, fs) {
  library(kernlab)
  data_fs_fold <- selectFeaturesWithFSMethods(data_train = data[-fold,], data_train_whc = data_withHighCorr[-fold,], number_features = n_feat, fs=fs)
  cols_fs <- colnames(data_fs_fold)
  # Train an SVM, excluding the fold
  fit <- kernlab::ksvm(Label ~ ., data = data_fs_fold,#data[-fold,],
                       kernel = "rbfdot", kpar = list(sigma=sigma), C = C, prob.model = TRUE, Label.weights = 1 / table(data$Label[-fold]))
  # Predict the fold
  yh <- predict(fit, newdata = data[fold,cols_fs], type = "probabilities")
  # Compare the predictions to the labels
  posneg <- split(yh[,1], data$Label[fold])
  # Return the AUC under the ROC
  PRROC::roc.curve(posneg[[1]], posneg[[2]])$auc
}

# Function for doing a k-fold cross-validation for each C in CC
cv_ksvm <- function( data, dataWHC, k, CC, SigmaList, NumFeaturesList, FeatureSelectionList, seed = NULL) {
  # Set the seed, if given
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  list_HP_NF_FS <- list(hpC = CC, hpS = SigmaList, nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_HP_NF_FS <- do.call(expand.grid, list_HP_NF_FS)
  
  # For each value of the hyperparameter C ...
  auc <- mapply(function(C, sigma,n_feat,fs) {
    folds <- createFolds(data$Label, k = k)
    # For each fold ...
    sapply(folds, function(fold) {
      # Train an SVM, and validate on the fold
      ksm_train_and_validate(data, dataWHC, fold, C, sigma, n_feat, fs)
    })
  }, C=grid_search_HP_NF_FS$hpC, sigma=grid_search_HP_NF_FS$hpS, n_feat=grid_search_HP_NF_FS$nFS, fs=grid_search_HP_NF_FS$FS)
  
  auc
}

ncv_svm <- function(data, dataWHC, k, CC, SigmaList, NumFeaturesList, FeatureSelectionList, n_reps = 1, seed = 123) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  list_HP_NF_FS <- list(hpC = CC, hpS = SigmaList, nFS = NumFeaturesList, FS = FeatureSelectionList)
  grid_search_HP_NF_FS <- do.call(expand.grid, list_HP_NF_FS)
  folds <- createFolds(data$Label, k = k)
  # For each fold ...
  auc <- replicate(n_reps, {
    sapply(folds, function(fold) {
      # Do a cross-validation for each C
      auc <- cv_ksvm(data[-fold,], dataWHC[-fold,], k, CC, SigmaList, NumFeaturesList, FeatureSelectionList, seed = seed)
      # Select the C with the highest AUC
      indx_max_C <- grid_search_HP_NF_FS$hpC
      indx_max_S <- grid_search_HP_NF_FS$hpS
      indx_max_nFS <- grid_search_HP_NF_FS$nFS
      indx_max_FS <- grid_search_HP_NF_FS$FS
      
      index_max_overall <- Reduce(intersect, list(indx_max_C, indx_max_S, indx_max_nFS, indx_max_FS))
      
      C <- grid_search_HP_NF_FS$hpC[index_max_overall]
      sigma <- grid_search_HP_NF_FS$hpS[index_max_overall]
      n_feat <- grid_search_HP_NF_FS$nFS[index_max_overall]
      fs <- grid_search_HP_NF_FS$FS[index_max_overall]
      # C <- CC[which.max(sapply(auc, mean))]
      # sigma <- SigmaList[which.max(sapply(auc, mean))]
      # n_feat <- NumFeaturesList[which.max(sapply(auc, mean))]
      # fs <- FeatureSelectionList[which.max(sapply(auc, mean))]
      # C  1, sample(C, 1), C)
      # Test this C on the test data
      ksvm_train_and_validate(data, dataWHC, fold, C, sigma, n_feat, fs)
    })
  })
  auc
}

results.model.cv <- function(model, roc_cv, ci_cv) {
  fold_indx_cv <- as.numeric(rownames(model$bestTune))
  
  auc_cv <- model$results$AUC[fold_indx_cv]
  kappa_cv <- model$results$Kappa[fold_indx_cv]
  acc_cv <- model$results$Accuracy[fold_indx_cv]
  
  auc_cv_SD <- model$results$AUCSD[fold_indx_cv]
  kappa_cv_SD <- model$results$KappaSD[fold_indx_cv]
  acc_cv_SD <- model$results$AccuracySD[fold_indx_cv]
  
  if (length(model$results$F1) != 0) {
    f1_cv <- model$results$F1[fold_indx_cv]
    sens_cv <- model$results$Sensitivity[fold_indx_cv]
    spec_cv <- model$results$Specificity[fold_indx_cv]
    ppv_cv <- model$results$Pos_Pred_Value[fold_indx_cv]
    npv_cv <- model$results$Neg_Pred_Value[fold_indx_cv]
    
    f1_cv_SD <- model$results$F1SD[fold_indx_cv]
    sens_cv_SD <- model$results$SensitivitySD[fold_indx_cv]
    spec_cv_SD <- model$results$SpecificitySD[fold_indx_cv]
    ppv_cv_SD <- model$results$Pos_Pred_ValueSD[fold_indx_cv]
    npv_cv_SD <- model$results$Neg_Pred_ValueSD[fold_indx_cv]
  } else {
    f1_cv <- model$results$Mean_F1[fold_indx_cv]
    sens_cv <- model$results$Mean_Sensitivity[fold_indx_cv]
    spec_cv <- model$results$Mean_Specificity[fold_indx_cv]
    ppv_cv <- model$results$Mean_Pos_Pred_Value[fold_indx_cv]
    npv_cv <- model$results$Mean_Neg_Pred_Value[fold_indx_cv]
    
    f1_cv_SD <- model$results$Mean_F1SD[fold_indx_cv]
    sens_cv_SD <- model$results$Mean_SensitivitySD[fold_indx_cv]
    spec_cv_SD <- model$results$Mean_SpecificitySD[fold_indx_cv]
    ppv_cv_SD <- model$results$Mean_Pos_Pred_ValueSD[fold_indx_cv]
    npv_cv_SD <- model$results$Mean_Neg_Pred_ValueSD[fold_indx_cv]
  }
  
  cv_results <- matrix(nrow = 8, ncol = 3, dimnames = list(c('AUC', 'Kappa', 'F1', 'Acc', 'Sens', 'Spec', 'PPV', 'NPV'), c('Mean', '95CI_Left', '95CI_Right')))
  error_auc_cv <- qnorm(0.975)*auc_cv_SD/sqrt(nrow(dataset)/3); left_auc_cv <- auc_cv-error_auc_cv; right_auc_cv <- auc_cv+error_auc_cv
  error_kappa_cv <- qnorm(0.975)*kappa_cv_SD/sqrt(nrow(dataset)/3); left_kappa_cv <- kappa_cv-error_kappa_cv; right_kappa_cv <- kappa_cv+error_kappa_cv
  error_f1_cv <- qnorm(0.975)*f1_cv_SD/sqrt(nrow(dataset)/3); left_f1_cv <- f1_cv-error_f1_cv; right_f1_cv <- f1_cv+error_f1_cv
  error_acc_cv <- qnorm(0.975)*acc_cv_SD/sqrt(nrow(dataset)/3); left_acc_cv <- acc_cv-error_acc_cv; right_acc_cv <- acc_cv+error_acc_cv
  error_sens_cv <- qnorm(0.975)*sens_cv_SD/sqrt(nrow(dataset)/3); left_sens_cv <- sens_cv-error_sens_cv; right_sens_cv <- sens_cv+error_sens_cv
  error_spec_cv <- qnorm(0.975)*spec_cv_SD/sqrt(nrow(dataset)/3); left_spec_cv <- spec_cv-error_spec_cv; right_spec_cv <- spec_cv+error_spec_cv
  error_ppv_cv <- qnorm(0.975)*ppv_cv_SD/sqrt(nrow(dataset)/3); left_ppv_cv <- ppv_cv-error_ppv_cv; right_ppv_cv <- ppv_cv+error_ppv_cv
  error_npv_cv <- qnorm(0.975)*npv_cv_SD/sqrt(nrow(dataset)/3); left_npv_cv <- npv_cv-error_npv_cv; right_npv_cv <- npv_cv+error_npv_cv
  
  # cv_results[1,]<-c(roc_cv$auc, ci_cv[1], ci_cv[3])
  cv_results[1,]<-c(auc_cv, left_auc_cv, right_auc_cv)
  cv_results[2,]<-c(kappa_cv, left_kappa_cv, right_kappa_cv)
  # if (length(c(f1_cv, left_f1_cv, right_f1_cv))!=0) {
  #   cv_results[3,]<-c(f1_cv, left_f1_cv, right_f1_cv)
  # }
  cv_results[3,]<-c(f1_cv, left_f1_cv, right_f1_cv)
  cv_results[4,]<-c(acc_cv, left_acc_cv, right_acc_cv)
  cv_results[5,]<-c(sens_cv, left_sens_cv, right_sens_cv)
  cv_results[6,]<-c(spec_cv, left_spec_cv, right_spec_cv)
  cv_results[7,]<-c(ppv_cv, left_ppv_cv, right_ppv_cv)
  cv_results[8,]<-c(npv_cv, left_npv_cv, right_npv_cv)
  
  return (cv_results)
}

multiClassSummaryOwn <- function(data, lev = NULL, model = NULL){
  #Check data
  if (!all(levels(data[, "pred"]) == levels(data[, "obs"])))
    stop("levels of observed and predicted data do not match")
  
  has_class_probs <- all(lev %in% colnames(data))
  
  if(has_class_probs) {
    ## Overall multinomial loss
    lloss <- mnLogLoss(data = data, lev = lev, model = model)
    
    caret::: requireNamespaceQuietStop("ModelMetrics")
    caret::: requireNamespaceQuietStop("MLmetrics")
    
    #Calculate custom one-vs-all ROC curves for each class
    prob_stats <-
      lapply(levels(data[, "pred"]),
             function(x){
               #Grab one-vs-all data for the class
               obs  <- ifelse(data[,"obs"] == x, 1, 0)
               prob <- data[,x]
               roc_auc <- try(ModelMetrics::auc(obs, data[,x]), silent = TRUE)
               # browser()
               pr_auc <- try(MLmetrics::PRAUC(y_pred = data[,x], 
                                              y_true = obs),
                             silent = TRUE)
               if(inherits(pr_auc, "try-error"))
                 pr_auc <- NA
               res <- c(ROC = roc_auc, AUC = pr_auc)
               return(res)
             })
    prob_stats <- do.call("rbind", prob_stats)
    prob_stats <- colMeans(prob_stats, na.rm = TRUE)
  }
  
  #Calculate confusion matrix-based statistics
  CM <- confusionMatrix(data[, "pred"], data[, "obs"],
                        mode = "everything", positive = positive_class)
  print(positive_class)
  
  #Aggregate and average class-wise stats
  #Todo: add weights
  # RES: support two classes here as well
  if (length(levels(data[, "pred"])) == 2) {
    class_stats <- CM$byClass
  } else {
    class_stats <- colMeans(CM$byClass)
    names(class_stats) <- paste("Mean", names(class_stats))
  }
  
  # Aggregate overall stats
  overall_stats <-
    if (has_class_probs)
      c(CM$overall,
        logLoss = as.numeric(lloss),
        AUC = unname(prob_stats["ROC"]),
        prAUC = unname(prob_stats["AUC"]))
  else
    CM$overall
  
  
  # Combine overall with class-wise stats and remove some stats we don't want
  stats <- c(overall_stats, class_stats)
  stats <- stats[! names(stats) %in% c('AccuracyNull', "AccuracyLower", "AccuracyUpper",
                                       "AccuracyPValue", "McnemarPValue",
                                       'Mean Prevalence', 'Mean Detection Prevalence')]
  
  # Clean names
  names(stats) <- gsub('[[:blank:]]+', '_', names(stats))
  
  # Change name ordering to place most useful first
  # May want to remove some of these eventually
  stat_list <- c("Accuracy", "Kappa", "Mean_F1", 
                 "Mean_Sensitivity", "Mean_Specificity",
                 "Mean_Pos_Pred_Value", "Mean_Neg_Pred_Value", 
                 "Mean_Precision", "Mean_Recall",
                 "Mean_Detection_Rate",
                 "Mean_Balanced_Accuracy")
  if(has_class_probs) stat_list <- c("logLoss", "AUC", "prAUC", stat_list)
  if (length(levels(data[, "pred"])) == 2) stat_list <- gsub("^Mean_", "", stat_list)
  
  stats <- stats[c(stat_list)]
  
  return(stats)
}
