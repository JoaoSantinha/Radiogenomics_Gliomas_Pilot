###############################################################
n_repeats <- 10
n_folds <- 3
highlyCorrelatedCutoff <- 1.1 #0.95 #
n_lambdas <- 100
min_lambda <- -10
max_lambda <- 3
positive_class <- "WT"
metric <- "MCC"
split_train <- 0.7
# label string sometimes change when we request the TCGA data from biolinks but its one of these two if you get an error around rows_na_label
label_string <- "subtype_IDH.status"
label_string <- "paper_IDH.status"
###############################################################

seed_value <- 9 
set.seed(seed_value)

# Paths vars
UTILS_PATHFILE <- "utils.r"

# Path radiomics file of GBMs and LGGs
# addapt this to your case
radiomics_file <- '<complete_dir_of_repository_folder>/Radiomics_File/radiomics_Glioma_All_Seqs_All_Labels_BW.csv'

# csv column separator
SEP <- ","

# Load doMC libraries
sysinf <- Sys.info()
if (!is.null(sysinf)){
  os <- sysinf['sysname']
  if (os == 'Darwin') {
    library(doMC)
    registerDoMC(cores=1)
  }
}

# Make debug easier 
options(show.error.locations = TRUE)
options(error=traceback)

# Helper functions
source(UTILS_PATHFILE)

# Load data from TCGA
library(TCGAbiolinks)
if (!exists("gbm_met_data_df")) {
  query.gbm.met <- GDCquery(project = "TCGA-GBM",
                            data.category = "DNA Methylation",
                            legacy = FALSE,
                            platform = c("Illumina Human Methylation 450"))
  
  
  gbm_met_data_df <- colDataPrepare(getResults(query.gbm.met)$cases)
}

if (!exists("lgg_met_data_df")) {
  query.lgg.met <- GDCquery(project = "TCGA-LGG",
                            data.category = "DNA Methylation",
                            legacy = FALSE,
                            platform = c("Illumina Human Methylation 450"))
  
  lgg_met_data_df <- colDataPrepare(getResults(query.lgg.met)$cases)
}

# read file with data

radiomics_df_hgg_lgg <- read.csv(radiomics_file, header=TRUE, sep=SEP)

radiomics_df_hgg_lgg <- radiomics_df_hgg_lgg[!duplicated(radiomics_df_hgg_lgg$ID), ]

clinical_data_hgg_lgg <- rbind(gbm_met_data_df, lgg_met_data_df)

clinical_data_hgg_lgg_primary_tumor <- clinical_data_hgg_lgg[clinical_data_hgg_lgg$shortLetterCode == "TP", ]

ids <- sort(intersect(radiomics_df_hgg_lgg$ID, clinical_data_hgg_lgg_primary_tumor$patient)) 

hgg_lgg_df_with_radiomics <- clinical_data_hgg_lgg_primary_tumor[ids, ]
radiomics_df_hgg_lgg_sorted_clinical <- radiomics_df_hgg_lgg[radiomics_df_hgg_lgg$ID %in% ids, ]
radiomics_df_hgg_lgg_sorted_clinical <- radiomics_df_hgg_lgg_sorted_clinical[order(radiomics_df_hgg_lgg_sorted_clinical$ID), ]

library(caret)
library(InvariantCausalPrediction)

# Get environments
ExpInd <- as.numeric(as.factor(as.numeric(hgg_lgg_df_with_radiomics$subtype_Tissue.source.site)))
ExpInd <- as.numeric(as.factor(as.numeric(hgg_lgg_df_with_radiomics$paper_Tissue.source.site)))

radiomics_df_hgg_lgg_sorted_clinical_noID <- radiomics_df_hgg_lgg_sorted_clinical[, -c(1)]

# Remove ZV, NZV, Center and Scale data
preprocessParams <- caret::preProcess(radiomics_df_hgg_lgg_sorted_clinical_noID, method=c("center", "scale", "zv", "nzv"), na.remove=TRUE, verbose=TRUE) # , "center", "scale"
trans <- stats::predict(preprocessParams, radiomics_df_hgg_lgg_sorted_clinical_noID)

rows_na_label <- which(is.na(hgg_lgg_df_with_radiomics[,label_string]))#$subtype_IDH.status))

if (length(rows_na_label) != 0) {
  X <- trans[-c(rows_na_label),]
  Y <- make.names(hgg_lgg_df_with_radiomics[-c(rows_na_label),label_string])
  ExpInd_final <- ExpInd[-c(rows_na_label)]
} else {
  X <- trans
  Y <- make.names(hgg_lgg_df_with_radiomics[,label_string])
  ExpInd_final <- ExpInd
}

# icp_lasso <- ICP(data.matrix(X), Y, ExpInd_final, alpha = 0.05, selection = "lasso")#"stability")#"all")#"boosting")#
# print(icp_lasso)
# plot(icp_lasso)
# 
# icp_stability <- ICP(data.matrix(X), Y, ExpInd_final, alpha = 0.05, selection = "stability")#"all")#"boosting")#"lasso")#
# print(icp_stability)
# plot(icp_stability)
# 
# icp_boosting <- ICP(data.matrix(X), Y, ExpInd_final, alpha = 0.05, selection = "boosting")#"lasso")#"stability")#"all")#
# print(icp_boosting)
# plot(icp_boosting)
# 
# library(nonlinearICP)
# nonlinearICP_results <- nonlinearICP(data.matrix(X), Y, ExpInd_final, maxSizeSets = 20)

multiClassSummaryOwnMCC <- function(data, lev = NULL, model = NULL){
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
  
  MCC <- (CM$table[2,2] * CM$table[1,1] - CM$table[2,1] * CM$table[1,2])/
    (sqrt((CM$table[2,2] + CM$table[2,1]) * 
            (CM$table[2,2] + CM$table[1,2]) * 
            (CM$table[1,1] + CM$table[2,1]) * 
            (CM$table[1,1] + CM$table[1,2])))
  previous_names <- names(stats)
  stats <- c(MCC, stats)
  names(stats) <- c("MCC", previous_names)
  
  # Change name ordering to place most useful first
  # May want to remove some of these eventually
  stat_list <- c("MCC", "Accuracy", "Kappa", "Mean_F1", 
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

# Get Results
train_test_results <- list()
for (i in sort(unique(ExpInd_final))) {
  indx_test <- ExpInd_final == i
  
  X_train <- X[!indx_test,]
  X_test <- X[indx_test,]
  
  Y_train <- Y[!indx_test]
  Y_test <- Y[indx_test]
  
  ExpInd_train <- as.numeric(as.factor(ExpInd_final[!indx_test]))
  
  icp_lasso <- ICP(data.matrix(X_train), Y_train, ExpInd_train, alpha = 0.05, selection = "lasso", showAcceptedSets = FALSE, showCompletion = FALSE)
  print(icp_lasso)
  
  icp_stability <- ICP(data.matrix(X_train), Y_train, ExpInd_train, alpha = 0.05, selection = "stability", showAcceptedSets = FALSE, showCompletion = FALSE)
  print(icp_stability)
  
  icp_boosting <- ICP(data.matrix(X_train), Y_train, ExpInd_train, alpha = 0.05, selection = "boosting", showAcceptedSets = FALSE, showCompletion = FALSE)
  print(icp_boosting)
  
  cv_results <- list()
  coefs_model <- list()
  test_results <- list()
  for (typeCausalityAnalysis in c("lasso", "stability", "boosting", "none")) { #
    LABEL_COL <- "Label"
    if (typeCausalityAnalysis == "lasso") {
      computeBool <- !icp_lasso$modelReject
      trans_red <- cbind(Y_train, as.data.frame(X_train[, icp_lasso$usedvariables]))
      # trans_red <- cbind(Y_train, as.data.frame(X_selected))
      colnames(trans_red) <- c(LABEL_COL, colnames(X_train[, icp_lasso$usedvariables]))
    } else if (typeCausalityAnalysis == "stability") {
      computeBool <- !icp_stability$modelReject
      trans_red <- cbind(Y_train, as.data.frame(X_train[, icp_stability$usedvariables]))
      colnames(trans_red) <- c(LABEL_COL, colnames(X_train[, icp_stability$usedvariables]))
    } else if (typeCausalityAnalysis == "boosting") {
      computeBool <- !icp_boosting$modelReject
      trans_red <- cbind(Y_train, as.data.frame(X_train[, icp_boosting$usedvariables]))
      colnames(trans_red) <- c(LABEL_COL, colnames(X_train[, icp_boosting$usedvariables]))
    } else {
      computeBool <- TRUE
      correlationMatrix <- stats::cor(X_train)
      highlyCorrelated <- caret::findCorrelation(abs(correlationMatrix), cutoff=highlyCorrelatedCutoff)
      if (length(highlyCorrelated) != 0) {
        lowlyCorrelated <- colnames(X_train)[-highlyCorrelated]
        # lowlyCorrelated <- append(LABEL_COL, lowlyCorrelated)
      } else {
        lowlyCorrelated <- colnames(X_train)
      }
      
      trans_red <- cbind(Y_train, as.data.frame(X_train[,lowlyCorrelated]))
      colnames(trans_red) <- c(LABEL_COL, colnames(X_train[,lowlyCorrelated]))
    }
    if (computeBool) {
    
      frm = as.formula(paste(LABEL_COL, "~", "."))
      control <- caret::trainControl(method="repeatedcv", number=n_folds, repeats = n_repeats, classProbs = TRUE, savePredictions=TRUE, summaryFunction=multiClassSummaryOwnMCC)#multiClassSummaryOwn) #seeds=seeds, index=createFolds(trans_red$Label), 
      
      fit.glmnet_lasso_train_test <- caret::train(frm, data=trans_red, preProcess = c("center", "scale"), method="glmnet", tuneGrid = data.frame(alpha = 1, lambda = 2^runif(n_lambdas, min = min_lambda, max_lambda)), tuneLength = n_lambdas, metric=metric, trControl=control)
      
      ## just so error doesn't fail
      dataset <- X_train
      
      library(pROC)
      roc_radiomics_birads_cv <- roc(fit.glmnet_lasso_train_test$pred$obs[fit.glmnet_lasso_train_test$pred$lambda == fit.glmnet_lasso_train_test$bestTune$lambda], fit.glmnet_lasso_train_test$pred[fit.glmnet_lasso_train_test$pred$lambda == fit.glmnet_lasso_train_test$bestTune$lambda, make.names(positive_class)])
      roc_radiomics_birads_cv$auc
      ci_radiomics_birads_cv <- ci.auc(roc_radiomics_birads_cv)
      
      ciobj_rb <- ci.se(roc_radiomics_birads_cv, specificities=seq(0, 1, l=length(roc_radiomics_birads_cv$sensitivities)))
      dat.ci_rb <- data.frame(x = as.numeric(rownames(ciobj_rb)),
                              lower = ciobj_rb[, 1],
                              upper = ciobj_rb[, 3])
      
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      n = 2
      cols = gg_color_hue(n)
      
      # Radiomics + BIRADS
      p_radiomics_birads <- ggroc(list("Radiomics CV" = roc_radiomics_birads_cv), aes = c("colour"), title = "ROC curve", legacy.axes = FALSE) + 
        geom_ribbon(data=dat.ci_rb, aes(x = x, ymin = lower, ymax = upper), alpha = 0.3, inherit.aes = FALSE,fill = cols[1])+
        geom_line() + labs(x = "Specificity", y = "Sensitivity") + labs(colour = "Models") +
        annotate("text", x = .35, y = .5, color=cols[1], size = 6, label = paste("AUC =", round(roc_radiomics_birads_cv$auc, 2), ' [', round(ci_radiomics_birads_cv[1], 2), ',', round(ci_radiomics_birads_cv[3], 2), ']', sep = '')) +
        theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12), 
              axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
              legend.title = element_text(size = 12), legend.text = element_text(size = 10))
      print(p_radiomics_birads)
      
      fold_indx_radiomics_birads_cv <- as.numeric(rownames(fit.glmnet_lasso_train_test$bestTune))
      
      aucs_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$AUC[fold_indx_radiomics_birads_cv]
      kappa_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$Kappa[fold_indx_radiomics_birads_cv]
      f1_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$F1[fold_indx_radiomics_birads_cv]
      mcc_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$MCC[fold_indx_radiomics_birads_cv]
      acc_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$Accuracy[fold_indx_radiomics_birads_cv]
      sens_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$Sensitivity[fold_indx_radiomics_birads_cv]
      spec_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$Specificity[fold_indx_radiomics_birads_cv]
      ppv_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$Pos_Pred_Value[fold_indx_radiomics_birads_cv]
      npv_radiomics_birads_cv <- fit.glmnet_lasso_train_test$results$Neg_Pred_Value[fold_indx_radiomics_birads_cv]
      
      aucs_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$AUCSD[fold_indx_radiomics_birads_cv]
      kappa_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$KappaSD[fold_indx_radiomics_birads_cv]
      f1_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$F1SD[fold_indx_radiomics_birads_cv]
      mcc_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$MCCSD[fold_indx_radiomics_birads_cv]
      acc_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$AccuracySD[fold_indx_radiomics_birads_cv]
      sens_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$SensitivitySD[fold_indx_radiomics_birads_cv]
      spec_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$SpecificitySD[fold_indx_radiomics_birads_cv]
      ppv_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$Pos_Pred_ValueSD[fold_indx_radiomics_birads_cv]
      npv_radiomics_birads_cv_SD <- fit.glmnet_lasso_train_test$results$Neg_Pred_ValueSD[fold_indx_radiomics_birads_cv]
      
      
      cv_radiomics_birads <- matrix(nrow = 9, ncol = 3, dimnames = list(c('AUC', 'Kappa', 'F1', 'Acc', 'Sens', 'Spec', 'PPV', 'NPV', 'MCC'), c('Mean', '95CI_Left', '95CI_Right')))
      print('CV Radiomics')
      error_kappa_radiomics_birads_cv <- qnorm(0.975)*kappa_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_kappa_radiomics_birads_cv <- kappa_radiomics_birads_cv-error_kappa_radiomics_birads_cv; right_kappa_radiomics_birads_cv <- kappa_radiomics_birads_cv+error_kappa_radiomics_birads_cv
      error_f1_radiomics_birads_cv <- qnorm(0.975)*f1_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_f1_radiomics_birads_cv <- f1_radiomics_birads_cv-error_f1_radiomics_birads_cv; right_f1_radiomics_birads_cv <- f1_radiomics_birads_cv+error_f1_radiomics_birads_cv
      error_mcc_radiomics_birads_cv <- qnorm(0.975)*mcc_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_mcc_radiomics_birads_cv <- mcc_radiomics_birads_cv-error_mcc_radiomics_birads_cv; right_mcc_radiomics_birads_cv <- mcc_radiomics_birads_cv+error_mcc_radiomics_birads_cv
      error_acc_radiomics_birads_cv <- qnorm(0.975)*acc_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_acc_radiomics_birads_cv <- acc_radiomics_birads_cv-error_acc_radiomics_birads_cv; right_acc_radiomics_birads_cv <- acc_radiomics_birads_cv+error_acc_radiomics_birads_cv
      error_sens_radiomics_birads_cv <- qnorm(0.975)*sens_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_sens_radiomics_birads_cv <- sens_radiomics_birads_cv-error_sens_radiomics_birads_cv; right_sens_radiomics_birads_cv <- sens_radiomics_birads_cv+error_sens_radiomics_birads_cv
      error_spec_radiomics_birads_cv <- qnorm(0.975)*spec_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_spec_radiomics_birads_cv <- spec_radiomics_birads_cv-error_spec_radiomics_birads_cv; right_spec_radiomics_birads_cv <- spec_radiomics_birads_cv+error_spec_radiomics_birads_cv
      error_ppv_radiomics_birads_cv <- qnorm(0.975)*ppv_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_ppv_radiomics_birads_cv <- ppv_radiomics_birads_cv-error_ppv_radiomics_birads_cv; right_ppv_radiomics_birads_cv <- ppv_radiomics_birads_cv+error_ppv_radiomics_birads_cv
      error_npv_radiomics_birads_cv <- qnorm(0.975)*npv_radiomics_birads_cv_SD/sqrt(nrow(dataset)/3); left_npv_radiomics_birads_cv <- npv_radiomics_birads_cv-error_npv_radiomics_birads_cv; right_npv_radiomics_birads_cv <- npv_radiomics_birads_cv+error_npv_radiomics_birads_cv
      
      cv_radiomics_birads[1,]<-c(roc_radiomics_birads_cv$auc, ci_radiomics_birads_cv[1], ci_radiomics_birads_cv[3])
      cv_radiomics_birads[2,]<-c(kappa_radiomics_birads_cv, left_kappa_radiomics_birads_cv, right_kappa_radiomics_birads_cv)
      cv_radiomics_birads[3,]<-c(f1_radiomics_birads_cv, left_f1_radiomics_birads_cv, right_f1_radiomics_birads_cv)
      cv_radiomics_birads[4,]<-c(acc_radiomics_birads_cv, left_acc_radiomics_birads_cv, right_acc_radiomics_birads_cv)
      cv_radiomics_birads[5,]<-c(sens_radiomics_birads_cv, left_sens_radiomics_birads_cv, right_sens_radiomics_birads_cv)
      cv_radiomics_birads[6,]<-c(spec_radiomics_birads_cv, left_spec_radiomics_birads_cv, right_spec_radiomics_birads_cv)
      cv_radiomics_birads[7,]<-c(ppv_radiomics_birads_cv, left_ppv_radiomics_birads_cv, right_ppv_radiomics_birads_cv)
      cv_radiomics_birads[8,]<-c(npv_radiomics_birads_cv, left_npv_radiomics_birads_cv, right_npv_radiomics_birads_cv)
      cv_radiomics_birads[9,]<-c(mcc_radiomics_birads_cv, left_mcc_radiomics_birads_cv, right_mcc_radiomics_birads_cv)
      
      coef(fit.glmnet_lasso_train_test$finalModel, fit.glmnet_lasso_train_test$bestTune$lambda)
      
      library(coefplot)
      library(RColorBrewer)
      best_lambda_radiomics_birads <- fit.glmnet_lasso_train_test$bestTune$lambda
      
      coeff_radiomics_birads <- coefplot::coefplot(fit.glmnet_lasso_train_test$finalModel, title = typeCausalityAnalysis, xlab = "Coefficients", ylab = "Model Variables", color = brewer.pal(n = 8, name = "Set2")[1], sort='magnitude')
      print(coeff_radiomics_birads)
      
      coefs_model[[typeCausalityAnalysis]] <- coeff_radiomics_birads[["data"]][["Coefficient"]]
      
      cv_results[[typeCausalityAnalysis]] <- cv_radiomics_birads 
      p_class_birads <- predict(fit.glmnet_lasso_train_test, X_test)
      p_prob_birads <- predict(fit.glmnet_lasso_train_test, X_test, type="prob")
      num_l <- as.numeric(as.factor(Y_test))-1
      
      roc_radiomics_birads_test <- pROC::roc(num_l, p_prob_birads[, positive_class])
      
      test_cm <- confusionMatrix(p_class_birads, as.factor(Y_test), positive = positive_class)
      MCC <- (test_cm$table[2,2] * test_cm$table[1,1] - test_cm$table[2,1] * test_cm$table[1,2])/
        (sqrt((test_cm$table[2,2] + test_cm$table[2,1]) * 
                (test_cm$table[2,2] + test_cm$table[1,2]) * 
                (test_cm$table[1,1] + test_cm$table[2,1]) * 
                (test_cm$table[1,1] + test_cm$table[1,2])))
      test_results[[typeCausalityAnalysis]] <- list(test_cm, MCC, roc_radiomics_birads_test)
    }
  }
  train_test_results[[i]] <- list(cv_results, coefs_model, roc_radiomics_birads_cv, test_results)
}
