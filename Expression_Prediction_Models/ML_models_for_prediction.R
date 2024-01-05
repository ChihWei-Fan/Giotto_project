# Load necessry libraries
library(tidyverse)
library(keras)
library(data.table)
library(terra)
library(tensorflow)
library(reticulate)
library(RColorBrewer)
library(Giotto)
library(e1071)
library(caTools)
library(caret)
library(scales)
library(ggpubr)
library(dgof)
source(file = "")

#Previous steps are done in MLmodel_inout_prep.R file
#load the imput_mat from RDS file
input_mat <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/input_mat.RDS")
tile_plot_df <- readRDS(file ="/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/tile_plot_df.RDS")
#input_mat <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/input_mat_vgg.RDS") 
#tile_plot_df <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/tile_plot_df_vgg.RDS") 

#split the dataset 
set.seed(123)
split = sample.split(input_mat[,1], SplitRatio = 0.65)
training_set = subset(input_mat, split == TRUE)
test_set = subset(input_mat, split == FALSE) 


#Prepare df for checking plots
traing_testdf <- merge(training_set, tile_plot_df[tile_plot_df$tile_name %in% training_set$tile_name, c("tile_name","x_cor","y_cor")], by = "tile_name")
testing_testdf <- merge(test_set, tile_plot_df[tile_plot_df$tile_name %in% test_set$tile_name, c("tile_name","x_cor","y_cor")], by = "tile_name")
train_n_test_testdf <- rbind(traing_testdf, testing_testdf)


#Model 1: #LASSO or elasticnet (glmnet) not using lm because linear regression model need more observations than variables(features) so the PC could be more reliable. -- use original features
#Put LASSO as a baseline to compare with other ML model, get to know if these models are better
model_n_preds <- model_train_predict("SPARC", training_set[,c(2:2049,2050)], "glmnet", test_set[,c(2:2049)] )
lasso_model <- model_n_preds$model
gene_preds <- model_n_preds$gene_preds
#Generate confusion matrix  
cm_gene_lso = data.frame( SPARC = test_set["SPARC"], 
                          SPARC_pred = gene_preds, 
                          tile_name = test_set[,1])

#Combine the original expression value and predict values with tile coordinates
cm_SPARC_lso <- merge(tile_plot_df[match(cm_gene_lso$tile_name,tile_plot_df$tile_name), c("cell_ID","tile_name","x_cor","y_cor")], 
                      cm_gene_lso, by = "tile_name")

#Calculate model evaluation matrics
calculate_metrics(cm_SPARC_lso$SPARC, cm_SPARC_lso$SPARC_pred)

#Plot predicted expression
plot_expression(cm_SPARC_lso, "SPARC", "predicted")
plot_expression(cm_SPARC_lso, "SPARC", "relation")

head(cm_SPARC_lso)

#Preleminary test to check the test assumptions
cor.test(cm_SPARC_lso$SPARC, cm_SPARC_lso$SPARC_pred, method=c("pearson"))
cor.test(cm_SPARC_lso$SPARC, cm_SPARC_lso$SPARC_pred, method=c("spearman"))

#Test the SPARC and predicted SPARC is normally distributed
# QQ-plot normality test for SPARC
ggqqplot(cm_SPARC_lso$SPARC, ylab = "SPARC")+
  ggtitle("Quantile-Quantile Plot of Observed SPARC Expression")

# QQ-plot normality test for predicted SPARC
ggqqplot(cm_SPARC_lso$SPARC_pred, ylab = "SPARC_pred")+
  ggtitle("Quantile-Quantile Plot of Predicted SPARC Expression")
#Kolmogorov-Smirnov test
ks.test(cm_SPARC_lso$SPARC, "pnorm")
ks.test(cm_SPARC_lso$SPARC_pred, "pnorm")



#Model 2: Random Forest model --COL1A1
model_n_preds <- model_train_predict("COL1A1", training_set[,c(2:2049,2051)], "rf", test_set[,c(2:2049)] )
rf_model <- model_n_preds$model
gene_preds <- model_n_preds$gene_preds
#Generate confusion matrix  
cm_gene_rf = data.frame( COL1A1 = test_set["COL1A1"], 
                         COL1A1_pred = gene_preds, 
                         tile_name = test_set[,1])

#Combine the original expression value and predict values with tile coordinates
cm_COL1A1_rf <- merge(tile_plot_df[match(cm_gene_rf$tile_name,tile_plot_df$tile_name), c("cell_ID","tile_name","x_cor","y_cor")], 
                      cm_gene_rf, by = "tile_name")

#Calculate model evaluation matrics
calculate_metrics(cm_COL1A1_rf$COL1A1, cm_COL1A1_rf$COL1A1_pred)

#Plot predicted expression
plot_expression(cm_COL1A1_rf, "COL1A1", "COL1A1_pred")

#Preleminary test to check the test assumptions
cor.test(cm_COL1A1_rf$COL1A1, cm_COL1A1_rf$COL1A1_pred, method=c("pearson"))
cor.test(cm_COL1A1_rf$COL1A1, cm_COL1A1_rf$COL1A1_pred, method=c("spearman"))
# normality check for COL1A1
ggqqplot(cm_COL1A1_rf$COL1A1, ylab = "COL1A1")
ks.test(cm_COL1A1_rf$COL1A1, "pnorm")
# normality check for predicted SPARC
ggqqplot(cm_COL1A1_rf$COL1A1_pred, ylab = "COL1A1_pred")
ks.test(cm_COL1A1_rf$COL1A1_pred, "pnorm")

#skewness of COL1A1 
skewness(cm_COL1A1_rf$COL1A1,na.rm=T) # is moderately skewed #not necessarily need to transformed
skewness(cm_COL1A1_rf$COL1A1_pred,na.rm=T)


#Model 3 GBM -LUM
model_n_preds <- model_train_predict("LUM", training_set[,c(2:2049,2052)], "gbm", test_set[,c(2:2049)] )
gbm_model <- model_n_preds$model
gene_preds <- model_n_preds$gene_preds
#Generate confusion matrix  
cm_gene_gbm = data.frame( LUM = test_set["LUM"], 
                          LUM_pred = gene_preds, 
                          tile_name = test_set[,1])

#Combine the original expression value and predict values with tile coordinates
cm_LUM_gbm <- merge(tile_plot_df[match(cm_gene_gbm$tile_name,tile_plot_df$tile_name), c("cell_ID","tile_name","x_cor","y_cor")], 
                    cm_gene_gbm, by = "tile_name")

#Calculate model evaluation matrics
calculate_metrics(cm_LUM_gbm$LUM, cm_LUM_gbm$LUM_pred)

#Plot predicted expression
plot_expression(cm_LUM_gbm, "LUM", "LUM_pred")


#Preleminary test to check the test assumptions
cor.test(cm_LUM_gbm$LUM, cm_LUM_gbm$LUM_pred, method=c("pearson"))
cor.test(cm_LUM_gbm$LUM, cm_LUM_gbm$LUM_pred, method=c("spearman"))
# normality check for LUM
ggqqplot(cm_LUM_gbm$LUM, ylab = "LUM")
ks.test(cm_LUM_gbm$LUM, "pnorm")
# normality check for predicted LUM
ggqqplot(cm_LUM_gbm$LUM_pred, ylab = "LUM_pred")
ks.test(cm_LUM_gbm$LUM_pred, "pnorm")

#skewness of COL1A1 
skewness(cm_LUM_gbm$LUM,na.rm=T) # is moderately skewed #not necessarily need to transformed
skewness(cm_LUM_gbm$LUM_pred,na.rm=T)

#Model4:  SVM --LUM
model_n_preds <- model_train_predict("LUM", training_set[,c(2:2049,2052)], "svmPoly", test_set[,c(2:2049)] )
svm_model <- model_n_preds$model
gene_preds <- model_n_preds$gene_preds
#Generate confusion matrix  
cm_gene_svm = data.frame( LUM = test_set["LUM"], 
                          LUM_pred = gene_preds, 
                          tile_name = test_set[,1])

#Combine the original expression value and predict values with tile coordinates
cm_LUM_svm <- merge(tile_plot_df[match(cm_gene_svm$tile_name,tile_plot_df$tile_name), c("cell_ID","tile_name","x_cor","y_cor")], 
                    cm_gene_svm, by = "tile_name")

#Calculate model evaluation matrics
calculate_metrics(cm_LUM_svm$LUM, cm_LUM_svm$LUM_pred)

#Plot predicted expression
plot_expression(cm_LUM_svm, "LUM", "LUM_pred")


#Model comparison
model_list <- list(LASSO = lasso_model, RF = rf_model, GBM= gbm_model, SVM = svm_model) 
all_results <- resamples(model_list)
summary(all_results)
bwplot(all_results, metric = "Rsquared")







###Check the correctness of the plotting result of LSSO model (training set + testing set)
combined_data <- rbind(data.frame(x_cor = cm_SPARC_lso$x_cor, y_cor = cm_SPARC_lso$y_cor, SPARC = cm_SPARC_lso$SPARC, dataset = "Lasso"),
                       data.frame(x_cor = traing_testdf$x_cor, y_cor = traing_testdf$y_cor, SPARC = traing_testdf$SPARC, dataset = "Training Set"))

# Plot combined data
ggplot(combined_data, aes(x = x_cor, y = y_cor, color = SPARC, shape = dataset)) +
  geom_point() +
  scale_colour_gradient2(low = "blue", high = "red" )


