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
#library(randomForest)
#library(glmnet)
#library(gbm) 
library(scales)
library(ggpubr)
library(dgof)

#Previous steps are done in expression_prediction.R file
#load the imput_mat from RDS file
input_mat <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/input_mat.RDS")
tile_plot_df <- readRDS(file ="/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/tile_plot_df.RDS")
#input_mat <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/input_mat_vgg.RDS") 
#tile_plot_df <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/tile_plot_df_vgg.RDS") 


#Prepare df for checking plots
traing_testdf <- merge(training_set, tile_plot_df[tile_plot_df$tile_name %in% training_set$tile_name, c("tile_name","x_cor","y_cor")], by = "tile_name")
testing_testdf <- merge(test_set, tile_plot_df[tile_plot_df$tile_name %in% test_set$tile_name, c("tile_name","x_cor","y_cor")], by = "tile_name")
train_n_test_testdf <- rbind(traing_testdf, testing_testdf)



# Define a function to calculate performance metrics
calculate_metrics <- function(observed, predicted) {
  mse <- mean((observed - predicted)^2)
  mae <- mean(abs(observed - predicted))
  rmse <- sum((observed - predicted)^2) / sum((observed - mean(observed))^2)
  r2 <- R2(observed, predicted, form = "corr")
  pearson <- cor(predicted, observed, method = "pearson")
  
  cat("MAE:", mae, "\n", "MSE:", mse, "\n", 
      "RMSE:", rmse, "\n", "R-squared:", r2, "\n", "Pearson:", pearson, "\n")
}

# Define a function to create expression mapping and relationship plots 
plot_expression <- function(data, gene_name, plot_type) {
  if (plot_type == "observed") {
    col = gene_name
    title = paste("Observed Expression of", gene_name)
  } else if (plot_type == "predicted") {
    col = paste(gene_name, "_pred", sep = "")
    title = paste("Predicted Expression of", gene_name)
  } else if (plot_type == "relation") {
    observed_col = gene_name
    predicted_col = paste(gene_name, "_pred", sep = "")
    return(ggscatter(data, x = observed_col, y = predicted_col, 
                     add = "reg.line", conf.int = TRUE, 
                     cor.coef = TRUE, cor.method = "pearson",
                     add.params = list(color = "blue", fill = "lightgray"),
                     xlab = paste("Observed", observed_col, "Expression"), 
                     ylab = paste("Predicted", predicted_col, "Expression"),
                     size = 1.5, cor.coef.size = 5))
  } else {
    stop("Invalid plot type. Choose 'observed', 'predicted', or 'relation'.")
  }
  
  ggplot(data, aes(x = x_cor, y = y_cor)) +
    geom_point(aes(colour = .data[[col]])) +
    scale_colour_gradient2(low = "blue", high = "red") +
    ggtitle(title)
}

#Check the expression of entire dataset of a single gene
#Example:
#plot_expression(train_n_test_testdf, "SPARC", "observed")
#plot_expression(traing_testdf, "SPARC", "observed")


#Define a train and predict model function
model_train_predict <- function(gene_name, train_set, methods, test_set){
  control <- trainControl(method="cv",                            
                          number = 10,                            
                          summaryFunction = defaultSummary,                         
                          savePredictions = 'all')
  model <- caret::train(
    formula(paste(gene_name, "~ .")),                      
    data = train_set,                    
    method = methods,                      
    metric = "RMSE",                      
    trControl = control) 
  
  # Make predictions on the test set
  gene_preds <- predict(model, newdata = test_set, type = "raw")
  return(list(model = model, gene_preds = gene_preds))
}


#Model 1: #LASSO or elasticnet (glmnet) not using lm because linear regression model need more observations than variables(features) so the PC could be more reliable. -- use original features
#Put LASSO as a baseline to compare with other ML model, get to know if these models are better
#Model 1 -- LASSO (another method)
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


#Model comparisom
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


