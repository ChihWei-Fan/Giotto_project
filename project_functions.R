
#Define a function that loop through all the tiles in the folder for feature extraction
extract_image_feats <- function(img_path_i) {
  print(img_path_i)
  img_path <- list_files[[img_path_i]]
  
  # load simple image by reading TIFF tile
  tile <- rast(img_path)
  
  # convert the image to an array, turn into raster object to array
  ff_img <- as.array(tile)
  
  # add an extra dimension (required by the model)
  ff_img <- array_reshape(ff_img, c(1, dim(ff_img)))
  
  # preprocess the image
  ff_img <- imagenet_preprocess_input(ff_img)
  
  # extract features
  ff_features <- resnet_model_shape %>% predict(ff_img)
  
  # return filename and features
  return(c(img_path, ff_features))
}



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

##Functions used in ML_models_for_prediction.R

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



