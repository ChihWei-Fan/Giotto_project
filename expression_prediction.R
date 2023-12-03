# Load necessry libraries
if(!require(keras)) {
  install.packages("keras")
}

library(tidyverse)
library(keras)
library(ggplot2)
library(data.table)
library(terra)
library(tensorflow)
library(reticulate)
library(RColorBrewer)
library(Giotto)
library(e1071)
library(caTools)
library(caret)
library(randomForest)
library(glmnet)
library(gbm) 

#select the conda environment 
reticulate::use_condaenv("/projectnb/rd-spat/HOME/ivycwf/.conda/envs/giotto_env_keras/bin/python")

#Load visium data
load(file ="/projectnb/rd-spat/HOME/ivycwf/project_1/sample_119B/visium_119B.RData" )


# Generates deep copy of SpatRaster
full_image <- Giotto::createGiottoLargeImage("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/119B.tif")


#put an large image in the Giotto object
visium_sample_119B@largeImages <- list(image = full_image)


#Get the full size image spatialRaster
fullsize_sr <- full_image@raster_object
#spatPlot2D(visium_sample_119B, show_image = T, largeImage_name = "image")


#Get all the coordinates of spots
cell_coords <-visium_sample_119B@spatial_locs$cell$raw@coordinates %>% as.data.frame()


#$fiducial_diameter_fullres -- 324.4803 (info from scale_json file)
cell_coords <- cell_coords %>%
  mutate(
    xmin = sdimx - 100,
    xmax = sdimx + 100,
    ymin = sdimy - 100,
    ymax = sdimy + 100
  )

# Train the resnet model in patch_run_resnet_model.R
#Get the PCA result of patch tiles # have "res_dfr", "tile_names", "s119B_patch_tiles_pca" variable in it
load(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_afterPCA.RData") 


#Get the patch number with corresponding extent #cells_order
load(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_after_maketiles.RData")

# Check features matrix from resnet50 model
#any(is.nan(image_mat))
#dim(image_mat) #2004 2048


# cells_order Get From make_patch_tiles_4tiles.R file
patches_info <- dplyr::inner_join(cell_coords, cells_order, by = c("xmin", "xmax", "ymin", "ymax"))


#Get the information of expression that links to the tiles
##Get the all the tile names
list_files = res_dfr[,1]

#Concat cell ID,  path number, and tilename in a tibble
#patch_info <- data.frame(cell_ID = cell_coords$cell_ID, patch_num = c(1:501))
tile_name <- data.frame( tile_name = sapply(list_files,basename), patch_number = as.numeric(gsub("s119B_(\\d+)_\\d+\\.tif", "\\1", sapply(list_files,basename))))
patch_tile_info <- dplyr::right_join(patches_info, tile_name, by = c("patch_number" = "patch_number"), multiple = "all")


#Get normalized visium gene expression of all genes
norm_expr_mt <-Giotto::getExpression(visium_sample_119B, values = "normalized",spat_unit = "cell", feat_type = "rna", output = "matrix")%>% 
  as.matrix()


#Get the expression values of specific spatial genes   #"APOE","APOC1","TYROBP","C1QB","C1QA","IFI27","C1QC","CD74","TREM2","CD52","MMP7"
spatial_genes <- t(norm_expr_mt[c("APOE","MMP7"), , drop = FALSE]) %>% 
  as.data.frame()
spatial_genes <-  mutate(spatial_genes, cell_ID = rownames(spatial_genes))
patch_tile_info <- dplyr::inner_join(patch_tile_info, spatial_genes, by = c("cell_ID" = "cell_ID"))


#Concat spatial gene's expression that is associating with those tiles to the df
tiles_df <- data.frame(
  tile_ID = unlist(list_files),
  tile_name = sapply(list_files,basename),
  x_cor = sapply(list_files, function(file) ext(rast(file))[1] + 50),
  y_cor = sapply(list_files, function(file) ext(rast(file))[3] + 50)
)
rownames(tiles_df)<- NULL


#Combine all the tile-related info together
tile_plot_df <- dplyr::inner_join(patch_tile_info, tiles_df, by = c("tile_name" = "tile_name"))

####################################################################################################
# run simple PCA on image feature matrix #image_mat was also scaled in this step 
# This step already been done in patch_run_resnet_model.R

#Get the all pca values for each patch tiles
pca_matrix <- s119B_patch_tiles_pca$x 
pca_matrix <-cbind(tile_names, pca_matrix)

#Plot PCA result from extracting feartures
explained_variance <- data.frame(PC= paste0("PC",1:2004), #used be 968
                                 var_explained= s119B_patch_tiles_pca$sdev^2/sum((s119B_patch_tiles_pca$sdev)^2))

explained_variance$PC <- factor(explained_variance$PC, levels = paste0("PC", 1:2004))

scree_plot_pc30 <-ggplot(explained_variance[1:500,], 
                         aes(x=PC,
                             y=var_explained, 
                             group=1))+
                  geom_point()+
                  geom_line()+
                  labs(title="Scree plot")+
                  theme(axis.text.x = element_text(size = 6))
scree_plot_pc30

# Decide how many PCs you want to filter out from the morphology features
# Using PC1-PC8 #Now use all the PCs
subset_pca_matrix <- pca_matrix[,1:500]

#subset_pca_matrix$tile_name <- gsub("/s119B_", "s119B_", subset_pca_matrix$tile_name)
input_mat <-data.frame()
input_mat <- dplyr::inner_join(subset_pca_matrix, tile_plot_df[,10:12], by = c("tile_name" = "tile_ID"))
#input_mat <- dplyr::inner_join(pca_matrix, tile_plot_df[,10:12], by = c("tile_name" = "tile_ID"))
input_mat[,1]<- sapply(input_mat$tile_name,basename)
###################################################################################################################################################################################################################################################################
#Instead of using PCs, this time we use original features from Resnet50 model







#Save the workspace #Skip this
#save.image(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_correct_tileplotdf.RData")
#get the tiles order first , set it to rownames, and mapping rownames to original expression values 
#load(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_correct_tileplotdf.RData")


set.seed(123)

split = sample.split(input_mat[,1], SplitRatio = 0.8)
training_set = subset(input_mat, split == TRUE)
test_set = subset(input_mat, split == FALSE)


#Create Contol here first

#Put a linear regression model here as a baseline to compare with other ML model, get to know if these models are better




#Use the SVM model 
svm_mod_APOE <- svm(APOE ~ .,
                    data = training_set[,-c(1,202)], #c(1,11)
                    type = 'eps-regression',
                    kernel = 'linear',
                    cost = 10)

# Predicting the Test set results (without giving it real expression value)
#APOE_pred = predict(svm_mod_APOE , newdata = test_set[,-c(1,10,11)])
APOE_pred = predict(svm_mod_APOE , newdata = test_set[,-c(1,201, 202)])

#Generate confusion matrix
cm_APOE = data.frame( APOE = test_set[201], APOE_pred = APOE_pred, tile_name = test_set[,1]) #APOE = test_set[10]


# Combine the original expression value and predict values with tile coordinates
cm_APOE <- merge(tile_plot_df[match(cm_APOE$tile_name, tile_plot_df$tile_name),c(3,9,13,14)], cm_APOE, by = "tile_name")


mse =  mean((test_set$APOE - APOE_pred)^2) 
mae = MAE(test_set$APOE, APOE_pred)
rmse = RMSE(test_set$APOE, APOE_pred)
r2 = R2(test_set$APOE, APOE_pred, form = "traditional")
cat(" MAE:", mae, "\n", "MSE:", mse, "\n", 
    "RMSE:", rmse, "\n", "R-squared:", r2)



#Prepare df for ploting
traing_testdf <- merge(training_set, tile_plot_df[tile_plot_df$tile_name %in% training_set$tile_name, c(9,13,14)], by = "tile_name")
testing_testdf <- merge(test_set, tile_plot_df[tile_plot_df$tile_name %in% test_set$tile_name, c(9,13,14)], by = "tile_name")
train_n_test_testdf <- rbind(traing_testdf, testing_testdf)


#Plot visium gene expression (training set)
ggplot(traing_testdf, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = APOE)) +
  scale_colour_gradient2()


#Polt visium gene expression (test set)
ggplot(cm_APOE, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = APOE)) +
  scale_colour_gradient2() #low = "blue", high = "red"

#or 
ggplot(testing_testdf, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = APOE)) +
  scale_colour_gradient2() 

#Plot visium gene expression (entire dataset)
ggplot(train_n_test_testdf, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = APOE)) +
  scale_colour_gradient2()



#Polt predict expression
ggplot(cm_APOE, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = APOE_pred)) +
  scale_colour_gradient2()


#Plot PCs on training data #To see if PCs showing some spatial patterns
ggplot(traing_testdf, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = PC2)) +
  scale_colour_gradient2()

#Plot PCs on whole data (2004 tiles)
ggplot(train_n_test_testdf, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = PC1)) +
  scale_colour_gradient2()





#Use the SVM model #MMP7
svm_mod_MMP7 <- svm(MMP7 ~ .,
                    data = training_set[,-c(1,10)],
                    type = 'eps-regression',
                    kernel = 'linear',
                    cost = 10)


MMP7_pred = predict(svm_mod_MMP7 , newdata = test_set[,-c(1,10,11)])

#Generate confusion matrix
cm_MMP7 = data.frame( MMP7 = test_set$MMP7, MMP7_pred = MMP7_pred, tile_name = test_set[,1])

# Combine the original expression value and predict values with tile coordinates
cm_MMP7 <- merge(tile_plot_df[match(cm_MMP7$tile_name, tile_plot_df$tile_name), c(3,9,13,14)], cm_MMP7, by = "tile_name")


mse =  mean((test_set$MMP7 - MMP7_pred)^2) 
mae = MAE(test_set$MMP7, MMP7_pred)
rmse = RMSE(test_set$MMP7, MMP7_pred)
r2 = R2(test_set$MMP7, MMP7_pred, form = "traditional")
cat(" MAE:", mae, "\n", "MSE:", mse, "\n", 
    "RMSE:", rmse, "\n", "R-squared:", r2)



#Plot visium gene expression (training set)
ggplot(traing_testdf, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = MMP7)) +
  scale_colour_gradient2()



#Polt visium gene expression (testing_data)
ggplot(cm_MMP7, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = MMP7)) +
  scale_colour_gradient2() #low = "blue", high = "red"

#or 
ggplot(testing_testdf, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = MMP7)) +
  scale_colour_gradient2() 


#Polt predict expression
ggplot(cm_MMP7, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = MMP7_pred)) +
  scale_colour_gradient2()


#save.image(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_correct_workspace.RData")
#load(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_correct_workspace.RData")


########################################################################################################################################################

linearModel <- glm(APOE ~ ., data= training_set[,-c(1,11)], 
                   family=gaussian(link="identity"))

linearPredicts <- predict(linearModel, newdata =  test_set[,-c(1,10,11)])

glm_result <- data.frame( APOE= test_set[,10], APOE_pred = linearPredicts, tile_name = test_set[,1])

summary(linearModel)

mse =  mean((glm_result$APOE - glm_result$APOE_pred)^2) 
mae = MAE(glm_result$APOE, glm_result$APOE_pred)
rmse = RMSE(glm_result$APOE, glm_result$APOE_pred)
r2 = R2(glm_result$APOE, glm_result$APOE_pred, form = "traditional")
cat(" MAE:", mae, "\n", "MSE:", mse, "\n", 
    "RMSE:", rmse, "\n", "R-squared:", r2)

########################################################################################################################################################

#Using a few spots and replicate it for multiple times 
test_sub_rep <- rbind(testing_testdf[1:20,],testing_testdf[1:20,],testing_testdf[1:20,],testing_testdf[1:20,],testing_testdf[1:20,])
test_sub_rep_pred <- testing_testdf[1:20,]

svm_mod_APOE <- svm(APOE ~ .,
                    data = test_sub_rep[,-c(1,8)],
                    type = 'eps-regression',
                    kernel = 'linear',
                    cost = 100)

# Predicting the Test set results (without giving it real expression value)
APOE_pred = predict(svm_mod_APOE , newdata = test_sub_rep_pred[,-c(1,7,8)])

#Generate confusion matrix
cm_APOE = data.frame( APOE = test_sub_rep_pred[7], APOE_pred = APOE_pred, tile_name = test_sub_rep_pred[,1])


# Combine the original expression value and predict values with tile coordinates
cm_APOE <- merge(tile_plot_df[match(cm_APOE$tile_name, tile_plot_df$tile_name),c(3,9,13,14)], cm_APOE, by = "tile_name")


#library(caret)
mse =  mean((test_sub_rep_pred$APOE - APOE_pred)^2) 
mae = MAE(test_sub_rep_pred[, 7], APOE_pred)
rmse = RMSE(test_sub_rep_pred[, 7], APOE_pred)
r2 = R2(test_sub_rep_pred[, 7], APOE_pred, form = "traditional")
cat(" MAE:", mae, "\n", "MSE:", mse, "\n", 
    "RMSE:", rmse, "\n", "R-squared:", r2)
#MAE: 0.4204569 
#MSE: 0.526939 
#RMSE: 0.7259056 
#R-squared: -1.587039

########################################################################################################################################################

#Mergeing the tiles that are red dots in PC1 
tile_names <- lapply(train_n_test_testdf[,1], function(x) paste0('/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/', x))
# Create an empty list to store the raster objects
fftiles <- list()
fftile_rows <- list()

# Load each raster into the list
for (tile_i in 1:length(tile_names)) {
  tile <- terra::rast(tile_names[[tile_i]])
  fftiles[[tile_i]] <- tile
  fftile_rows[[tile_i]] <- nrow(terra::values(fftiles[[tile_i]]))
}

# Stack the rasters in the list on top of each other
fftiles_ls <-terra::sprc(fftiles)
merge_image <- terra::merge(fftiles_ls)

# Plot the merged raster
terra::plot(merge_image)

vis_spots <-Giotto::spatPlot2D(visium_sample_119B, largeImage_name = "image", point_alpha = 0.3)





############################################################################################
#This section does not work
library(grid)

# Load required libraries
library(ggplot2)
library(here)      #Path management
library(prismatic) #Approximation to hexcode colors in console

# Read the raster image using the terra package
raster_image <-raster::stack("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/119B.tif")
#raster_image <- rast("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/119B.tif")

# Convert the raster image to a data frame for ggplot
raster_df <- as.data.frame(raster_image, xy = TRUE)

colnames(raster_df) <- c("x", "y", "Red", "Green", "Blue")

df_subset <- raster_df[1:1000000,]

# Plot the raster image using ggplot
ggplot(data = raster_df, aes(x = x, y = y)) + 
  geom_raster(fill = rgb(r = raster_df$Red, g = raster_df$Green, b = raster_df$Blue, maxColorValue = 255), show.legend = FALSE) +
  scale_fill_identity() + 
  ggtitle("Plot .tif rgb") +
  theme_minimal()

ggsave(Map,                                            #save map
       filename = paste0(here(), "/satellite_img.jpg"), dpi = 200)



#######################################################################################################

#Testing
subset_tile_plot_df <- grep('s119B_[0-9]{1}_[1-4]\\.tif$', list_files, value = TRUE)  

#plot first 9 patches
ggplot(tile_plot_df[tile_plot_df$tile_ID %in% subset_tile_plot_df,], aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = APOE)) +
  scale_colour_gradient2()

#Plot all patches APOE
ggplot(tile_plot_df, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = APOE)) +
  scale_colour_gradient2()

#Plot all patches MMP7
ggplot(tile_plot_df, aes(x=x_cor, y=y_cor)) +
  geom_point(aes(colour = MMP7)) +
  scale_colour_gradient2()


original_expression = merge(cell_coords[,1:3], spatial_genes[1:3], by = "cell_ID")
#all(cell_coords$cell_ID == spatial_genes$cell_ID) #Check if the cell IDs match
ggplot(original_expression, aes(x=sdimx, y=sdimy)) +
  geom_point(aes(colour = APOE), size=2) +
  scale_colour_gradient2()


ggplot(original_expression, aes(x=sdimx, y=sdimy)) +
  geom_point(aes(colour = MMP7), size=2) +
  scale_colour_gradient2()


###################################################################################################################
original_APOE_expression[1:9,]
tile_plot_df[tile_plot_df$tile_ID %in% subset_tile_plot_df,]
original_APOE_expression
tile_plot_df
tiles_df
spatial_genes
cell_coords
patch_tile_info
tile_name # This may be wrong, patch number may not match to the cell ID
##########################################################
#Tiles have the coordinate and the spots also have the coordinate, 
#using this to combine the tiles to the expression value 
#also go back to check the order of the patch/spot that use to generate smaller tiles.

save.image(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_08_06.RData")

