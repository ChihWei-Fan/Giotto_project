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
library(scales)
library(ggpubr)

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
load(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_afterPCA.RData") 
#Have "res_dfr", "tile_names", "s119B_patch_tiles_pca" variable in it
#res_dfr -- contain all the features from each spot-covered tiles before performing PCA


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
patch_tile_info <- data.frame()
patch_tile_info <- dplyr::right_join(patches_info, tile_name, by = c("patch_number" = "patch_number"), multiple = "all")


#Get normalized visium gene expression of all genes
scaled_expr_mt <-Giotto::getExpression(visium_sample_119B, values = "scaled",spat_unit = "cell", feat_type = "rna", output = "matrix")%>% 
  as.matrix()


#Get the expression values of specific spatial genes (Breast cancer related marker genes)
spatial_genes <-data.frame()
spatial_genes <- t(scaled_expr_mt[c("SPARC","COL1A1","LUM","SFRP2","COL3A1","SULF1","COL1A2","VCAN","IGFBP7","COL18A1","THY1"), , drop = FALSE]) %>% 
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
tile_plot_df <- data.frame()
tile_plot_df <- dplyr::inner_join(patch_tile_info, tiles_df, by = c("tile_name" = "tile_name"))


# Use original features from Resnet50 model would be reliable than using PCs
# Convert results to matrix (image x features)
tile_names <- data.frame(tile_name = res_dfr[,1])
image_mat <- matrix(as.numeric(res_dfr[,-1]), ncol = 2048) %>% as.data.frame()
dim(image_mat) #Check how many tiles are not empty

# modify the column names avoid using number as column names
colnames(image_mat) <- paste0("f", seq_along(image_mat))

#scale the features
#scal_image_mat <- scale(image_mat, center = T, scale = T) #scaled by col

#Combine the features with the corresponding tile name 
features_matrix <-cbind(tile_names, image_mat)
#scal_features_matrix <-cbind(tile_names, scal_image_mat)

#Create input matrix #Do not need to do this every time.
input_mat <- data.frame()
input_mat <- dplyr::inner_join(features_matrix, tile_plot_df[,10:21], by = c("tile_name" = "tile_ID")) #select target gene and tile_ID columns
input_mat[,1] <- sapply(input_mat$tile_name, basename) #input_mat include tile_name, original 2048 features, and target genes' expression

#scal_input_mat <- data.frame()
#scal_input_mat <- dplyr::inner_join(scal_features_matrix, tile_plot_df[,10:13], by = c("tile_name" = "tile_ID"))
#scal_input_mat[,1] <- sapply(scal_input_mat$tile_name, basename) #input_mat include tile_name, scaled 2048 features, and target genes' expression


#Don't need to do this every time (unless you made some changes in input_mat)
#saveRDS(input_mat, file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/input_mat.RDS") 
#saveRDS(scal_input_mat, file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/scal_input_mat.RDS") 


#load the imput_mat from RDS file
input_mat <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/input_mat.RDS")
#scal_input_mat <- readRDS(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/scal_input_mat.RDS")


#split the dataset 
set.seed(123)
split = sample.split(input_mat[,1], SplitRatio = 0.65)
training_set = subset(input_mat, split == TRUE)
test_set = subset(input_mat, split == FALSE) 

#split the scaled dataset
#set.seed(123)
#split = sample.split(scal_input_mat[,1], SplitRatio = 0.8)
#training_set = subset(scal_input_mat, split == TRUE)
#test_set = subset(scal_input_mat, split == FALSE)

#Model 1: #LASSO or elasticnet (glmnet) not using lm because linear regression model need more observations than variables(features) so the PC could be more reliable. -- use original features
#Put LASSO as a baseline to compare with other ML model, get to know if these models are better








#########################################################################################################################################################################################
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








