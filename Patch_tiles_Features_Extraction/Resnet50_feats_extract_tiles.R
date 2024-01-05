# Load necessry libraries
library(tidyverse)
library(keras)
library(ggplot2)
library(data.table)
library(terra)
library(tensorflow)
library(reticulate)
library(RColorBrewer)
library(Giotto)

#select the conda environment 
reticulate::use_condaenv("/projectnb/rd-spat/HOME/ivycwf/.conda/envs/giotto_env_keras/bin/python")

#Set the resnet_model
resnet_model_shape <- keras::application_resnet50(weights = "imagenet", include_top = FALSE, 
                                                  pooling = "max", input_shape = c(100, 100, 3)) 

#Set the vgg16_model
#vgg16_model_shape <- keras::application_vgg16(weights = "imagenet", include_top = FALSE, 
                                              #pooling = "max", input_shape = c(100, 100, 3))


#Get the filnames and folder
pilot_folder <- "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles"
list_files = list.files(pilot_folder, pattern = 's119B_\\d+_\\d+\\.tif$', full.names = T)


#Createing empty list for storeing filenames and features
res_list = list() # for features
res_list <- lapply(seq_along(list_files), extract_image_feats) #extract_image_feats is a self_defined function in project_functions.R


#Combine all features of the rows (tiles) 
res_dfr = do.call('rbind', res_list)


# Convert results to matrix (image x features)
image_mat <- matrix(as.numeric(res_dfr[,-1]), ncol = 2048) #vgg16 -- ncol = 512
dim(image_mat) #Check how many tiles are not empty


#save(res_dfr, image_mat, file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_resnet50_extracted_feats.RData")
#save(res_dfr, tile_names, image_mat, file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_vgg16_extracted_feats.RData")


# run simple PCA on image feature matrix
#s119B_patch_tiles_pca <- prcomp(image_mat, center = T, scale. = T)
#pdf("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_pca_barplot.pdf") 
#plot(s119B_patch_tiles_pca)
#dev.off()
#save(res_dfr, tile_names, s119B_patch_tiles_pca, file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_afterPCA.RData")