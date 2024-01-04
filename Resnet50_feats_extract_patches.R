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
tile_list = list() # for filenames

#Loop through all the tiles in the folder, and seperate the empty tiles from normal tiles
for(img_path_i in 1:length(list_files)) {
  print(img_path_i)
  img_path = list_files[[img_path_i]]
  
  # load simple image #by reading TIFF tile
  tile <- rast(img_path)
  
  # convert the image to an array # turn into raster object to array
  ff_img <- as.array(tile)
  
  # add an extra dimension (required by the model)
  ff_img <- array_reshape(ff_img, c(1, dim(ff_img)))
  # preprocess the image
  ff_img <- imagenet_preprocess_input(ff_img)
  
  # extract features
  ff_features <- resnet_model_shape %>% predict(ff_img)  #vgg16_model_shape #resnet_model_shape
  res_list[[img_path_i]] = c(img_path,ff_features) #associating filename should match to the xempty_tile_ls
}


#Combine all features of the rows (tiles) 
res_dfr = do.call('rbind', res_list)

# Convert results to matrix (image x features)
#tile_names <- data.frame(tile_name = res_dfr[,1])
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