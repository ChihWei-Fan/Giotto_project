# Load necessry libraries
if(!require(keras)) {
  install.packages("keras")
}

if(!require(raster)) {
  install.packages("raster")
}

library(keras)
library(ggplot2)
library(data.table)
library(terra)
library(rgdal)
library(tensorflow)
library(reticulate)

#Modify the time for Save.image() to save the all things in workspace
options(eval.sec = 600)

#select the conda environment 
reticulate::use_condaenv("/projectnb/rd-spat/HOME/ivycwf/.conda/envs/giotto_env_keras/bin/python")


# Load resnet50 model for featurization and specify input_shape # the input_shape is the shape of the individual tile -- the resolution (ex: dimension of ff1)
resnet_model_shape <- keras::application_resnet50(weights = "imagenet", include_top = FALSE, pooling = "max", input_shape = c(105, 105, 3)) 

pilot_folder <- "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_88x93_tiles"
list_files = list.files(pilot_folder, pattern = 's119B_\\d+\\.tif$', full.names = T)


#Createing empty list for storeing filenames 
res_list = list()
xempty_tile_ls = list()


#Loop through all the tiles in the folder, and seperate the empty tiles from normal tiles
for(img_path_i in 1:length(list_files)) {
  print(img_path_i)
  img_path = list_files[[img_path_i]]
  
  # load simple image #by reading TIFF tile
  tile <- rast(img_path)
  minmax_RGB <- terra::minmax(tile)
  
  #Filter out those empty tiles
  if (all(!is.nan(minmax_RGB))){
    #print(list_files[[img_path_i]])
    xempty_tile_ls[[img_path_i]] = list_files[[img_path_i]]
    
    #Finding those tiles are portion being mask that cause NaN
    #And fill those NaN with mean of RBG values from background tiles (select seperately)
    if (is.na(any(terra::values(tile) > 255) || any(terra::values(tile) < 0))) {
      tile_values <- terra::values(tile)
      NaN_rows <- which(apply(tile_values, 1, function(row) all(is.nan(row))))
      tile_values[NaN_rows, 1] <- 179
      tile_values[NaN_rows, 2] <- 184
      tile_values[NaN_rows, 3] <- 168
      terra::values(tile) <- tile_values
    } 
    
    # convert the image to an array # turn into raster object to array
    ff_img <- as.array(tile)
    
    # add an extra dimension (required by the model)
    ff_img <- array_reshape(ff_img, c(1, dim(ff_img)))
    # preprocess the image
    ff_img <- imagenet_preprocess_input(ff_img)
    
    # extract features
    ff_features <- resnet_model_shape %>% predict(ff_img)
    res_list[[img_path_i]] = ff_features #associating filename should match to the xempty_tile_ls
  }
}

# Remove all the NULL in the xempty_tile_ls
xempty_ls <- Filter(Negate(is.null), xempty_tile_ls)
#Combine all features of the rows (tiles) 
res_dfr = do.call('rbind', res_list) # at this point the tile will match to the xempty_ls


# Convert results to matrix (image x features)
image_mat <- as.matrix(res_dfr)
dim(image_mat) #Check how many tiles are not empty

# run simple PCA on image feature matrix
s119B_pca <- prcomp(image_mat, center = T, scale. = T)
pdf("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_88x93_tiles/s119B_pca_barplot.pdf") 
plot(s119B_pca)
dev.off()

#save.image(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_88x93_tiles/s119B_afterPCA.RData")




