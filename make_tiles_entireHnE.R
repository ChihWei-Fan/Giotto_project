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

# import the H&E image of the sample
he_image_119 = terra::rast(x = '/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/119B.tif')
#he_image_masked_crop <- terra::unwrap(pack_spatRaster)
##plot(he_image_119)

#To check if the he_image_119 object is not empty # Get to know the resolution of the image
dim(he_image_119)


# Cropped the image
#Method 1 
#e <-ext(c(6538, 15778, 6059, 15929))
he_image_masked_crop <- terra::rast("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/keras_s119B_ext_n_cropped.tif", snap="in")
#window(he_image_masked_crop) <-e
terra::plot(he_image_masked_crop)

# create the raster image and specify the number of tiles you want by setting ncols and nrows # set xmin, xmax, ymin, and ymax to the same as the extent of your he_image_masked_crop 
# Checking xmin, xmax, ymin, ymax, and dim first
he_image_masked_crop #
#low_res-- dimensions : 818, 775, 3  (nrow, ncol, nlyr) 
#high_res-- dimensions: 9855, 9225, 3  (nrow, ncol, nlyr) #extent : 6538, 15763, 6059, 15914  (xmin, xmax, ymin, ymax)

#Extend the extent of the image 
#he_image_masked_crop_extend <- extend(he_image_masked_crop, c(6538, 15778, 6059, 15929), fill=NA)

#Decide how many tiles you want to generate from the image # Base on the info from the chuck above
sub <- terra::rast(ncols=88, nrows=93, xmin = 6591, xmax = 15831, ymin = 6127, ymax = 15892)
values(sub) <- 1:ncell(sub) # add values for visualization purpose
plot(sub)


# tile he_image_masked_crop by created raster image (sub)
filename <- paste0('/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_88x93_tiles/s119B_', ".tif")
he_tiles <- makeTiles(x = he_image_masked_crop, y = sub, filename = filename, overwrite= TRUE, extend=TRUE)
#he_tiles



