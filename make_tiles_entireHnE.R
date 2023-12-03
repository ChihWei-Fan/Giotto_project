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
library(tensorflow)
library(reticulate)

#Modify the time for Save.image() to save the all things in workspace
options(eval.sec = 600)

#select the conda environment 
reticulate::use_condaenv("/projectnb/rd-spat/HOME/ivycwf/.conda/envs/giotto_env_keras/bin/python")

# import the H&E image of the sample
he_image_119 = terra::rast(x = '/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/119B.tif')


#To check if the he_image_119 object is not empty # Get to know the resolution of the image
dim(he_image_119)


# Import the cropped image
he_image_masked_crop <- terra::rast("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/keras_s119B_ext_n_cropped_for100.tif", snap="in")
terra::plot(he_image_masked_crop)



# Checking xmin, xmax, ymin, ymax, and dim first
he_image_masked_crop 



#Decide how many tiles you want to generate from the image # Base on the info from the chuck above
# specify the number of tiles you want by setting ncols and nrows
# set xmin, xmax, ymin, and ymax to the same as the extent of your he_image_masked_crop 
sub <- terra::rast(ncols=89, nrows=103, xmin = 6800, xmax = 15700, ymin = 5700, ymax = 16000)
values(sub) <- 1:ncell(sub) # add values for visualization purpose
plot(sub)



# tile he_image_masked_crop by created raster image (sub)
filename <- paste0('/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_89x103_tiles/s119B_', ".tif")
he_tiles <- makeTiles(x = he_image_masked_crop, y = sub, filename = filename, overwrite= TRUE, extend=TRUE)
#he_tiles



