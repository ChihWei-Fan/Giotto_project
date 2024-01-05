#This works on local computer
# Load necessry libraries
library(keras)
library(ggplot2)
library(data.table)
library(terra)
library(rgdal)
library(tensorflow)
library(reticulate)

#select the conda environment 
reticulate::use_condaenv("/Users/fanzhiwei/Library/r-miniconda-arm64/envs/r-reticulate/bin/python")

# import the H&E image of the sample
#he_image_119 = terra::rast(x = '/Users/fanzhiwei/Desktop/Bioinformatics-intern/Project/improve_res/s119B_tissue_hires_image.png') # rast - A SpatRaster represents a spatially referenced surface divided into three dimensional cells (rows, columns, and layers).
he_image_119 = terra::rast(x = '/Users/fanzhiwei/Desktop/Bioinformatics-intern/Project/improve_res/119B.tif')
plot(he_image_119)

#To check if the he_image_119 object is not empty # Get to know the resolution of the image
dim(he_image_119)


# launch and external plot viewer to work with terra interactive drawing ('draw')
dev.new(noRStudioGD = TRUE)
# create a mask around the region you want to tile
terra::plot(he_image_119)
whole <- draw(x = 'polygon', col = 'blue') #draw() = Create a SpatExtent by drawing it on top of a map (plot)


#Select the region inside you want to remove
inside1 <- draw(x = 'polygon', col = 'black') 


#Select the region inside you want to remove
inside2 <- draw(x = 'polygon', col = 'lightblue') 

# Erase the unwanted region from the whole sample
sample_polygon <- terra::erase(whole,inside1)%>%terra::erase(inside2)
#dev.new(noRStudioGD = TRUE)
#plot(sample_polygon) # check if it erase correctly
he_image_masked = terra::mask(x = he_image_119, mask = sample_polygon )


# draw a new extent #get to know the extent of the polygon (to remove the white space in the ploygon we just draw?)
dev.new(noRStudioGD = TRUE)
terra::plot(he_image_masked)
s119B_extent = draw(x = 'extent', col = 'green')
s119B_extent
#Extend the extent before cropping the image
s119B_extent <- ext(6591, 15831, 6127, 15892) #xmin= 6591, xmax= 15831, ymin=5917, ymax= 15892) 
# ext(he_image_masked) <- c(6591, 15831, 5917, 15892) # Change the extent of the SpatRast directly


# Cropping the image
he_image_masked_crop = terra::crop(he_image_masked, s119B_extent)
dev.new(noRStudioGD = TRUE)
terra::plot(he_image_masked_crop)

#Save the Cropped SpatRaster
#Method3 # This work for me -- M1 CHIP
writeRaster(he_image_masked_crop, "/Users/fanzhiwei/Desktop/keras_s119B_ext_n_cropped.tif", filetype = "GTiff", overwrite = TRUE)

#s119B_ext_n_cropped <-terra::rast("/Users/fanzhiwei/Desktop/Bioinformatics-intern/Project/improve_res/s119B_final/keras_s119B_ext_n_cropped.tif")



