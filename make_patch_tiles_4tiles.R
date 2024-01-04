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

#Load visium data
load(file ="/projectnb/rd-spat/HOME/ivycwf/project_1/sample_119B/visium_119B.RData" )

# Generates deep copy of SpatRaster
full_image <- Giotto::createGiottoLargeImage("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/119B.tif")

#put an large image in the Giotto object
visium_sample_119B@largeImages <- list(image = full_image)
#spatPlot2D(visium_sample_119B, show_image = T, largeImage_name = "image")

#know the information of the images store in the Giotto object 
#showGiottoImageNames(visium_sample_119B)

#Get all the spatial location
sl= getSpatialLocations(visium_sample_119B)

#Get the scale factors that is used here
scale_json_path = '/projectnb/rd-spat/HOME/ivycwf/GSE210616/119B/GSM6433605_119B_scalefactors_json.json'

#Read in the scale factor file
scalefactors = jsonlite::read_json(scale_json_path)

#Get the full size image spatialRaster
fullsize_sr <- full_image@raster_object

#Get all the coordinates of spots
cell_coords <-visium_sample_119B@spatial_locs$cell$raw@coordinates %>% as_tibble()

#Get the Xmin Xmax Ymin Ymax for each cell (visium spot)
#(info from scale_json file) #$spot_diameter_fullres -- 200.8687 # the dimension of each tile = 100 *100
cell_coords <- cell_coords %>%
  mutate(
    xmin = sdimx - 100,
    xmax = sdimx + 100,
    ymin = sdimy - 100,
    ymax = sdimy + 100
  )


# Create a list of SpatExtent objects using the extent regions and include the cell_ID
patches_size <- lapply(1:nrow(cell_coords), function(i) {
  extent <- terra::ext(cell_coords$xmin[i], cell_coords$xmax[i], cell_coords$ymin[i], cell_coords$ymax[i])
  attr(extent, "cell_ID") <- cell_coords$cell_ID[i]
  extent
})
# Get the cell_ID for each SpatExtent object
patches_size_df <- data.frame(SpatExtent = character(),
                        cell_ID = character(),
                        stringsAsFactors = FALSE)

# Loop through each element in the patches_size list and extract the "cell_ID" attribute
for (i in seq_along(patches_size)) {
  extent <- patches_size[[i]]
  cell_ID <- attr(extent, "cell_ID")
  patches_size_df <- rbind(patches_size_df, data.frame(SpatExtent = as.character(extent),
                                           cell_ID = cell_ID,
                                           stringsAsFactors = FALSE))
}

# Cropping the image for all the patches
patch_crops <- lapply(patches_size, function(patch_i) {
  terra::crop(fullsize_sr, patch_i)
})


cell_IDs = list()
#Making tiles (cutting one patch(visium spot) into 4 tiles)
for(patch_i in 1:length(patch_crops)){
  #print(patch_i)
  cell_IDs[[patch_i]] = terra::ext(patch_crops[[patch_i]])
  sub <- terra::rast(ncols=2, nrows=2, terra::ext(patch_crops[[patch_i]]))
  values(sub) <- 1:ncell(sub)
  
  filename <- paste0('/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_',patch_i, "_", ".tif")
  #he_tiles <- makeTiles(x = patch_crops[[patch_i]], y = sub, filename = filename, overwrite= TRUE, extend=TRUE)
}


cell_order <- lapply(1:length(cell_IDs), function(patch_i) {
  extents <- terra::ext(cell_IDs[[patch_i]])
  c(patch_number = patch_i, xmin = extents$xmin, xmax = extents$xmax, ymin = extents$ymin, ymax = extents$ymax)
})

cells_order <- do.call(rbind, cell_order)%>% as.data.frame()
colnames(cells_order) <- c("patch_number", "xmin", "xmax", "ymin", "ymax")

#save(cells_order, file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_patch_after_maketiles.RData")
#saveRDS(cell_coords, file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/s119B_cell_coords.RDS")


