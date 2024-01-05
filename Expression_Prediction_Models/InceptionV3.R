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
library(abind)
#library(Giotto) #Giotto is not available for R 4.3
tile_plot_df <- readRDS(file ="/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles/tile_plot_df.RDS")


#Set the InceptionV3 model
inception_model_base <- application_inception_v3(weights = 'imagenet', include_top = FALSE, 
                                                  pooling = "max", input_shape = c(100, 100, 3))

# Add custom layers for your task
# incept_model <- inception_model_base %>%
#                   layer_global_average_pooling_2d() %>%
#                   layer_dense(units = 1024, activation = 'relu') %>%
#                   layer_dense(units = 1, activation = 'sigmoid')

# Compile the model
inception_model_base %>% compile(
  loss = 'mean_squared_error', #'binary_crossentropy'
  optimizer = optimizer_adam(learning_rate = 0.0001),
  metrics = c('mean_absolute_error') #'accuracy'
)

#Load and Preprocess Images
#Get the filnames and folder
pilot_folder <- "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/patch_tiles_4tiles"
list_files = list.files(pilot_folder, pattern = 's119B_\\d+_\\d+\\.tif$', full.names = T)


# image_list <- lapply(list_files[1:2], function(img_path) {
#   tile <- terra::rast(img_path)
#   ff_img <- as.array(tile)
#   ff_img <- array_reshape(ff_img, c(1, dim(ff_img)))
#   img_array <- inception_v3_preprocess_input(ff_img)
#   return(img_array)
# })

image_list <- setNames(lapply(seq_along(list_files), function(i) {
                        cat(i,"\n")
                        img_path <- list_files[i]
                        tile <- terra::rast(img_path)
                        ff_img <- as.array(tile)
                        ff_img <- array_reshape(ff_img, c(1, dim(ff_img)))
                        img_array <- inception_v3_preprocess_input(ff_img)
                        return(img_array)
                        })
                       , 
                       # Extract partial file name link to image array
                       sapply(list_files, function(x) sub("^.*/(s119B_\\d+_\\d+\\.tif)$", "\\1", x))
                      )

#Access the array corresponding to the first image path
#first_image_array <- image_list[[list_files[1]]]

# Combine all preprocessed images into one array
whole_data <- abind(image_list, along = 1) #array is stacked base on first dimension

#Get the image path 
#first_img_path = dimnames(train_data)[[1]][1]
lists_img_path = dimnames(whole_data)[[1]]


# Placeholder for training labels
# train_labels <- ...

#split the dataset 
set.seed(123)
split = sample.split(lists_img_path, SplitRatio = 0.70)
train_data <- whole_data[split, , , ]
test_data <- whole_data[!split, , , ]

#Get the expression value of 
train_labels <- tile_plot_df[tile_plot_df$tile_name %in% dimnames(train_data)[[1]], c("tile_name","SPARC","COL1A1","LUM","SFRP2","COL3A1","SULF1","COL1A2","VCAN","IGFBP7","COL18A1","THY1")] 


# Train the model
history <- inception_model_base %>% fit(
  x = train_data,
  y = train_labels,
  epochs = 25,
  batch_size = 32,
  #validation_data = list(validation_data, validation_labels) 
)


# Add custom layers for your specific task
model <- inception_model %>%
  keras::layer_global_average_pooling_2d() %>%
  keras::layer_dense(units = 1024, activation = 'relu') %>%
  keras::layer_dense(units = 1, activation = 'sigmoid')  # Adjust the units and activation based on your output

# Compile the model
model %>% compile(
  loss = 'binary_crossentropy',  # Adjust loss function based on your task
  optimizer = optimizer_adam(lr = 0.0001),
  metrics = 'accuracy'
)

# Train the model
# Note: You need to provide your training and validation datasets.
history <- model %>% fit(
  x = train_data,  # Replace with your training data
  y = train_labels,  # Replace with your training labels
  epochs = 25,
  batch_size = 32,
  validation_data = list(validation_data, validation_labels)  # Replace with your validation data and labels
)
