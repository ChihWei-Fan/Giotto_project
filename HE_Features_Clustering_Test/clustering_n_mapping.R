# Load necessry libraries
if(!require(keras)) {
  install.packages("keras")
}

library(keras)
library(ggplot2)
library(data.table)
library(terra)
library(tensorflow)
library(reticulate)
library(RColorBrewer)

#Modify the time for Save.image() to save the all things in workspace
options(eval.sec = 600)

#select the conda environment 
reticulate::use_condaenv("/projectnb/rd-spat/HOME/ivycwf/.conda/envs/giotto_env_keras/bin/python")

#Load the workspace
load(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_89x103_tiles/entireHnE_resnet_feats_afterPCA.RData")

# Check if I remove and fill all the NaN values 
any(is.nan(image_mat))
dim(image_mat) # 4854 2048

# run simple PCA on image feature matrix
#s119B_pca <- prcomp(image_mat, center = T, scale. = T)


#Scree Plot (y-axis explained_var) #Notice the PC values need to change depend on the tiles numbers that are not empty 
### When the observation (row) is more than features, prcomp() will only take 2048 PCs
explained_variance <- data.frame(PC= paste0("PC",1:2048),
                                 var_explained= (s119B_pca$sdev^2/sum((s119B_pca$sdev)^2))*100
                                 )

explained_variance$PC <- factor(explained_variance$PC, levels = paste0("PC", 1:2048))

scree_plot_pc15 <-ggplot(explained_variance[1:15,], 
       aes(x=PC,
           y=var_explained, 
           group=1))+
  geom_point()+
  geom_line()+
  ylab("var_explained (%)")
  labs(title="Scree plot")+
  theme(axis.text.x = element_text(size = 8))

pdf("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_89x103_tiles/entireHnE_scree_plot_pc15.pdf")
scree_plot_pc15
dev.off()

#Get the top N pca results
pca_matrix <- s119B_pca$x 

# Decide how many PCs you want to be as input to do kmeans clustering
subset_pca_matrix <- pca_matrix[,1:10]            #centers =5, nstart = 10

# Decide how many clusters you want #
kmeans_clust8 = stats::kmeans(subset_pca_matrix, centers = 8, iter.max = 500, nstart = 100)

#Plot pc1 vs pc2 matrix color by diff clusters
#dev.new(noRStudioGD = TRUE) 
# Decide how many clusters you want

# Create a data frame with PCA results and cluster labels
df_PCA_clust_result <- data.frame(PC1 = pca_matrix[, 1], PC2 = pca_matrix[, 2], cluster = as.factor(kmeans_clust8$cluster))


# Get the number of clusters
#num_clusters <- unique(df$cluster) %>% length()

# Define a combined palette with colors from different palettes
combined_palette <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"))


# Define a ggplot with x and y aesthetics and cluster as color
pdf("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_88x93_tiles/entireHnE_pca10_cen8_it500_star100.pdf")
ggplot(df_PCA_clust_result, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) + scale_color_manual(values = combined_palette[1:num_clusters])
dev.off()


# Combine filename, associating coordinates, pca and cluster information together in a dataframe

#First get filenames(ID), PC1, PC2, cluster info together
pca_res = data.table::data.table(ID = unlist(lapply(xempty_ls, basename)),
                                 PC1 = pca_matrix[,1], PC2 = pca_matrix[,2], kmeans = kmeans_clust8$cluster)


# Extract x and y coordinates
# Create an empty dataframe to store the results
tile_coordin <- data.frame(ID = character(), x = numeric(), y = numeric())

# Iterate over the filenames and their coordinate
for (file_i in 1:length(pca_res$ID)) {
  filename <- pca_res$ID[[file_i]]
  
  # Extract x and y coordinates # the value follow after %% or %/% number is decides by the ncols
  numeric_part <- as.numeric(gsub("^s119B_(\\d+)\\.tif$", "\\1", filename))
  x <- ifelse((numeric_part %% 89) == 0, 89, (numeric_part %% 89))
  y <- (numeric_part - 1) %/% 89 + 1
  
  # Add the coordinates to the dataframe
  tile_coordin <- rbind(tile_coordin, data.frame(ID = filename, x = x, y = y))
}

# Print the resulting dataframe
print(tile_coordin)


# combine based on file name to form a table with all info we need
comb_res = data.table::merge.data.table(x = tile_coordin, y = pca_res, by = "ID")

# visualize the result 
#dev.new(noRStudioGD = TRUE) 
pl <- ggplot(comb_res, aes(x = as.factor(x), y = as.factor(-y), color = as.factor(kmeans))) +
  geom_point(size = 2) + scale_color_brewer(palette = "Set2") + theme_dark() +
  theme(axis.text.x = element_text(size = 6),   # Change the font size of x-axis tick labels
        axis.text.y = element_text(size = 6))   
  
pl
pdf("/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_89x103_tiles/entireHnE_kmean_pca10_cen8_it500_star100_map.pdf", width=11.2, height=10)
pl
dev.off()


save.image(file = "/projectnb/rd-spat/HOME/ivycwf/project_1/resolution/s119B_89x103_tiles/entireHnE_kmean8_pca10_it500_star100.RData")





