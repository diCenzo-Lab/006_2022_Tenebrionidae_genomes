# Load libraries
library("ape")
library("vegan")
library("philentropy")
library("ggplot2")
library("ggfortify")
library("data.table")

# Load the data
data <- read.table(file = "parsed_Orthogroups.csv", stringsAsFactors = TRUE, 
                   header = TRUE, sep = ',')
data2 <- t(data[,-1]);

# Prepare distance tree based on the pangenome
r_names <- row.names(data2)
matrix_dist <- distance(data2, method = "jaccard")
row.names(matrix_dist) <- r_names
tree_bionj <- bionj(matrix_dist)
write.tree(tree_bionj, file = "distance.tre")

