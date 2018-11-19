rm(list = ls())
library(igraph) # visualization of trees and solution. 
library(amap)   # compute several distance metrics https://www.rdocumentation.org/packages/amap/versions/0.8-16/topics/Dist

# set your working directory
setwd("/home/vanhoan310/server/alignProject/Trajan")
# load Trajan.R file
source("Trajan.R")

######### Read the inputs from files: CHOOSE one of the following 3 options. 
# OPTION 1: Given rooted trees and distance matrix. Every row in input tree in the following format: [child node] [parent node]
t1_s = read.table("input/t1.tree", header = FALSE, sep = " ", stringsAsFactors = FALSE)
t2_s = read.table("input/t2.tree", header = FALSE, sep = " ", stringsAsFactors = FALSE)
distance_matrix_s = read.table("input/distance_matrix.csv", header = TRUE, row.names = 1, sep = ",")

# OPTION 2: Given undirected trees with roots specified and distance matrix. 
t1_s = read.table("input/t1.tree", header = FALSE, sep = " ", stringsAsFactors = FALSE)
t1_root_s = toString(read.table("input/t1.root", header = FALSE, sep = " ")[1,1])  # might be given
t2_s = read.table("input/t2.tree", header = FALSE, sep = " ", stringsAsFactors = FALSE)
t2_root_s = toString(read.table("input/t2.root", header = FALSE, sep = " ")[1,1])  # might be given
distance_matrix_s = read.table("input/distance_matrix.csv", header = TRUE, row.names = 1, sep = ",")

# OPTION 3: Given rooted trees (or undirected trees with roots) and single cell data. 
# Distance metric (method) must be speficied to compute the distance matrix. 
t1_s = read.table("input/t1.tree", header = FALSE, sep = " ", stringsAsFactors = FALSE)
t1_data_s = read.table("input/t1_data.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
t2_s = read.table("input/t2.tree", header = FALSE, sep = " ", stringsAsFactors = FALSE)
t2_data_s = read.table("input/t2_data.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
# Declare the method to compute distance matrix, e.g,  euclidean, manhattan, pearson, correlation, spearman, etc, 
# for detail distance metrics please refer to https://www.rdocumentation.org/packages/amap/versions/0.8-16/topics/Dist
method_s <- "correlation" 

######### Declare the penalty scheme
penalty <- "avg"  # use avarege scheme, other: max, dtw

######### Creatre a Trajan object with three options (choose one of the three options, respectively). 
trajan <- Trajan(t1 = t1_s, t2 = t2_s, distance_matrix = distance_matrix_s, penalty = penalty)   # option 1
trajan <- Trajan(t1 = t1_s, t1_root = t1_root_s, t2 = t2_s, t2_root = t2_root_s, distance_matrix = distance_matrix_s, penalty = penalty)   # option 2
trajan <- Trajan(t1 = t1_s, t1_data = t1_data_s, t2 = t2_s, t2_data = t2_data_s, method = method_s, penalty = penalty)   # option 3

######### Export the inputs for binary Trajan.
export(trajan) # with the default name, or with a speficic file names. 
export(trajan, t1_treefName = "t1.tree", t1_mapfName = "t1.map", t2_treefName = "t2.tree", 
       t2_mapfName = "t2.map", distance_matrixfName = "distance_matrix.csv")

######### Plot and visulization. 
plot_first_tree(trajan)  # this will export a first_tree.pdf file plot
plot_second_tree(trajan) # this will export a second_tree.pdf file plot

######### Run binary trajan, e.g, ./trajan t1.tree t1.map t2.tree t2.map optimal_solution.csv 2 e 0 0 0 2

# read the optimal solution.
# optimal_path_s = read.table("optimal_solution.csv", header = FALSE, sep = " ", stringsAsFactors = FALSE)
# optimal_path_s <- optimal_path_s[,1:2]
# # plot the optimal solution.
# plot_solution(trajan, optimal_path = optimal_path_s)


