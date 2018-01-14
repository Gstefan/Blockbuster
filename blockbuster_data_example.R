########################################################################################################
#                                                                                                      #
# Example of using Blockbuster on simulated data.                                                      #
#                                                                                                      #
########################################################################################################

# This script loads data (that was simulated with 5 groups) from a csv file and applies the Blockbuster 
# clustering algorithm from Brownlees, Gudmundsson & Lugosi - Community detection in partial correlation 
# networks to it.

########################################################################################################
# Part 1: Preliminaries                                                                                #
########################################################################################################

# External libraries:
library(MASS)     # For multivariate normal.
library(igraph)   # For plotting the graph.

# Number of blocks:
k <- 5

########################################################################################################
# Part 1.1: Functions                                                                                  #
########################################################################################################

# The community detection algorithm. The theoretical results may be found in the paper.
blockbuster <- function(Y, k, f = 0){
  
  # Libraries:
  library(LICORS)   # For k-means++.
  
  # Estimate the SCM:
  S   <- cov(Y)
  
  # Calculate the top k eigenvectors, taking factors into account:
  U   <- eigen(S)$vectors[,(f+1):(k+f)]
  
  # Row normalisation:
  N <- diag(1/(sqrt(rowSums(U^2))) )
  X   <- N %*% U 
  
  # k-means clustering on the top k eigenvectors:
  ind <- kmeanspp(Re(X), k, start ='random', nstart = 50)$cluster
}

########################################################################################################
# Part 2: Reading in the data                                                                          #
########################################################################################################

# Input here the directory of the folder that contains the data file:
dir <- 'C:/Users/Gummi/Documents/Nám/UPF/PhD/Blockbuster/Code/Demo/'

# Reading in the data:
Y   <- as.matrix(read.table(paste(c(dir,'sim_data_example.csv'), collapse = ''), sep=','))

# Removing column and row headers:
colnames(Y) <- NULL
Y <- Y[2:nrow(Y), 2:ncol(Y)]

# Converting the data from strings to numeric:
class(Y) <- "numeric"

# Dimensions of the data:
n   <- ncol(Y)
T   <- nrow(Y) 

########################################################################################################
# Blockbuster clustering                                                                               #
########################################################################################################

# Obtaining the estimated partition indices from the Blockbuster algorithm:
ind.hat       <- blockbuster(Y, k)