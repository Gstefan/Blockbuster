########################################################################################################
#                                                                                                      #
# Example of using Blockbuster on simulated data.                                                      #
#                                                                                                      #
########################################################################################################

# This script simulates a stochastic block model and uses it to create a concentration matrix as in the
# paper (Brownlees, Gudmundsson & Lugosi - Community detection in partial correlation networks). It 
# then draws data from a multivariate normal and applies the Blockbuster clustering algorithm. The 
# results are plotted against the actual clustering for comparison.

########################################################################################################
# Part 1: Preliminaries                                                                                #
########################################################################################################

# External libraries:
library(MASS)     # For multivariate normal.
library(igraph)   # For plotting the graph.

# General parameters:
n          <- 100               # Number of vertices.
T          <- 300              # Number of simulated observations.

# Core stochastic block model parameters:
k          <- 5                # Number of blocks.
sc 		     <- rep(1,k)*(n/k)   # Sizes of communities.

# The probabilities depend on n. They are set so that they are 25% and 1% when n = 100.
p          <- 0.25 * (log(n)/log(100))^(1.01) / (n/100) # Probability of an edge.
q          <- 0.01 * (log(n)/log(100))^(1.01) / (n/100) # Probability of edges between blocks.

# Partial correlation network parameters:
sigma2     <- 1                # Reference variance of the data.
phi        <- 5                # Network dependence, in the interval: [0, inf).

########################################################################################################
# Part 1.1: The functions                                                                              #
########################################################################################################

# A function to simulate a stochastic block model random graph:
sbm <- function(Z, B){
  
  # Preliminaries:
  n  <- dim(Z)[1]
  rs <- 0
  
  # We draw a graph until there are no zero degree vertices:
  while(any(rs == 0)){
    
    # Population adjacency matrix:
    A_cal <- Z %*% B %*% t(Z)
    
    # Random adjacency matrix:
    A <- (replicate(n, runif(n, 0, 1)) < A_cal) * 1 
    
    # Making the adjacency matrix symmetric:
    A[lower.tri(A)] <- t(A)[lower.tri(A)]
    
    # Removing loops:
    diag(A) <- 0
    
    # The row sums to check for zero degree vertices:
    rs <- rowSums(A)
  }
  
  # The degree matrix:
  D <- diag(colSums(A))
  
  # The output as a list:
  output <- list("A" = A, "D" = D)
}

# A function to calculate a simple covariance matrix from the random graph. For more information,
# see the paper.
pcnm <- function(A, D, sigma2, phi){
  
  # Creating the normalised Laplacian:
  Dsq <- diag(1/sqrt(diag(D)))
  L   <- Dsq%*%(D-A)%*%Dsq  
  
  # The concentration matrix:
  K   <- (1/sigma2) * diag(replicate(dim(A)[1], 1)) + (phi/sigma2) * L
  
  # The covariance matrix:
  S   <- solve(K)
}

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
# Part 2: Simulating data                                                                              #
########################################################################################################

# Community membership matrix:
Z          <- matrix(,n,k)        
csc        <- c(0, cumsum(sc))
for(i in 1:k){
  z                     <- rep(0,k)
  z[i]                  <- 1
  Z[(csc[i]+1):csc[i+1],] <- t(replicate(sc[i],z))
}

# Matrix of edge probabilities:
B       <- replicate(k, rep(q, k)) 
diag(B) <- p

# Drawing a stochastic block model random graph:
block   <- sbm(Z, B)

# Creating the covariance matrix from the random graph:
S       <- pcnm(block$A, block$D, sigma2, phi)

# Simulating from a multivariate Gaussian:
Y       <- mvrnorm(T, replicate(n, 0), S)

########################################################################################################
# Part 3: Plotting the actual network                                                                  #
########################################################################################################

# Creating a graph object with the igraph library:
g    <- graph.adjacency(block$A, mode = 'undirected')
d    <- degree(g)

# Colouring the actual clusters:
colour <- rainbow(k)

clust <- 0
for(i in 1:k){
  clust[(csc[i]+1):csc[i+1]] <- i
}
for(i in 1:k){
  V(g)$color[clust == i] <- colour[i]
}

# Removing labels, setting vertex sizes and layout for the plot:
V(g)$label <- NA
V(g)$size  <- 5#round((d/max(d))*5+2)
loc        <- layout.fruchterman.reingold(g)

# Plotting:
par(mfrow=c(1,2))
plot(g, edge.color = "#808080", edge.width = 1.5, layout = loc, main = "Actual partition")

########################################################################################################
# Blockbuster clustering                                                                               #
########################################################################################################

# Obtaining the estimated partition indices from the Blockbuster algorithm:
ind.hat       <- blockbuster(Y, k)

# Colouring by the clustering results:
for( i in 1:length(unique(ind.hat))){
  V(g)$color[ind.hat==i] <- colour[i]
}

# Plotting with the same layout:
plot(g, layout = loc, main = "Estimated partition")