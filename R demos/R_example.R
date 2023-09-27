library(igraph) 
p <- 2000 # problem dimension

# Build a BA graph
BA_graph <- barabasi.game(p,  directed = FALSE)

# Generate data matrix, if you've already obtained sample covariance
# matrix, please skip this step

adjacency_matrix <- as_adjacency_matrix(BA_graph,  type = c("both"))
max_eig          <- eigen(adjacency_matrix)$values[1]
A                <- 1.05*max_eig*diag(p) - adjacency_matrix
inv_A            <- solve(A)
D                <- diag(sqrt(diag(inv_A)))
Mtrue            <- D %*% A %*% D
Ratio            <- 5
X                <- MASS::mvrnorm(Ratio * p , mu = rep(0, p), Sigma = solve(Mtrue))
 
# Compute the sample covariance matrix
S <- cov(X)
  
# Compute the regularization matrix
Theta0_mat <- matrix(0, nrow = p, ncol = p)
for (i in 1:p) {
  for (j in 1:p) {
    if (S[i, j] > 0 && i != j) {
      Theta0_mat[i, j] <- -S[i, j] / (S[i, i] * S[j, j] - S[i, j] * S[i, j])
    }
  }
}

chi          <- 0.01
Lambda       <- chi / (abs(Theta0_mat) + 0.0001)
diag(Lambda) <- 0
S_lambda     <- S - Lambda
 

# Conduct bridge-block decomposition

## Get the thresholded graph

S_res               <- S_lambda
S_res[S_lambda < 0] <- 0
diag(S_res)         <- 0
S_supp              <- as.matrix(S_res > 0)
G_thresholded       <- graph_from_adjacency_matrix(S_supp, mode = "undirected")

## Compute set of bridges in the thresholded graph
 
comps <- biconnected_components(G_thresholded)
bridges <- comps$component_edges[ lapply(comps$component_edges, length) == 1 ]

## Removing bridges

G_reduced <- G_thresholded
G_reduced <- delete_edges(G_reduced, unlist(bridges))

## Get the indexes for sub-graphs containing more than one node

Final_components <- clusters(G_reduced)$membership
Final_components_csize <- table(Final_components)

Subgraphs_index <- which(Final_components_csize > 1)
Subgraph_list <- vector("list", length(Subgraphs_index))

kit <- 1
for (k_index in Subgraphs_index) {
  nodes_indexes <- which(Final_components == k_index)
  Subgraph_list[[kit]] <- nodes_indexes
  kit <- kit + 1
}

# Solve sub-problems individually

library(R.matlab)

solve_using_matlab <- function(S, lmd) {
  n <- nrow(S)
  cur_wd <- getwd()
  setwd('FPN/functions')
  write.table(S, file='S.csv', sep=",", row.names=FALSE, col.names=FALSE)
  write.table(lmd, file='lmd.csv', sep=",", row.names=FALSE, col.names=FALSE)
  
  matlabCode <- "
  S = csvread('S.csv');
  lmd = csvread('lmd.csv');
  opts_FPN.max_iter = 1e4;
  opts_FPN.tol = 1e-10;
  out_FPN = solver_fpn(S, lmd, opts_FPN, 0);
  csvwrite('Theta_opt.csv', out_FPN.X_est);
  "
  
  writeLines(matlabCode, "matlab_code.m")
  
  system("matlab -nodisplay -r \"run('matlab_code.m'); exit\"")  # Run MATLAB CLI command
  
  while (!file.exists("Theta_opt.csv")) {
    Sys.sleep(1)  # Sleep for 1 second before checking again
  }
  
  Theta_opt <- read.csv("Theta_opt.csv", header = FALSE)
  Theta_opt <- matrix(unlist(Theta_opt), n, n)
  file.remove("S.csv", "lmd.csv", "Theta_opt.csv")
  setwd(cur_wd)
  return(Theta_opt)
}

out_FPN_sub_opt <- vector("list", length(Subgraphs_index))
for (i in 1:length(Subgraphs_index)) {
  sub_index <- Subgraph_list[[i]]
  S_sub <- S[sub_index, sub_index]
  Lambda_sub <- Lambda[sub_index, sub_index]
  out_FPN_sub_opt[[i]] <- solve_using_matlab(S_sub, Lambda_sub)
}

# Solve with FPN with bridge-block decomposition 
 
Theta_hat <- matrix(0, p, p)

Edge_array <- as_edgelist(G_thresholded)

for (e in 1:nrow(Edge_array)) {
  i <- Edge_array[e, 1]
  j <- Edge_array[e, 2]
  if (Final_components[i] != Final_components[j]) {
    Theta_hat[i, j] <- -(S_lambda[i, j]) / (S[i, i] * S[j, j] - (S_lambda[i, j] * S_lambda[i, j]))
    Theta_hat[j, i] <- Theta_hat[i, j]
  }
}

for (i in 1:p) {
  if (Final_components_csize[Final_components[i]] == 1) {
    Theta_ii <- 1
    for (j in neighbors(G_thresholded, i)) {
      Theta_ii <- Theta_ii + (S_lambda[i, j] * S_lambda[i, j]) / (S[i, i] * S[j, j] - S_lambda[i, j] * S_lambda[i, j])
    }
    Theta_ii <- Theta_ii / S[i, i]
    Theta_hat[i, i] <- Theta_ii
  } else {
    Theta_hat[i, i] <- 1 / S[i, i]
  }
}

for (k in which(Final_components_csize > 1)) {
  k_id <- which(Subgraphs_index == k)
  sub_i <- Subgraph_list[[k_id]]
  Theta_sub <- out_FPN_sub_opt[[k_id]] 
  Theta_hat[sub_i, sub_i] <- Theta_sub
  
  for (i in sub_i) {
    for (j in neighbors(G_thresholded, i)) {
      if (Final_components[i] != Final_components[j]) {
        Theta_hat[i, i] <- Theta_hat[i, i] + (1 / S[i, i]) * (S_lambda[i, j] * S_lambda[i, j]) / (S[i, i] * S[j, j] - S_lambda[i, j] * S_lambda[i, j])
      }
    }
  }
}

# Final objective value

print(-log(det(Theta_hat)) +sum((Theta_hat * S)) +sum(Lambda * abs(Theta_hat)))




# Solve with FPN without bridge-block decomposition

Theta_opt <- solve_using_matlab(S, Lambda)

print(-log(det(Theta_opt)) +sum((Theta_opt * S)) +sum(Lambda * abs(Theta_opt)))

