import igraph as ig
import numpy as np
from scipy.stats import multivariate_normal
from FPN import solver_fpn, objective_function

# Set problem dimension
p = 2000

# Build a BA graph
# If you've already have the data matrix, skip this step

BA_graph = ig.Graph.Barabasi(n=p, m=1, directed=False)

adjacency_matrix = BA_graph.get_adjacency()
adjacency_matrix = np.array(adjacency_matrix.data)

max_eig = np.max(np.linalg.eigvals(adjacency_matrix)).real
A = 1.05 * max_eig * np.eye(p) - adjacency_matrix
inv_A = np.linalg.inv(A)
D = np.diag(np.sqrt(np.diag(inv_A)))
Mtrue = D @ A @ D
Ratio = 5
X = multivariate_normal.rvs(np.zeros(p), np.linalg.inv(Mtrue), size=Ratio * p)

# Compute the sample covariance matrix
S = np.cov(X, rowvar=False)

# Compute the regularization matrix
Theta0_mat = np.zeros((p, p))
for i in range(p):
    for j in range(p):
        if S[i, j] > 0 and i != j:
            Theta0_mat[i, j] = -S[i, j] / (S[i, i] * S[j, j] - S[i, j] * S[i, j])

chi = 0.02
Lambda = chi / (np.abs(Theta0_mat) + 0.0001)
np.fill_diagonal(Lambda, 0)

def solve_using_bbd(S, Lambda):
    S_lambda = np.maximum(0, S - Lambda)  # Thresholded matrix

    # Conduct bridge-block decomposition

    # Get the thresholded graph
    S_res = np.copy(S_lambda)
    S_res[S_lambda < 0] = 0
    np.fill_diagonal(S_res, 0)
    S_supp = S_res > 0
    G_thresholded = ig.Graph.Adjacency(S_supp.tolist(), mode="undirected")

    # Compute set of bridges in the thresholded graph
    bridges = G_thresholded.bridges()

    # Remove bridges from the graph
    G_reduced = G_thresholded.copy()
    for bridge in reversed(bridges):
        G_reduced.delete_edges(bridge)

    # Get the indexes for sub-graphs containing more than one node
    Final_components = G_reduced.connected_components().membership
    Final_components_csize = np.bincount(Final_components)

    Subgraphs_index = np.where(Final_components_csize > 1)[0]
    Subgraph_list = [np.where(Final_components == k_index)[0] for k_index in Subgraphs_index]

    # Solve sub-problems individually using FPN

    out_FPN_sub_opt = []
    for i in range(len(Subgraphs_index)):
        sub_index = Subgraph_list[i]
        S_sub = S[sub_index][:, sub_index]
        Lambda_sub = Lambda[sub_index][:, sub_index]
        out_FPN_sub_opt.append(solver_fpn(S_sub, Lambda_sub))

    # Ontain optimal solution

    p = S.shape[0]
    Theta_hat = np.zeros((p, p))

    Edge_array = np.array(G_thresholded.get_edgelist())

    for e in range(Edge_array.shape[0]):
        i = Edge_array[e, 0]
        j = Edge_array[e, 1]
        if Final_components[i] != Final_components[j]:
            Theta_hat[i, j] = -(S_lambda[i, j]) / (S[i, i] * S[j, j] - (S_lambda[i, j] ** 2))
            Theta_hat[j, i] = Theta_hat[i, j]

    for i in range(p):
        if Final_components_csize[Final_components[i]] == 1:
            Theta_ii = 1
            for j in G_thresholded.neighbors(
                    i):  # Assuming neighbors is a function returning neighbors of node i in graph G_res
                Theta_ii += (S_lambda[i, j] ** 2) / (S[i, i] * S[j, j] - S_lambda[i, j] ** 2)
            Theta_ii /= S[i, i]
            Theta_hat[i, i] = Theta_ii
        else:
            Theta_hat[i, i] = 1 / S[i, i]

    k_indices = np.where(Final_components_csize > 1)[0]
    for k in k_indices:
        k_id = np.where(Subgraphs_index == k)[0][0]
        sub_i = Subgraph_list[k_id]
        Theta_sub = out_FPN_sub_opt[k_id]["X_est"]
        Theta_hat[np.ix_(sub_i, sub_i)] = Theta_sub

        for i in sub_i:
            for j in G_thresholded.neighbors(i):
                if Final_components[i] != Final_components[j]:
                    Theta_hat[i, i] += (1 / S[i, i]) * (S_lambda[i, j] ** 2) / (S[i, i] * S[j, j] - S_lambda[i, j] ** 2)

    return Theta_hat


Theta_bbd = solve_using_bbd(S, Lambda)
obj_FPN_bbd = objective_function(Theta_bbd, S - Lambda)["value"]
print(obj_FPN_bbd)


# Use FPN without birgde-block decomposition

opts_FPN = {'max_iter': 1e4, 'tol': 1e-10}
out_FPN = solver_fpn(S, Lambda, opts_FPN, 0)
print(objective_function(out_FPN["X_est"], S - Lambda)["value"])

