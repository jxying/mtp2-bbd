# mtp2-bbd
Matlab implementation of fast projected Newton-like (FPN) method [1] for learning large-scale MTP2 Gaussian graphical models and its accleration by introducing bridge-block decomposition [2]. The problem can be formulated as

$$
\underset{\boldsymbol{\Theta}}{\mathsf{minimize}}  -\log\det\left(\boldsymbol{\Theta}\right)+\left\langle \boldsymbol{\Theta},\mathbf{S}\right\rangle +\sum_{i\neq j}\Lambda_{ij}\left|\Theta_{ij}\right|, 
$$

subject to  

$$ 
	\boldsymbol{\Theta}\succ\mathbf{0}, \text{ and } \Theta_{ij}\leq0,\forall i\neq j
$$ 

The codes contain following procedures.

(1) Generating the data.

(2) Computing thresholded graph and bridge-block decomposition.

(3) Solving sub-problems individually using FPN solver [1].

(4) Obtaining optimal solution using methods in [2].

Please skip first step if you have data matrix or sample covariance matrix provided. The methods could significantly accelerate the convergence of existing algroithms and reduce memory cost when the thresholded graphs are sparse. 

## Simple Usage

Use fast projected Newton-like method:

```
opts_FPN = {'max_iter': 1e4, 'tol': 1e-10}
Theta = solver_fpn(S, Lambda, opts_FPN)
```

Use bridge-block decomposition approach:

```
Theta = solver_bbd(S, Lambda)
```
 
## References

[1] J.-F. Cai, J. V. de Miranda Cardoso, D. P. Palomar, and J. Ying, "Fast Projected Newton-like Method for Precision Matrix Estimation under Total Positivity", Neural Information Processing Systems (NeurIPS), New Orleans, LA, USA, Dec. 2023.

[2] X. Wang, J. Ying, and D. P. Palomar, "Learning Large-Scale MTP2 Gaussian Graphical Models via Bridge-Block Decomposition," Neural Information Processing Systems (NeurIPS), New Orleans, LA, USA, Dec. 2023.







