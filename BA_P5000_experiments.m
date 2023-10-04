%% Add paths
addpath('FPN/functions');

%% Read EdgeList from local file 
Edgelist = readtable('Data/EdgeList_BA5000.csv');
Edgelist = Edgelist{:,:};

%% Generating the data 
% If you have data generated, please skip this step
p = 5000;
adjacency_matrix = sparse(Edgelist(:,1), Edgelist(:,2), ones(size(Edgelist,1),1), p, p);
adjacency_matrix = adjacency_matrix + adjacency_matrix';
adjacency_matrix = full(adjacency_matrix);
eigs    = eig(adjacency_matrix);
max_eig = eigs(p);
A       = 1.05 * max_eig * diag(ones(p,1)) - adjacency_matrix;
inv_A   = inv(A);
D       = diag(sqrt(diag(inv_A)));
Mtrue   = D * A * D;
Ratio   = 10; 
X       = mvnrnd(zeros(p,1), inv(Mtrue), Ratio * p); % Data matrix
 
%% Compute the sample covariance matrix
S = cov(X); 
 
%% Compute the regularization matrix 
Theta0_mat = zeros(p,p);
for i = 1:p
    for j = 1:p
        if S(i,j) > 0 && i ~= j
            Theta0_mat(i,j) = -S(i,j)/(S(i,i) * S(j,j) - S(i,j) * S(i,j));
        end
    end
end

chi = 0.01; 
Lambda = chi./((abs(Theta0_mat) + 0.0001));
for i =1:p
    Lambda(i,i) = 0;
end 
S_lambda = S - Lambda;
 
%% Conduct bridge-block decomposition   
S_res      = max(0, S_lambda);
S_res      = S_res - diag(diag(S_res));
S_supp     = boolean(S_res); 
G_res      = graph(S_supp);
Edge_array = table2array(G_res.Edges); 
bins       = biconncomp(G_res);
comp_count = histc(bins, unique(bins)); 

s = zeros(1, length(bins));
t = zeros(1, length(bins));
counter = 0;

for i = 1:length(bins)
    if comp_count(bins(i)) == 1
        counter = counter + 1;
        s(counter) = Edge_array(i, 1);
        t(counter) = Edge_array(i, 2);
    end
end

s = s(1:counter);
t = t(1:counter);
G_res_reduced = rmedge(G_res,s,t);
Final_components = conncomp(G_res_reduced);
Final_components_csize =  histc(Final_components, unique(Final_components));

Subgraphs_index = find(Final_components_csize>1);
Subgraph_list = cell(length(Subgraphs_index),1);

kit = 1;
for k_index = Subgraphs_index
    nodes_indexes = find(Final_components == k_index);
    Subgraph_list{kit} = nodes_indexes; % Get the sub-problem with dimention greater than 1
    kit = kit+1;
end

%% Compute optimal solution of sub-problems

Theta_hat = zeros(p,p);

for e = 1:size(Edge_array,1)
   i = Edge_array(e,1);
   j = Edge_array(e,2);
   if Final_components(i) ~= Final_components(j)
       Theta_hat(i,j) = -(S_lambda(i,j))/(S(i,i) * S(j,j) - (S_lambda(i,j) * S_lambda(i,j)));
       Theta_hat(j,i) = Theta_hat(i,j); 
   end 
end

for i = 1:p
    if Final_components_csize(Final_components(i)) == 1
        Theta_ii = 1; 
        for j = neighbors(G_res,i)' 
            Theta_ii = Theta_ii + (S_lambda(i,j) * S_lambda(i,j))/(S(i,i) * S(j,j) - S_lambda(i,j) * S_lambda(i,j));
        end 
        Theta_ii = Theta_ii / S(i,i);
        Theta_hat(i,i) = Theta_ii;
    else 
        Theta_hat(i,i) = 1/S(i,i); 
    end
end

out_FPN_sub_opt = cell(length(Subgraphs_index),1);
for i = 1:length(Subgraphs_index)
    sub_index = Subgraph_list{i}; 
    S_sub = S(sub_index,sub_index);
    Lambda_sub = Lambda(sub_index,sub_index);
    opts_FPN.max_iter = 1e4; 
    opts_FPN.tol      = 1e-12; 
    out_FPN_sub_opt{i} = solver_fpn(S_sub, Lambda_sub, opts_FPN, 0);
end

%% Solve with FPN with bridge-block decomposition 
for k = find(Final_components_csize>1)
    k_id = find(Subgraphs_index == k);
    sub_i = Subgraph_list{k_id};
    Theta_sub =  out_FPN_sub_opt{k_id}.X_est;
    Theta_hat(sub_i,sub_i) = Theta_sub; 
        
    for i = sub_i 
        for j = neighbors(G_res,i)' 
            if Final_components(i) ~= Final_components(j) 
                Theta_hat(i,i) = Theta_hat(i,i) +(1/S(i,i))* (S_lambda(i,j) * S_lambda(i,j))/(S(i,i) * S(j,j) - S_lambda(i,j) * S_lambda(i,j));
            end
        end
    end
end

obj_FPN_bbd = objective_function(Theta_hat, S_lambda).value;
     
%% Solve with FPN without bridge-block decomposition

opts_FPN.max_iter = 1e4; 
opts_FPN.tol      = 1e-10;
out_FPN           = solver_fpn(S, Lambda, opts_FPN, 0);
disp(objective_function(out_FPN.X_est, S_lambda).value);

