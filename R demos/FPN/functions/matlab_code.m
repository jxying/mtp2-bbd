
  S = csvread('S.csv');
  lmd = csvread('lmd.csv');
  opts_FPN.max_iter = 1e4;
  opts_FPN.tol = 1e-10;
  out_FPN = solver_fpn(S, lmd, opts_FPN, 0);
  csvwrite('Theta_opt.csv', out_FPN.X_est);
  
