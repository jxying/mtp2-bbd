library(Matrix)

check_stopping_fpn <- function(opts, x, y, iter, t0){
  if(!'max_time' %in% names(opts)){
    opts$max_time <- 1e4
  }
  if(!'max_iter' %in% names(opts)){
    opts$max_iter <- 1e3
  }
  if(!'tol' %in% names(opts)){
    opts$tol <- 1e-6
  }
  
  stop <- 0
  converge <- 0
  
  if (Sys.time() - t0 > opts$max_time){
    stop <- 1  # Maximum CPU time exceeded
  }
  if (iter > opts$max_iter){
    stop <- 1  # Maximum iterations exceeded
  }
  if (norm(x - y, "F") / norm(y, "F") < opts$tol){
    stop <- 1  # Condition on successive iterations holds
    converge <- 1
  }
  
  return(list(stop = stop, converge = converge))
}

chol_decomposition <- function(X){
  if(all(eigen(X)$values > 0)){
    L <- chol(X)
    flag_pd <- 0
  } else {
    svd_decomp <- svd(X)
    rank <- sum(svd_decomp$d > 1e-10)
    flag_pd <- rank + 1
    L <- NULL
  }
  return(list(L = L, flag_pd = flag_pd))
}

objective_function <- function(T, S){
  tryCatch({
    A <- chol(T)
    FLAG <- 0
    lg_det <- 2 * sum(log(diag(A)))
    fun <- -1 * lg_det + sum(T * S)
    return(list(value = fun, flag = FLAG))
  }, error = function(e) {
    return(list(value = Inf, flag = 1))
  })
}

solver_fpn <- function(S, lmbd, opts=NULL, hist_num=NULL) {
  t0 <- Sys.time()
  
  p <- ncol(S)
  iter <- 0
  delta <- 1e-15
  alpha <- 0.5
  
  if(length(dim(lmbd)) > 1) {
    Lamb <- lmbd
  } else {
    Lamb <- lmbd * (matrix(1, nrow=p, ncol=p) - diag(p))
  }
  
  if(is.null(hist_num)) {
    hist_num <- 0
  }
  
  if(is.null(opts)) {
    opts <- list()
  }
  if(!"max_iter" %in% names(opts)) opts$max_iter <- 1e4
  if(!"max_time" %in% names(opts)) opts$max_time <- 1e4
  if(!"tol" %in% names(opts)) opts$tol <- 1e-12
  beta <- ifelse("beta" %in% names(opts), opts$beta, 0.5)
  display <- ifelse("display" %in% names(opts), opts$display, 1)
  
  flag_edge <- "edge" %in% names(opts)
  edgeset <- which(!as.logical(opts$edge))
  
  flag_opt <- "X_opt" %in% names(opts)
  if(flag_opt) {
    X_opt <- opts$X_opt
    f_opt <- objective_function(X_opt, S - Lamb)$value
    relobj_iter <- vector()
    relerr_iter <- vector()
  }
  
  proj <- function(T) {
    pmin(0, T - diag(diag(T))) + diag(diag(T))
  }
  
  X <- diag(1/diag(S))
  Fcur <- objective_function(X, S - Lamb)
  if(Fcur$flag) {
    X <- diag(1 / (diag(S) + rep(1, p) * 1e-3))
    Fcur <- objective_function(X, S - Lamb)
  }
  objcur <- Fcur$value
  if(hist_num > 0) {
    Theta_hist <- array(0, c(p, p, 200))
    Theta_hist[,,1] <- X
  }
  obj_iter <- c(objcur)
  time_iter <- c(Sys.time() - t0)
  
  if(flag_opt) {
    rel_object <- abs(f_opt - objcur) / abs(f_opt)
    relobj_iter <- c(relobj_iter, rel_object)
    rel_err <- norm(X_opt - X, type="F") / norm(X_opt, type="F")
    relerr_iter <- c(relerr_iter, rel_err)
  }
  
  grad <- function(T) {
    -solve(T) + S - Lamb
  }
  gradf <- grad(X)
  
  check <- check_stopping_fpn(opts, X, X + matrix(1e16, nrow=p, ncol=p), iter, t0)
  
  if(check$stop) {
    X_new <- X
    objnew <- objcur
  }
  
  while(!check$stop) {
    iter <- iter + 1
    
    rstset <- (X - diag(diag(X)) - diag(rep(1e3, p)) > -delta) & (gradf < 0)
    
    if(flag_edge) {
      rstset <- union(rstset, edgeset)
    }
    
    X_up <- X
    X_up[rstset] <- 0
    
    grady <- gradf
    grady[rstset] <- 0
    
    descent <- X_up %*% grady %*% X_up
    descent[rstset] <- 0
    
    step_size <- 1
    
    Theta_f <- function(gamma) {
      proj(X_up - gamma * descent)
    }
    X_new <- Theta_f(step_size)
    
    list_A_flag_pd <- chol_decomposition(X_new)
    A <- list_A_flag_pd[[1]]
    flag_pd <- list_A_flag_pd[[2]]
    if(!flag_pd) {
      lg_det <- 2 * sum(log(diag(A)))
      objnew <- -1 * lg_det + sum(X_new * (S - Lamb))
    } else {
      objnew <- 1e8
    }
    
    gd <- abs(sum(gradf * descent))
    gdI <- abs(sum(gradf[rstset] * X[rstset]))
    
    while ((objnew > objcur - alpha * step_size * gd - alpha * gdI) && step_size > .Machine$double.eps) {
      step_size <- step_size * beta
      X_new <- Theta_f(step_size)
      
      list_A_flag_pd <- chol_decomposition(X_new)
      A <- list_A_flag_pd[[1]]
      flag_pd <- list_A_flag_pd[[2]]
      if(!flag_pd) {
        lg_det <- 2 * sum(log(diag(A)))
        objnew <- -1 * lg_det + sum(X_new * (S - Lamb))
      } else {
        objnew <- 1e8
      }
    }
    
    check <- check_stopping_fpn(opts, X_new, X, iter, t0)
    
    obj_iter <- c(obj_iter, objnew)
    time_iter <- c(time_iter, Sys.time() - t0)
    
    if(hist_num > 0 && iter < hist_num) {
      Theta_hist[,,length(obj_iter)] <- X_new
    }
    
    if(flag_opt) {
      rel_object <- abs(f_opt - objnew) / abs(f_opt)
      relobj_iter <- c(relobj_iter, rel_object)
      rel_err <- norm(X_opt - X_new, type="F") / norm(X_opt, type="F")
      relerr_iter <- c(relerr_iter, rel_err)
    }
    
    if(display && iter %% 5 == 0) {
      cat(paste("iter:", iter, "objective:", sprintf("%.11f", obj_iter[length(obj_iter)]), "cpu time:", sprintf("%.11f", time_iter[length(time_iter)])), "\n")
    }
    
    gradf <- grad(X_new)
    X <- X_new
    objcur <- objnew
  }
  
  run_time <- Sys.time() - t0
  
  out <- list()
  if(flag_opt) {
    out$relobj_itr <- relobj_iter
    out$relerr_itr <- relerr_iter
  }
  
  out$time <- run_time
  out$X_est <- X_new
  out$objective <- objnew
  out$obj_itr <- obj_iter
  out$time_itr <- time_iter
  out$iterate <- iter
  out$converge <- check$converge
  
  if(hist_num > 0) {
    out$Theta_hist <- Theta_hist[,,1:length(obj_iter)]
  }
  return(out)
}


