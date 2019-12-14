amp <- function(A, threshold, iter_max = 1000, verbose = T){
  x_mat <- matrix(0, nrow = nrow(A), ncol = iter_max)
  x_mat[,1] <-.amp_initialize_x(A)
  
  for(iter in 2:iter_max){
    if(verbose && iter_max >= 10 && iter %% floor(iter_max/10) == 0) cat('*')
    x_mat[,iter] <- .amp_update_x(A, x = x_mat[,iter-1], x_prev = x_mat[,max(iter-2,1)], threshold)
  }
  
  x_mat
}

state_evolution <- function(A, threshold, lambda, x_param, iter_max = 1000, trials = 5000, verbose = T){
  param_mat <- matrix(0, nrow = 2, ncol = iter_max)
  rownames(param_mat) <- c("mu", "sigma")
  res_old <- .state_initialize_param(A, lambda)
  param_mat[,1] <- c(res_old$mu, res_old$sigma)
  
  for(iter in 2:iter_max){
    if(verbose && iter_max >= 10 && iter %% floor(iter_max/10) == 0) cat('*')
    
    res_new <- .state_update(res_old, threshold, x_param, trials)
    param_mat[,iter] <- c(res_new$mu, res_new$sigma)
  }
  
  param_mat
}

l2_evaluation <- function(vec1, vec2){
  (vec1 - vec2)^2
}

amp_evaluation_function <- function(x_est, x_pop, func = l2_evaluation){
  mean(func(x_est, x_pop))
}

state_evaluation_function <- function(x_param, mu, sigma, trials = 5000, func = l2_evaluation){
  x <- .sample_x(trials, x_param)
  g <- stats::rnorm(trials)
  mean(func(x, mu*x + sigma*g))
}

######

# initialize x
.amp_initialize_x <- function(A){
  n <- nrow(A)
  psi <- as.numeric(mgcv::slanczos(A, k = 1)$vectors[,1])
  x <- sqrt(n) * psi

  x
}

# update x
.amp_update_x <- function(A, x, x_prev, threshold, tol = 1e-6){
  stopifnot(length(x) == nrow(A), nrow(A) == ncol(A))
  n <- nrow(A)
  y <- .soft_threshold(x, threshold)
  y_prev <- .soft_threshold(x_prev, threshold)
  
  b <- length(which(abs(y) >= tol))/n
  
  as.numeric(A %*% y - b * y_prev)
}

# soft threshold, where x can be a vector and val is a scalar
.soft_threshold <- function(x, val){
  stopifnot(val >= 0)
  
  sapply(x, function(z){
    if(z >= val) return(z-val)
    if(z <= -val) return(z+val)
    0
  })
}

####

.state_initialize_param <- function(A, lambda){
  mu = sqrt(1-lambda^(-2))
  sigma = 1/lambda
  
  list(lambda = lambda, mu = mu, sigma = sigma)
}

.state_update <- function(param, threshold, x_param, trials = 5000){
  stopifnot(all(names(param) == c("lambda", "mu", "sigma")))
  mu <- .state_evaluate_mu_exp(param$lambda, param$mu, param$sigma, threshold, x_param, trials)
  sigma <- .state_evaluate_sigma_exp(param$lambda, param$mu, param$sigma, threshold, x_param, trials)
  
  list(lambda = param$lambda, mu = mu, sigma = sigma)
}

.state_evaluate_mu_exp <- function(lambda, mu, sigma, threshold, x_param, trials = 5000){
  x <- .sample_x(trials, x_param)
  g <- stats::rnorm(trials)
  lambda * mean(x * .soft_threshold(mu*x + sigma*g, threshold))
}

.state_evaluate_sigma_exp <- function(lambda, mu, sigma, threshold, x_param, trials = 5000){
  x <- .sample_x(trials, x_param)
  g <- stats::rnorm(trials)
  mean((.soft_threshold(mu*x + sigma*g, threshold))^2)
}

