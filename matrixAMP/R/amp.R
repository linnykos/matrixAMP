# initialize x
.amp_initialize_x <- function(A){
  n <- nrow(A)
  n*as.numeric(mgcv::slanczos(A, k = 1)$vectors[,1])
}

# update x
.amp_update_x <- function(A, x, threshold, tol = 1e-6){
  stopifnot(length(x) == nrow(A), nrow(A) == ncol(A))
  n <- nrow(A)
  y <- .soft_threshold(x, threshold)
  b <- length(which(abs(x) >= tol))/n
  
  as.numeric(A %*% y - b * y)
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

.state_initialize_param <- function(A){
  lambda_max <- abs(mgcv::slanczos(A, k = 1)$values)
  lambda <- .5*(lambda_max + sqrt(lambda_max^2 - 4))
  
  mu = sqrt(1-lambda^(-2))
  sigma = 1/lambda
  
  list(lambda = lambda, mu = mu, sigma = sigma)
}

.state_update <- function(param, threshold, x_mat, trials = 5000){
  stopifnot(all(names(param) == c("lambda", "mu", "sigma")))
  mu <- .state_evaluate_mu_exp(param$lambda, param$mu, param$sigma, threshold, x_mat, trials)
  sigma <- .state_evaluate_sigma_exp(param$lambda, param$mu, param$sigma, threshold, x_mat, trials)
  
  list(lambda = param$lambda, mu = mu, sigma = sigma)
}

.state_evaluate_mu_exp <- function(lambda, mu, sigma, threshold, x_mat, trials = 5000){
  x <- .sample_x(trials, x_mat)
  g <- stats::rnorm(trials)
  lambda * mean(x * .soft_threshold(mu*x + sigma*g, threshold))
}

.state_evaluate_sigma_exp <- function(lambda, mu, sigma, threshold, x_mat, trials = 5000){
  x <- .sample_x(trials, x_mat)
  g <- stats::rnorm(trials)
  mean((.soft_threshold(mu*x + sigma*g, threshold))^2)
}