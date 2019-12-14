# epsilon is the sparsity of the target distribution x_0
.generate_A <- function(n, epsilon, lambda){
  x_param <- .generate_x(epsilon)
  x <- .sample_x(n, x_param)
  A <- (lambda/n)*x %*% t(x) + .generate_W(n)
  
  list(A = A, x = x)
}

# generate a symmetric distribution of x taking values -a, 0, a such that the probability of 0
## is 1-epsilon and the second moment of the distribution is 1.
## this is represented as a 2x3 matrix
.generate_x <- function(epsilon = 0.1){
  x_param <- matrix(0, nrow = 2, ncol = 3)
  a <- sqrt(1/epsilon)
  x_param[1,] <- c(-a, 0, a)
  x_param[2,] <- c(epsilon/2, 1-epsilon, epsilon/2)
  
  rownames(x_param) <- c("value", "probability")
  x_param
}

.sample_x <- function(n, x_param){
  stopifnot(all(dim(x_param) == c(2,3)))
  sample(x_param[1,], size = n, replace = T, prob = x_param[2,])
}

.generate_W <- function(n){
  w_mat <- matrix(0, n, n)
  upper_tri_idx <- upper.tri(w_mat)
  
  w_mat[upper_tri_idx] <- stats::rnorm(sum(upper_tri_idx), mean = 0, sd = 1/sqrt(n))
  w_mat <- w_mat + t(w_mat)
  diag(w_mat) <- stats::rnorm(n, mean = 0, sd = sqrt(2/n))
  
  w_mat
}

.l2norm <- function(x){sqrt(sum(x^2))}

