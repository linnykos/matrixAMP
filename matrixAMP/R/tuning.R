.s_star <- function(x_param, gamma, theta, trials = 5000){
  x1 <- .sample_x(trials, x_param)
  x2 <- .sample_x(trials, x_param)
  g1 <- stats::rnorm(trials)
  g2 <- stats::rnorm(trials)
  
  num <- (mean(x1 * .soft_threshold(sqrt(gamma)*x1 + g1, theta)))^2
  denom <- mean((.soft_threshold(sqrt(gamma)*x2 + g2, theta))^2)
  
  num/denom
}

# we're deviating a bit 
.s <- function(gamma, theta, p_num = 21, trials = 5000){
  p_vec <- seq(0, 1, length.out = p_num+2)
  p_vec <- p_vec[c(2:(length(p_vec)-1))]
  val <- sapply(p_vec, function(p){
    x_param <- .generate_x(epsilon = p)
    .s_star(x_param, gamma, theta, trials)
  })
  
  min(val)
}

.gamma_update <- function(lambda, gamma, theta, p_num = 21, trials = 5000){
  lambda^2 * .s(gamma, theta, p_num, trials)
}

.theta_update <- function(gamma, max_val = 5, theta_num = 11, p_num = 21, trials = 5000){
  theta_val <- seq(0, max_val, length.out = theta_num)
  val <- sapply(theta_val, function(theta){
    .s(gamma, theta, p_num, trials)
  })
  
  max(val)
}