rm(list=ls())
set.seed(10)
res <- .generate_A(1000, epsilon = 0.1, lambda = 2)

A <- res$A
threshold <- max(res$x)/5
amp_res <- amp(A, threshold, iter_max = 100, verbose = T)

plot(amp_res[,1], amp_res[,100], asp = T)
plot(res$x, amp_res[,1], asp = T)
plot(1:ncol(amp_res), colSums(abs(amp_res)), ylim = c(0, 10*sum(abs(res$x))))

######################

param <- .state_initialize_param(A, lambda = 2)
x_param <- .generate_x(epsilon = 0.1)
state_res <- state_evolution(A, threshold, lambda = 2, x_param)

state_evaluation_function(x_param, state_res["mu",1000], state_res["sigma",1000])
