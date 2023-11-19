library(rstan)
library(MCMCpack)

M <- matrix(c(1,1,1,0,0,0,0,0,0,
              0,0,0,1,1,1,0,0,0,
              0,0,0,0,0,0,1,1,1,
              1,0,0,1,0,0,1,0,0,
              0,1,0,0,1,0,0,1,0,
              1,0,0,0,0,0,0,0,0,
              0,1,0,0,0,0,0,0,0,
              0,0,0,1,0,0,0,0,0,
              0,0,0,0,0,0,0,1,0), nrow = 9, ncol = 9, byrow = TRUE)

M_inv <- solve(M)

theta_to_prob <- function(theta){
    return(c(theta[1], theta[4], theta[6]/4,
             theta[5], theta[2], theta[6]/4, 
             theta[6]/4, theta[6]/4, theta[3]))
}

prob_to_theta <- function(prob){
    return(c(prob[1], prob[5], prob[9], prob[2], prob[4], prob[3]+prob[6]+prob[7]+prob[8]))
}

Gibbs_Z <- function(n, Y, warmup = 1000, seed = NULL){
    # Z = [n_11, n_12, n_21, n_32]
    # Y = [n_1., n_2., n_3., n_.1, n_.2]
    Z <- matrix(0, nrow = warmup + n, ncol = 4)
    set.seed(seed)
    iter <- warmup + n
    for (i in 2:iter){
        Z[i, 1] <- rhyper(1, m = Y[1] - Z[i-1, 2], 
                          n = Y[3] - Z[i-1, 4],
                          k = Y[4] - Z[i-1, 3])
        Z[i, 2] <- rhyper(1, m = Y[1] - Z[i, 1],
                          n = Y[2] - Z[i-1, 3],
                          k= Y[5] - Z[i-1, 4])
        Z[i, 3] <- rhyper(1, m = Y[2] - Y[5] + Z[i, 2] + Z[i-1, 4], 
                          n = Y[3] - Z[i-1, 4],
                          k = Y[4] - Z[i, 1])
        Z[i, 4] <- rhyper(1, m = Y[3] - Y[4] + Z[i, 1] + Z[i, 3],
                          n = Y[2] - Z[i, 3],
                          k = Y[5] - Z[i, 2])
        
    }
    return(Z[warmup+1:n,])
}

single_survey_sampler <- function(m, Y, seed = 1234){
    # Y = [n_1., n_2., n_3., n_.1, n_.2]
    set.seed(seed)
    theta_sample <- matrix(0, nrow = m, ncol = 6) # (\rho_1, \rho_2, \rho_3, \rho_{12}, \rho_{21}, 4\gamma)
    
    n <- sum(Y[1:3])
    
    z_sample <- Gibbs_Z(m, Y)
    
    for(i in 1:m){
        u_new <- M_inv %*% c(Y, z_sample[i,])
        theta_new <- rdirichlet(1, rep(1, 6) + n * prob_to_theta(u_new/n))
        theta_sample[i, ] <- theta_new
    }
    colnames(theta_sample) <- c("rho1", "rho2", "rho3", "rho12", "rho21", "gamma")
    return(theta_sample)
}

multi_survey_sampler <- function(Y, seed = 1234){
    # Y is a k x 5 matrix where k is the number of surveys
    k <- nrow(Y)
    
    n_z <- 10
    set.seed(seed)
    z_sample <- array(0, dim = c(n_z, k, 4))
    U_sample <- array(0, dim = c(n_z, k, 9))
    for(i in 1:k){
        z_sample[, i, ] <- Gibbs_Z(n_z, Y[i, ])
        for(j in 1:n_z){
            U_sample[j, i, ] <- M_inv %*% c(Y[i, ], z_sample[j, i, ])
        }
    }
    
    model <- stan_model("multi_survey.stan")
    
    samples <- vector("list", n_z)
    for(i in 1:n_z){
        fit <- sampling(model,list(J = k, U = U_sample[i, , ]),
                        iter=2000, warmup=1000,
                        chains=4, cores = 4)
        samples[[i]] <- extract(fit, pars = c("theta", "nu", "alpha"))
    }
    return(samples)
}






