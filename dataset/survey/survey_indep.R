library(rstan)
library(bayesplot)
library(shinystan)


model <- stan_model("dataset/survey/survey_indep.stan")

survey_data <- read.csv("dataset/survey/survey_data.csv",header = TRUE)

attach(survey_data)
N1 <- cbind(round(n*p_1), round(n*p_2))
N1 <- cbind(N1, n - N1[,1] - N1[,2])
N2 <- cbind(round(n*q_1), round(n*q_2))
N2 <- cbind(N2, n - N2[,1] - N2[,2])
detach(survey_data)

fit <- sampling(model,list(J = 6, N1 = N1, N2 = N2),
               iter=5000, warmup=2000,
               chains=4, cores = 4)

launch_shinystan(fit)

summary(fit, pars = c("nu", "rho"))

p1_sample <- extract(fit, pars="p1")[[1]]
p1_post_mean <- apply(p1_sample, c(2,3), mean)
p1_post_sd <- apply(p1_sample, c(2,3), sd)
p2_sample <- extract(fit, pars="p2")[[1]]
p2_post_mean <- apply(p2_sample, c(2,3), mean)
p2_post_sd <- apply(p2_sample, c(2,3), sd)

# difference in proportions

T_1 <- matrix(0, nrow = nrow(p1_sample), ncol = 6)
for(i in 1:6){
    T_1[,i] <- p1_sample[,i,1] - p2_sample[,i,1]
}

plot(0, type="n", xlim = c(-0.2,0.2), ylim = c(0, 30))
for(i in 1:6){
    lines(density(T_1[,i]), col = i) 
}

dip_summary <- tibble(post_mean = 100*apply(T_1, 2, mean),
                      post_sd = 100*apply(T_1, 2, sd),
                      CI_lower = 100*apply(T_1, 2, quantile, 0.025),
                      CI_upper = 100*apply(T_1, 2, quantile, 0.975)) %>%
    round(.,2)

# difference in differences
T_2 <- matrix(0, nrow = nrow(p1_sample), ncol = 6)
for(i in 1:6){
    T_2[,i] <- (p1_sample[,i,1] - p1_sample[,i,2]) - 
        (p2_sample[,i,1] - p2_sample[,i,2])
}

plot(0, type="n", xlim = c(-0.15,0.25), ylim = c(0, 15))
for(i in 1:6){
    lines(density(T_2[,i]), col = i) 
}

did_summary <- tibble(post_mean = 100*apply(T_2, 2, mean),
                      post_sd = 100*apply(T_2, 2, sd),
                      CI_lower = 100*apply(T_2, 2, quantile, 0.025),
                      CI_upper = 100*apply(T_2, 2, quantile, 0.975)) %>%
    round(.,2)

## overall dip
nu_sample <- extract(fit, pars = "nu")[[1]]

T_3 <- nu_sample[,1,1] - nu_sample[,2,1]
T_4 <- (nu_sample[,1,1] - nu_sample[,1,2]) - (nu_sample[,2,1] - nu_sample[,2,2])
