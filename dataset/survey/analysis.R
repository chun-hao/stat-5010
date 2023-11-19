survey_data <- read.csv("../dataset/survey/survey_data.csv",header = TRUE)

attach(survey_data)
N1 <- cbind(round(n*p_1), round(n*p_2))
N1 <- cbind(N1, n - N1[,1] - N1[,2])
N2 <- cbind(round(n*q_1), round(n*q_2))
N2 <- cbind(N2, n - N2[,1] - N2[,2])
detach(survey_data)
Y <- cbind(N1, N2[,1:2])

out <- multi_survey_sampler(Y)

## test
set.seed(123)
theta <- rdirichlet(1, rep(1, 6))
U <- rmultinom(1, 2000, prob = theta_to_prob(theta))
Y <- (M %*% U)[1:5]
Z <- (M %*% U)[6:9]
#Y <- c(988, 802, 256, 943, 851)


theta_sample <- single_survey_sampler(m, Y)

apply(theta_sample, 2, mean) |> 
    theta_to_prob() |> 
    matrix(nrow = 3, ncol = 3, byrow = TRUE) |> 
    round(3)

theta_to_prob(theta) |> 
    matrix(nrow = 3, ncol = 3, byrow = TRUE) |> 
    round(3)

model <- stan_model("survey.stan")

theta_sample <- array(0, dim = c(40000, 6, 6))
nu_sample <- array(0, dim = c(40000, 6))

for (i in 1:length(out)){
    theta_sample[((i-1)*4000+1):(i*4000), , ] <- out[[i]]$theta
    nu_sample[((i-1)*4000+1):(i*4000), ] <- out[[i]]$nu
}


apply(theta_sample, c(2,3), mean) %>%
    apply(., 2, theta_to_prob) |> 
    matrix(nrow = 3, ncol = 3, byrow = TRUE) |> 
    round(3)
apply(theta_sample, c(2,3), sd)

apply(nu_sample, 2, mean) |> 
    theta_to_prob() |>
    matrix(nrow = 3, ncol = 3, byrow = TRUE) |>
    round(3)
