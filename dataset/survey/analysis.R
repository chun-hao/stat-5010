source("dataset/survey/survey_Missing.R")
survey_data <- read.csv("dataset/survey/survey_data.csv",header = TRUE)

attach(survey_data)
N1 <- cbind(round(n*p_1), round(n*p_2))
N1 <- cbind(N1, n - N1[,1] - N1[,2])
N2 <- cbind(round(n*q_1), round(n*q_2))
N2 <- cbind(N2, n - N2[,1] - N2[,2])
detach(survey_data)
Y <- cbind(N2, N1[,1:2])

# single survey
theta_sample <- single_survey_sampler(m, Y[1,])

apply(theta_sample, 2, mean) |> 
    theta_to_prob() |> 
    matrix(nrow = 3, ncol = 3, byrow = TRUE) |> 
    round(3)


#multiple surveys
model <- stan_model("dataset/survey/multi_survey.stan")
samples <- multi_survey_sampler(Y, model)
apply(samples$theta, c(2,3), mean)

apply(samples$nu, 2, \(x) mean(100*x)) |> 
    theta_to_prob() |>
    matrix(nrow = 3, ncol = 3, byrow = TRUE) |>
    round(2)

