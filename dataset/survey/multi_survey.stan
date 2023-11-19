functions {
  vector theta_to_prob(vector theta) {
    vector[9] prob;
    prob[1] = theta[1];
    prob[2] = theta[4];
    prob[3] = theta[6]/4;
    prob[4] = theta[5];
    prob[5] = theta[2];
    prob[6] = theta[6]/4;
    prob[7] = theta[6]/4;
    prob[8] = theta[6]/4;
    prob[9] = theta[3];
    return prob;
  }
}

data {
    int<lower=0> J;   
    int<lower=0> U[J, 9];
}

parameters {
    simplex[6] theta[J];
    simplex[6] nu;
    real<lower=0> alpha;
}

transformed parameters {
    vector<lower=0>[9] prob[J];
    for (j in 1:J) {
        prob[j] = theta_to_prob(theta[j]);
    }
    
    vector<lower=0>[6] xi;
    xi = alpha * nu;
}

model {
    alpha ~ exponential(0.001);
    nu ~ dirichlet(rep_vector(1, 6));
    for (j in 1:J) {
        theta[j] ~ dirichlet(xi);
    }
    for (j in 1:J) {
        U[j] ~ multinomial(prob[j]);
    }
}
