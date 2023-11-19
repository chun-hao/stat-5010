data {
    int<lower=0> J;   
    int<lower=0> N1[J, 3];
    int<lower=0> N2[J, 3];
}

parameters {
  simplex[3] p1[J];
  simplex[3] p2[J];
  simplex[3] nu[2];
  real<lower=0> alpha;
}
transformed parameters {
  vector<lower=0>[3] xi[2];
  xi[1] = alpha * nu[1];
  xi[2] = alpha * nu[2];
}
model {
    //alpha ~ cauchy(0,1);
    alpha ~ exponential(0.001);
    //nu ~ dirichlet([4/15, 4/15, 2/15, 2/15, 1/5]');
    nu[1] ~ dirichlet(rep_vector(1, 3));
    nu[2] ~ dirichlet(rep_vector(1, 3));
    for (j in 1:J) {
        p1[j] ~ dirichlet(xi[1]);
        p2[j] ~ dirichlet(xi[2]);
    }
    for (j in 1:J) {
        N1[j] ~ multinomial(p1[j]);
        N2[j] ~ multinomial(p2[j]);
    }
}
