data {
  int<lower=1> n_replicates; // number of replicates for each cytosine, name is kept compact
  int<lower=1> n_predictors; // number of predictors

  real<lower=0,upper=1> bsEff[n_replicates];
  real<lower=0,upper=1> bsBEff[n_replicates];
  real<lower=0,upper=1> seqErr[n_replicates];

  int<lower=0> bsC[n_replicates];
  int<lower=0> bsTot[n_replicates];

  matrix[n_replicates, n_predictors] X;

  real<lower=0> sigmaB2;

  real<lower=0> alpha;
  real<lower=0> beta;

}
parameters {
  vector[n_predictors] B;
  real<lower=0> sigmaE2;
  vector[n_replicates] Y;
  real<lower=0> sigmaR2;
  
}

transformed parameters {

  vector<lower=0,upper=1>[n_replicates] theta;

  for (i in 1:n_replicates){
    theta[i] = inv_logit(Y[i]);
  }
}

model {
  vector[n_replicates] mu;
  sigmaE2 ~ gamma(alpha,beta);
  B ~ normal(0,sigmaB2);
  mu = X * B;

  Y ~ normal(mu, sigmaE2);

  for (i in 1:n_replicates){
    bsC[i] ~ binomial(bsTot[i],
            theta[i] * ((1.0 - seqErr[i]) * (1.0 - bsEff[i]) + seqErr[i] * bsEff[i]) +
            (1-theta[i]) * ((1.0 - bsBEff[i]) * (1.0 - seqErr[i]) + seqErr[i] * bsBEff[i]));
  }
}
