data {
  int<lower=0> n;
  int<lower=0, upper=1> y[n];
  vector[n] x;
}
parameters {
  real beta0;
  real beta1;
}
transformed parameters {
  vector[n] prob;
  vector[n] eta = beta0 + beta1*x;

  for (i in 1:n) {
    if(eta[i] > 0){
      prob[i] = 1 - (0.5/ (1 + eta[i]));
    }else{
      prob[i] = (0.5 / (1 - eta[i]));
    }
  }
}
model {
  beta0 ~ normal(0,100);
  beta1 ~ normal(0,100);
  y ~ bernoulli(prob);
}

generated quantities {
  vector[n] log_lik;
  for (i in 1:n) {
    log_lik[i] = bernoulli_lpmf(y[i] | prob[i]);
  }
}
