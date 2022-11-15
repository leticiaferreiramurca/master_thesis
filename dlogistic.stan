
data {
  int<lower=0> n;
  int<lower=0, upper=1> y[n];
  vector[n] x;
}
parameters {
  real beta0;
  real beta1;
}
model {
  beta0 ~ normal(0,100);
  beta1 ~ normal(0,100);
  y ~ bernoulli_logit(beta0 + beta1*x);
}
generated quantities {
  vector[n] log_lik;
  for (ni in 1:n) {
    log_lik[ni] = bernoulli_logit_lpmf(y[ni] | beta0 + x[ni] * beta1);
  }
}
