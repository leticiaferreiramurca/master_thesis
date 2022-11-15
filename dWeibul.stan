//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n;
  int<lower=0, upper=1> y[n];
  vector[n] x;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.

parameters {
  real<lower=-3, upper=3> logc;
  real beta0;
  real beta1;
}


transformed parameters {
  vector[n] eta = beta0 + beta1*x;
  vector[n] prob;
  real c = exp(logc);
  
  for (i in 1:n) {
    if(eta[i] > 0){
      prob[i] = 1 - 0.5 * exp(- eta[i]^c);
    }else{
       prob[i] = 0.5*exp(- fabs(eta[i])^c);
    }
  }
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  logc ~ uniform(-3,3);
  beta0 ~ normal(0.0,100);
  beta1 ~ normal(0.0,100);
  y ~ bernoulli(prob);
}
