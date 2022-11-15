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
  real<lower=-2, upper=2> logc;
  real<lower=-2, upper=2> loglambda;
  real beta0;
  real beta1;
}


transformed parameters {
  vector[n] prob;

  for (i in 1:n) {
    real lambda = exp(loglambda);
    real c = exp(logc);
    vector[n] eta = beta0 + x*beta1;
    
    if(eta[i] > 0){
      prob[i] = pow(1 - 0.5 * exp(- (eta[i]) ^c), lambda);
    }else{
      prob[i] = pow(0.5*exp(- fabs(eta[i])^c), lambda);
    }
  }
}




// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  loglambda ~ uniform(-2,2);
  logc ~ uniform(-2,2);
  beta0 ~ normal(0.0,100);
  beta1 ~ normal(0.0,100);
  y ~ bernoulli(prob);
}
