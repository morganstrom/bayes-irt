data{
  int<lower=2, upper=5> K;                  //number of categories
  int<lower=0> N;                           //number of students
  int<lower=0> I;                           //number of items
  int<lower=1,upper=K> X[N, I];             //data matrix
}
parameters {
  vector[N] theta;                          //person ability
  real<lower=0> alpha[I];                   //item discrimination
  ordered[K-1] kappa[I];                    //category difficulty
  real mu_kappa;                            //mean of the prior distribution of category difficulty
  real<lower=0> sigma_kappa;                //sd of the prior distribution of category difficulty
}
model {
  alpha ~ lognormal(1,1);                   //item discrimination prior
  theta ~ normal(0,1);                      //person ability prior
  mu_kappa ~ normal(0,5);                   //kappa hyperprior mean
  sigma_kappa ~ cauchy(0,5);                //kappa hyperprior sd
  
  //category diffuculties, one per category boundary per item
  for (i in 1:I){
    for (k in 1:(K-1)){
      kappa[i,k] ~ normal(mu_kappa, sigma_kappa); 
    }
  }
  
  //model each rating as ordered logistic
  for (n in 1:N){
    for (i in 1:I){
      X[n,i] ~ ordered_logistic(theta[n]*alpha[i], kappa[i]); 
    }
  }
}
