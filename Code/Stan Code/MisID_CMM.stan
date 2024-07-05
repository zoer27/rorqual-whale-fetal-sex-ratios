//Model for fetal sex mis-ID
//this one is specific for common minke whales
//limit L50 and delta to be small

data {
  int<lower =0> N; //number of observations
  array[N] int<lower = 0> DatBin; //binomial data of fetal sex 
  vector<lower = 0>[N] fetalLengths; //fetal lengths data
  int<lower = 0> correction; //is there a correction switch
  int<lower =0> correctionMale; //switch to turn on male correction
  int<lower = 0> correctionFemale; //switch to turn on female correction
  real maxL; //maximum length for correction
}



parameters {
 //linear model
    real a10; //intercept (sex ratio at 0 length)
    real b; //slope
    array[correction] real<upper = 2> L50; //if correction is 0 this has no dimension
    array[correction] real<lower = 0, upper = 2> delta;
}

transformed parameters{
  vector<lower = 0, upper = 1>[N*correction] pcor; //probability correct
  vector<lower = 0, upper = 1>[N] predm; //predicted proportion male with correction
  vector<lower = 0, upper = 1>[N] pmale; //linear proportion male on logit scale
  
  //linear mod
  pmale = inv_logit(a10 + b*fetalLengths); //vetorized
  
  //probability correct
  if(correction == 1)
    for(l in 1:N){
      if(fetalLengths[l] <= maxL)
        pcor[l] = (1.0 + exp((-log(19.0)*(fetalLengths[l] - L50[1]))/delta[1]))^-1.0;
      
      else
        pcor[l] = 1.0;
        
    }
  
  
  //prediction
    if(correctionMale == 1){
      predm = pmale .* pcor; //vectorized
    }
    else if(correctionFemale == 1){
      predm = pmale + (1-pmale).*(1-pcor); //vectorized
    }
    else{
      predm = pmale;
    }
}

model {
  //priors
    a10 ~ normal(0, 10); //logit transformation
    b ~ normal(0, 10);
  
  if(correction == 1){
    L50 ~ normal(0, 10); //truncated because of lower limit above
    delta ~ normal(0, 10); //truncated because of lower limit above
  }
  
  //likelihood
  DatBin ~ bernoulli(predm);
  
}

generated quantities{
  array[N] real loglike; //loglikelihood for WAIC
  array[correction] real L99; //estimated L99 
  for(n in 1:N){
    loglike[n] = bernoulli_lpmf(DatBin[n] | predm[n]);
  }
  
  if(correction == 1){
     L99[correction] = ((delta[1]*log((1/0.99) - 1))/-log(19)) + L50[1];
  }
  
}


