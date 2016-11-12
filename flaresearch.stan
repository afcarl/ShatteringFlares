data {
  int nt;
  int nobs;
  int cts[nobs, nt];
}

parameters {
  real<lower=0> bg[nobs];
  real<upper=0> tcent[nobs];
  real<lower=0> amp[nobs];
  real<upper=0> mu_t;
  real<lower=0> sigma_t;
  real<lower=0> mu_a;
  real<lower=0> sigma_a;
}

model {
  // Priors
  bg ~ normal(0.0, 10.0); // Broad prior on bg rates, since we will measure this well
  
  mu_t ~ normal(-5.0, 2.0);  // "Physics" tells us that flares should
			     // be ~five samples +/- 2 before burst
  sigma_t ~ cauchy(2.0, 1.0); // Broadly, scatter should be about 2

  mu_a ~ normal(0.0, 10.0);  // Broad prior on amplitude location.
  sigma_a ~ cauchy(2.0, 1.0); // Broad prior on amplitude scatter

  tcent ~ normal(mu_t, sigma_t);
  amp ~ normal(mu_a, sigma_a);

  // Likelihood
  for (i in 1:nobs) {
    vector[nt] rate;

    for (j in 1:nt) {
      real tlow;
      real thigh;

      tlow = j - (nt + 1);
      thigh = j - nt;

      rate[j] = bg[i];
      rate[j] = rate[j] + 0.5*amp[i]*(erf((tcent[i]-tlow)/sqrt(2.0)) - erf((tcent[i]-thigh)/sqrt(2.0)));
    }

    cts[i] ~ poisson(rate);
  }
}

