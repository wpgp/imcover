// imcover - Bayesian model of immunisation coverage
// multi-likelihood, shared latent mean model
// country + vacc + country x time + vacc x time + country x vacc x time

data {
  int<lower=1> N;  // number of obs

  int<lower=1> N_i;  // number of countries
  int<lower=1> N_j;  // number of vaccines
  int<lower=1> N_t;  // number of timepoints

  int<lower=1, upper=N_i> i[N];  // country
  int<lower=1, upper=N_j> j[N];  // vaccine
  int<lower=1, upper=N_t> t[N];  // year

  vector[N] y;  // logit-coverage (0, 1), corrected for >100%

  int<lower=1> nsources;
  int<lower=1, upper=nsources> source[N];
  int<lower=1> sizes[nsources];

  // vector[nsources] U_lambda;  // upper-bounds
  vector[nsources] U_sigma;
  // vector[nsources] L_lambda;  // lower-bounds
  vector[nsources] L_sigma;

  real prior_lambda[nsources]; // priors
  real prior_sigma[nsources];

  // shared mean identifiers
  int<lower=1, upper=N_i * N_j * N_t> mu_lookup[N];  // which record?
  int<lower=1, upper=N_i> ii[N_i * N_j * N_t];
  int<lower=1, upper=N_j> jj[N_i * N_j * N_t];
  int<lower=1, upper=N_t> tt[N_i * N_j * N_t];
}


parameters {
  // regression effects
  vector[nsources] lambda; // intercepts
  real beta_i[N_i == 1 ? 0 : N_i];  // i-th country random effect
  real alpha_j[N_j];  // j-th vaccine random effect
  real gamma_t [N_t];  // t-th time effect

  // interactions
  real phi_it[N_i == 1 ? 0 : N_i, N_t];  // t-th time point, replicated per country
  real delta_jt[N_j, N_t];  // t-th time point, replicated per vaccine
  real psi_ij[N_i == 1 ? 0 : N_i, N_j];  // country-vaccine random effect
  matrix[N_i == 1 ? 0 : N_j, N_t] omega_ijt[N_i == 1 ? 0 : N_i];

  // AR(1) correlations
  real<lower=-1, upper=1> rho_i[N_i == 1 ? 0 : 1];
  real<lower=-1, upper=1> rho_j;
  real<lower=-1, upper=1> rho_t;
  real<lower=-1, upper=1> rho_ij[N_i == 1 ? 0 : 1];

  // standard deviation
  vector<lower=0, upper=1>[nsources] sigma_raw;  // source-specific

  real<lower=0> sigma_i[N_i == 1 ? 0 : 1];
  real<lower=0> sigma_j;
  real<lower=0> sigma_t;
  real<lower=0> sigma_it[N_i == 1 ? 0 : 1];
  real<lower=0> sigma_jt;
  real<lower=0> sigma_ij[N_i == 1 ? 0 : 1];
  real<lower=0> sigma_ijt[N_i == 1 ? 0 : 1];
}


transformed parameters {
  // shared mean
  vector[N_t * N_j * N_i] mu;

  // varying bounds;
  vector[nsources] sigma = L_sigma + (U_sigma - L_sigma) .* sigma_raw;

  // shared mean
  if(N_i > 1){
    for(idx in 1:num_elements(mu)){
      mu[idx] = beta_i[ii[idx]] + alpha_j[jj[idx]] + gamma_t[tt[idx]] + phi_it[ii[idx], tt[idx]] + delta_jt[jj[idx], tt[idx]] + psi_ij[ii[idx], jj[idx]] + omega_ijt[ii[idx], jj[idx], tt[idx]];
    }
  } else{
    for(idx in 1:num_elements(mu)){
      mu[idx] = alpha_j[jj[idx]] + gamma_t[tt[idx]] + delta_jt[jj[idx], tt[idx]];
    }
  }
}


model {
  alpha_j ~ normal(0, sigma_j);

  for(s in 1:nsources){
    segment(lambda, s, 1) ~ normal(0, prior_lambda[s]);
    segment(sigma_raw, s, 1) ~ cauchy(0, prior_sigma[s]);
  }

  // Autoregressive models (AR(1))
  // time
  gamma_t[1] ~ normal(0, sqrt(pow(sigma_t, 2) / (1-pow(rho_t, 2))));

  for(time in 2:N_t){
	  gamma_t[time] ~ normal(rho_t * gamma_t[time-1], sigma_t);
  }

  // vaccine - time
  for(vax in 1:N_j){
    delta_jt[vax, 1] ~ normal(0, sqrt(pow(sigma_jt, 2) / (1-pow(rho_j, 2))));

    for(time in 2:N_t){
      delta_jt[vax, time] ~ normal(rho_j * delta_jt[vax, time-1], sigma_jt);
    }
  }

  // for multi-country (regional) models
  if(N_i > 1){
    beta_i ~ normal(0, sigma_i[1]);

    // country - time
    for(ctry in 1:N_i){
      phi_it[ctry, 1] ~ normal(0, sqrt(pow(sigma_it[1], 2) / (1-pow(rho_i[1], 2))));

      for(time in 2:N_t){
        phi_it[ctry, time] ~ normal(rho_i[1] * phi_it[ctry, time-1], sigma_it[1]);
      }
    }

    // country - vaccine
    for(ctry in 1:N_i){
    	for(vax in 1:N_j){
    		psi_ij[ctry, vax] ~ normal(0, sigma_ij[1]);
    	}
    }

    // country - vaccine - time
    for(ctry in 1:N_i){
      for(vax in 1:N_j){
        omega_ijt[ctry, vax, 1] ~ normal(0, sqrt(pow(sigma_ijt[1], 2) / (1-pow(rho_ij[1], 2))));

        for(time in 2:N_t){
          omega_ijt[ctry, vax, time] ~ normal(rho_ij[1] * omega_ijt[ctry, vax, time-1], sigma_ijt[1]);
        }
      }
    }
  }

  // likelihoods
  {
    int start = 1;
    // loop over the different data sources
    for(idx in 1:nsources){
      int end = start + sizes[idx] - 1;
      y[start:end] ~ normal(lambda[source[start:end]] + mu[mu_lookup[start:end]], sigma[source[start:end]]);

      start = end + 1;
    }
  }

  target += sum(sigma_raw);
}


generated quantities{
  vector[N] log_lik;

  for(idx in 1:N){
    log_lik[idx] = normal_lpdf(y[idx] | lambda[source[idx]] + mu[mu_lookup[idx]], sigma[source[idx]]);
  }
}

