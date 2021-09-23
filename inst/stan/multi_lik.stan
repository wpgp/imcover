// imcover - Bayesian model of immunisation coverage
// multi-likelihood, shared latent mean model
// country + vacc + country x time + vacc x time + country x vacc x time

data {
  int<lower=0> N;  // number of obs
  int<lower=0> N_a;  // number of admin records
  int<lower=0> N_o;  // number of official records
  int<lower=0> N_s;  // number of survey records

  int<lower=0> N_i;  // number of countries
  int<lower=0> N_j;  // number of vaccines
  int<lower=0> N_t;  // number of timepoints

  int<lower=1, upper=N_i> i[N];  // country
  int<lower=1, upper=N_j> j[N];  // vaccine
  int<lower=1, upper=N_t> t[N];  // year

  int<lower=1> start_o;
  int<lower=1> start_s;

  real y[N];  // logit-coverage (0, 1), corrected for >100%
}


parameters {
  real lambda_a;  // source-specific intercepts
  real lambda_o;
  real lambda_s;

  real beta_i[N_i];  // i-th country random effect
  real alpha_j[N_j];  // j-th vaccine random effect
  real gamma_t [N_t];  // t-th time effect

  // interactions
  real phi_it[N_i, N_t];  // t-th time point, replicated per country
  real delta_jt[N_j, N_t];  // t-th time point, replicated per vaccine
  real psi_ij[N_i, N_j];  // country-vaccine random effect
  matrix[N_j, N_t] omega_ijt[N_i];

  // AR(1) correlations
  real<lower=-1, upper=1> rho_i;
  real<lower=-1, upper=1> rho_j;
  real<lower=-1, upper=1> rho_t;
  real<lower=-1, upper=1> rho_ij;

  // standard deviation
  real<lower=0> sigma_a; // source-specific SD
  real<lower=0> sigma_o;
  real<lower=0, upper=0.4> sigma_s;

  real<lower=0> sigma_i;
  real<lower=0> sigma_j;
  real<lower=0> sigma_t;
  real<lower=0> sigma_it;
  real<lower=0> sigma_jt;
  real<lower=0> sigma_ij;
  real<lower=0> sigma_ijt;
}


transformed parameters {
 matrix[N_t, N_j] mu[N_i];

  // shared mean
  for(ii in 1:N_i){  // each country
    for(jj in 1:N_j){  // each vaccine
      for(tt in 1:N_t){  // each time
        mu[ii, tt, jj] = beta_i[ii] + alpha_j[jj] + gamma_t[tt] + phi_it[ii, tt] + delta_jt[jj, tt] + psi_ij[ii, jj] + omega_ijt[ii, jj, tt];
      }
    }
  }
}


model {
  beta_i ~ normal(0, sigma_i);
  alpha_j ~ normal(0, sigma_j);

  lambda_a ~ normal(0, 1);
  lambda_o ~ normal(0, 1);
  lambda_s ~ normal(0, 1);

  sigma_a ~ cauchy(0, 2);
  sigma_o ~ cauchy(0, 2);
  sigma_s ~ cauchy(0, 0.2);

  // Autoregressive models (AR(1))
  // time
  gamma_t[1] ~ normal(0, sqrt(pow(sigma_t, 2) / (1-pow(rho_t, 2))));

  for(tt in 2:N_t){
	  gamma_t[tt] ~ normal(rho_t * gamma_t[tt-1], sigma_t);
  }

  // country - time
  for(ii in 1:N_i){
    phi_it[ii, 1] ~ normal(0, sqrt(pow(sigma_it, 2) / (1-pow(rho_i, 2))));

    for(tt in 2:N_t){
      phi_it[ii, tt] ~ normal(rho_i * phi_it[ii, tt-1], sigma_it);
    }
  }

  // vaccine - time
  for(jj in 1:N_j){
    delta_jt[jj, 1] ~ normal(0, sqrt(pow(sigma_jt, 2) / (1-pow(rho_j, 2))));

    for(tt in 2:N_t){
      delta_jt[jj, tt] ~ normal(rho_j * delta_jt[jj, tt-1], sigma_jt);
    }
  }

  // country-vaccine
  for(ii in 1:N_i){
  	for(jj in 1:N_j){
  		psi_ij[ii, jj] ~ normal(0, sigma_ij);
  	}
  }

  // country - vaccine - time
  for(ii in 1:N_i){
    for(jj in 1:N_j){
      omega_ijt[ii, jj, 1] ~ normal(0, sqrt(pow(sigma_ijt, 2) / (1-pow(rho_ij, 2))));

      for(tt in 2:N_t){
        omega_ijt[ii, jj, tt] ~ normal(rho_ij * omega_ijt[ii, jj, tt-1], sigma_ijt);
      }
    }
  }

  // likelihoods
  for(n_a in 1:N_a)
    y[n_a] ~ normal(lambda_a + mu[i[n_a], t[n_a], j[n_a]], sigma_a);

  for(n_o in start_o:(start_o + N_o - 1))
    y[n_o] ~ normal(lambda_o + mu[i[n_o], t[n_o], j[n_o]], sigma_o);

  for(n_s in start_s:(start_s + N_s - 1))
    y[n_s] ~ normal(lambda_s + mu[i[n_s], t[n_s], j[n_s]], sigma_s);
}


generated quantities{
  vector[N] log_lik;

  for(n_a in 1:N_a)
    log_lik[n_a] = normal_lpdf(y[n_a] | lambda_a + mu[i[n_a], t[n_a], j[n_a]], sigma_a);

  for(n_o in start_o:(start_o + N_o - 1))
    log_lik[n_o] = normal_lpdf(y[n_o] | lambda_o + mu[i[n_o], t[n_o], j[n_o]], sigma_o);

  for(n_s in start_s:(start_s + N_s - 1))
    log_lik[n_s] = normal_lpdf(y[n_s] | lambda_s + mu[i[n_s], t[n_s], j[n_s]], sigma_s);
}

