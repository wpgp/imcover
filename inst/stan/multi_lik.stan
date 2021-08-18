// Bayesian model of vaccination coverage
// country + vacc + country x time + vacc x time

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
  real lambda_a[N_j];  // vaccine-source-specific intercepts
  real lambda_o[N_j];
  real lambda_s[N_j];

  real beta_i[N_i];  // i-th country random effect
  real alpha_j[N_j];  // j-th vaccine random effect
  real phi_t [N_t]; //t-th time effect
  real gamma_it[N_i, N_t];  // t-th time point, replicated per country
  real gamma_jt[N_j, N_t];  // t-th time point, replicated per vaccine
  real gamma_ij[N_i, N_j]; //country-vaccine random effect
  matrix[N_j, N_t] delta_ijt[N_i];

  real<lower=-1, upper=1> rho_i; // AR(1) corr
  real<lower=-1, upper=1> rho_j;
  real<lower=-1, upper=1> rho_t;
  real<lower=-1, upper=1> rho_ij;

  // standard deviation
  real<lower=0> sigma_a; // source-specific intercepts
  real<lower=0> sigma_o;
  //real<lower=0> sigma_s;
  real<lower=0, upper=0.4> sigma_s;  //Note this

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
  for(ii in 1:N_i){
    for(jj in 1:N_j){
      for(tt in 1:N_t){
        mu[ii, tt, jj] = beta_i[ii] + alpha_j[jj] + phi_t[tt] + gamma_it[ii, tt] + gamma_jt[jj, tt] + gamma_ij[ii, jj] + delta_ijt[ii, jj, tt];
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

  sigma_a ~ cauchy(0, 2);  //Note that truncated cauchy or normal priors worked well for sigma_s          0.4
  sigma_o ~ cauchy(0, 2);  //For sigma_o and sigma_a, half-normal or half-cauchy priors also worked well  0.4
  sigma_s ~ cauchy(0, 0.2); //0.002

  //Autoregressive model (AR(1))
  phi_t[1] ~ normal(0, sqrt(pow(sigma_t, 2) / (1-pow(rho_t, 2))));
  for(tt in 2:N_t){
	phi_t[tt] ~ normal(rho_t * phi_t[tt-1], sigma_t);
  }

  // country - time
  for(ii in 1:N_i){ //N_i countries
    gamma_it[ii, 1] ~ normal(0, sqrt(pow(sigma_it, 2) / (1-pow(rho_i, 2))));

    for(tt in 2:N_t){
      gamma_it[ii, tt] ~ normal(rho_i * gamma_it[ii, tt-1], sigma_it);
    }
  }

  // vaccine - time
  for(jj in 1:N_j){
    gamma_jt[jj, 1] ~ normal(0, sqrt(pow(sigma_jt, 2) / (1-pow(rho_j, 2))));

    for(tt in 2:N_t){
      gamma_jt[jj, tt] ~ normal(rho_j * gamma_jt[jj, tt-1], sigma_jt);
    }
  }

  //country-vaccine
  for(ii in 1:N_i){
  	for(jj in 1:N_j){
  		gamma_ij[ii, jj] ~ normal(0, sigma_ij);
  	}
  }

  // country - vaccine - time
  for(ii in 1:N_i){
    for(jj in 1:N_j){
      delta_ijt[ii, jj, 1] ~ normal(0, sqrt(pow(sigma_ijt, 2) / (1-pow(rho_ij, 2))));

      for(tt in 2:N_t){
        delta_ijt[ii, jj, tt] ~ normal(rho_ij * delta_ijt[ii, jj, tt-1], sigma_ijt);
      }
    }
  }


  // likelihoods
  for(n_a in 1:N_a)
    y[n_a] ~ normal(mu[i[n_a], t[n_a], j[n_a]], sigma_a);

  for(n_o in start_o:(start_o + N_o - 1))
    y[n_o] ~ normal(mu[i[n_o], t[n_o], j[n_o]], sigma_o);

  for(n_s in start_s:(start_s + N_s - 1))
    y[n_s] ~ normal(mu[i[n_s], t[n_s], j[n_s]], sigma_s);
}


generated quantities{
  vector[N] log_lik;

  for(n_a in 1:N_a)
    log_lik[n_a] = normal_lpdf(y[n_a] | mu[i[n_a], t[n_a], j[n_a]], sigma_a);

  for(n_o in start_o:(start_o + N_o - 1))
    log_lik[n_o] = normal_lpdf(y[n_o] | mu[i[n_o], t[n_o], j[n_o]], sigma_o);

  for(n_s in start_s:(start_s + N_s - 1))
    log_lik[n_s] = normal_lpdf(y[n_s] | mu[i[n_s], t[n_s], j[n_s]], sigma_s);
}

