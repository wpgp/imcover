// Bayesian model of vaccination coverage
// Source-specific RE like INLA
// country + vacc + country x time + vacc x time

data {
  int<lower=0> N;  // number of obs
  
  int<lower=0> N_i;  // number of countries
  int<lower=0> N_j;  // number of vaccines
  int<lower=0> N_t;  // number of timepoints
  int<lower=0> N_s;  // number of sources
  
  int<lower=1, upper=N_i> i[N];  // country 
  int<lower=1, upper=N_j> j[N];  // vaccine 
  int<lower=1, upper=N_t> t[N];  // year
  int<lower=1, upper=N_s> s[N];  // source
  
  real y[N];  // logit-coverage (0, 1), corrected for >100%
}


parameters {
  real lambda; //overall intercept
  real beta_i[N_i];  // i-th country random effect
  real alpha_j[N_j];  // j-th vaccine random effect
  real phi_t[N_t];  // t-th time random effect
  real psi_s[N_s]; // s-th source random effect
  real gamma_it[N_i, N_t];  // t-th time point, replicated per country
  real gamma_ij[N_i, N_j];  // country-vaccine effect
  real gamma_jt[N_j, N_t];  // t-th time point, replicated per vaccine
  matrix[N_j, N_t] delta_ijt[N_i];
  
  real<lower=-1, upper=1> rho_i; // AR(1) corr
  real<lower=-1, upper=1> rho_j;
  real<lower=-1, upper=1> rho_t;
  real<lower=-1, upper=1> rho_ij;
  
  // standard deviation
  real<lower=0> sigma;
  real<lower=0> sigma_lam;
  real<lower=0> sigma_s;
  real<lower=0> sigma_s3;
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
        mu[ii, tt, jj] = lambda + beta_i[ii] + alpha_j[jj] + phi_t[tt] + gamma_it[ii, tt] + gamma_jt[jj, tt] + gamma_ij[ii, jj] + delta_ijt[ii, jj, tt];
      }
    }
  }
}


model {
  beta_i ~ normal(0, sigma_i);
  alpha_j ~ normal(0, sigma_j);
  lambda ~ normal(0, sigma_lam);
  psi_s[1] ~ normal(0, sigma_s); //admin
  psi_s[2] ~ normal(0, sigma_s); //official
  psi_s[3] ~ normal(0, sigma_s3); //survey
  
  sigma ~ cauchy(0, 2);
  sigma_i ~ cauchy(0, 2);
  sigma_j ~ cauchy(0, 2);
  sigma_s ~ cauchy(0, 2); //10
  sigma_s3 ~ cauchy(0, 2); //0.002
  
  // autoregressive model (AR(1))
  phi_t[1] ~ normal(0, sqrt(pow(sigma_t, 2) / (1-pow(rho_t, 2))));
  for (tt in 2:N_t){
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
  
    // country - vaccine 
  for (ii in 1:N_i){
	for (jj in 1:N_j){
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
  
  
  // likelihood
  for(n in 1:N){
    y[n] ~ normal(psi_s[s[n]] + mu[i[n], t[n], j[n]], sigma);
  }
}


generated quantities{
  vector[N] log_lik;
  
  for(n in 1:N)
    log_lik[n] = normal_lpdf(y[n] | psi_s[s[n]] + mu[i[n], t[n], j[n]], sigma);
}

