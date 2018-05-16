data {
  int<lower=1> N; // num subj
  int<lower=1> T; // num trials per subj
  int<lower=0, upper=1> group_n_cond[N]; // indicator nocebo group
  int<lower=0, upper=1> group_p_cond[N]; // indicator placebo group
  int<lower=0, upper=1> group_n_cntrl[N]; // indicator control+negative task manipulation group
  int<lower=0, upper=1> group_p_cntrl[N]; // indicator control+positive task manipulation group

  int<lower=1, upper=3> pair[N,T]; // easy, medium, hard
  int<lower=0, upper=1> day2[N,T];  // indicator for day2
  int<lower=0, upper=1> outcome[N, T]; // reward?
  int<lower=0, upper=1> accuracy[N, T]; // correct object chosen?
}
transformed data {
  vector[6] Q0;
  Q0  = rep_vector(0.0, 6);
}
parameters {
  // declare as vectors for vectorizing
  vector[3] mu_p;  
  vector<lower=0>[3] sigma; 
  real<lower=0> sigma_b; // constrain bs
  vector[N] alpha_g_pr;          // learning rate gains
  vector[N] alpha_l_pr;          // learning rate losses
  vector[N] beta_pr;          // inv temp
  real b_alg_day2; // fixed effects for day2
  real b_alg_pcond; // group effects for day2
  real b_alg_ncond; 
  real b_alg_pctrl; 
  real b_alg_nctrl; 
  
  real b_all_day2; // fixed effects for day2
  real b_all_pcond; // group effects for day
  real b_all_ncond; 
  real b_all_pctrl; 
  real b_all_nctrl; 

  real b_be_day2;
  real b_be_pcond; // group effects for day2
  real b_be_ncond; 
  real b_be_pctrl; 
  real b_be_nctrl; 
  
}
transformed parameters{
  vector<lower=0,upper=1>[N] alpha_g;
  vector<lower=0,upper=1>[N] alpha_l;
  vector<lower=0>[N] beta;
  vector<lower=0,upper=1>[N] alpha_g_day2;
  vector<lower=0,upper=1>[N] alpha_l_day2;
  vector<lower=0>[N] beta_day2;
  
  for (i in 1:N) {
    alpha_g[i]  = Phi_approx( mu_p[1] + sigma[1] * alpha_g_pr[i] );
    alpha_g_day2[i] = Phi_approx( mu_p[1] + // grand mean
                      b_alg_day2 +          // intercept day2
                      b_alg_pcond*group_p_cond[i]  + b_alg_ncond*group_n_cond[i] +  // group x day IAs
                      b_alg_pctrl*group_p_cntrl[i] + b_alg_nctrl*group_n_cntrl[i] + 
                      sigma[1] * alpha_g_pr[i] ); // random intercept
    
    alpha_l[i]  = Phi_approx( mu_p[2] + sigma[2] * alpha_l_pr[i] );
    alpha_l_day2[i] = Phi_approx( mu_p[2] + // grand mean
                      b_all_day2 +          // intercept day2
                      b_all_pcond*group_p_cond[i]  + b_all_ncond*group_n_cond[i] +  // group x day IAs
                      b_all_pctrl*group_p_cntrl[i] + b_all_nctrl*group_n_cntrl[i] + 
                      sigma[2] * alpha_l_pr[i] ); // random intercept

    beta[i]  = exp( mu_p[3] + sigma[3] * beta_pr[i] );
    beta_day2[i]  = exp( mu_p[3] + 
                    b_be_day2 +          // intercept day2
                    b_be_pcond*group_p_cond[i]  + b_be_ncond*group_n_cond[i] +  // group x day IAs
                    b_be_pctrl*group_p_cntrl[i] + b_be_nctrl*group_n_cntrl[i] + 
                    sigma[3] * beta_pr[i] );
  }
}

model {  
  // weakly informative prior
  mu_p[1] ~ normal(-0.5, 0.6);
  sigma[1] ~ cauchy(0, 0.01);

  mu_p[2] ~ normal(-0.5, 0.6);
  sigma[2] ~ cauchy(0, 0.01);
  
  mu_p[3] ~ normal(-1.5, 0.8);
  sigma[3] ~ normal(0, 0.3);
  
  // individual parameters w/ Matt trick see https://groups.google.com/forum/#!msg/stan-users/4gv3fNCqSNk/J6ZItL2ZJ-IJ
  alpha_g_pr  ~ normal(0, 1.0);   
  alpha_l_pr  ~ normal(0, 1.0);   
  beta_pr  ~ normal(0, 1.0);   
  
  sigma_b ~ normal(0,1);
  b_alg_day2 ~ normal(0,sigma_b);
  b_all_day2 ~ normal(0,sigma_b);
  b_be_day2 ~ normal(0,sigma_b);
  b_alg_pcond ~ normal(0,sigma_b); // group effects for day2
  b_alg_ncond ~ normal(0,sigma_b); 
  b_alg_pctrl ~ normal(0,sigma_b); 
  b_alg_nctrl ~ normal(0,sigma_b); 
  
  b_all_pcond ~ normal(0,sigma_b); // group effects for day2
  b_all_ncond ~ normal(0,sigma_b); 
  b_all_pctrl ~ normal(0,sigma_b); 
  b_all_nctrl ~ normal(0,sigma_b); 

  b_be_pcond ~ normal(0,sigma_b); // group effects for day2
  b_be_ncond ~ normal(0,sigma_b); 
  b_be_pctrl ~ normal(0,sigma_b); 
  b_be_nctrl ~ normal(0,sigma_b); 

  for (i in 1:N) {
    
    vector[6] Q;  // Q value for all actions
    real pcorrect;   // prob of correct choice
    int choice; // index of chosen action
    real pe; // prediction error
    Q = Q0;

    for (t in 1:T)  {
      if(t==241){
        Q=Q0; // reset before start of day2
      }
      pcorrect = inv_logit( (Q[2*pair[i,t]-1]-Q[2*pair[i,t]])/( beta[i]*(1-day2[i,t])+beta_day2[i]*day2[i,t]) );
      accuracy[i,t] ~ bernoulli( pcorrect );
      
      // update action values
      choice=2*pair[i,t];
      if(accuracy[i,t]==1){
        choice=choice-1;
      } 
      
      if(outcome[i,t]>0){ // gain
        Q[choice] = Q[choice] + (alpha_g[i]*(1-day2[i,t])+alpha_g_day2[i]*day2[i,t]) * (outcome[i,t]-Q[choice]);
      } else { // loss
        Q[choice] = Q[choice] + (alpha_l[i]*(1-day2[i,t])+alpha_l_day2[i]*day2[i,t]) * (outcome[i,t]-Q[choice]);
      }
    } // end of t loop (T trials)
  } // end of i loop (N subjects)
}


generated quantities {
  real<lower=0, upper=1> mu_alpha_g;
  real<lower=0, upper=1> mu_alpha_l;  
  real<lower=0> mu_beta;
  real log_lik[N];
  
  mu_alpha_g  = Phi_approx(mu_p[1]);
  mu_alpha_l  = Phi_approx(mu_p[2]);
  mu_beta  = exp(mu_p[3]);

  { // local section, this saves time and space
    for (i in 1:N) {
      vector[6] Q;  // Q value for all actions
      real pcorrect;   // prob of correct choice
      int choice; // index of chosen action
      
      log_lik[i] = 0;

      Q = Q0;
  
      for (t in 1:T)  {
        if(t==241){
          Q=Q0; // reset before start of day2
        }
        
        pcorrect = inv_logit( (Q[2*pair[i,t]-1]-Q[2*pair[i,t]])/( beta[i]*(1-day2[i,t])+beta_day2[i]*day2[i,t]) );
        log_lik[i] = log_lik[i] + bernoulli_lpmf( accuracy[i,t] | pcorrect );

        // update action values
        choice=2*pair[i,t];
        if(accuracy[i,t]==1){
          choice=choice-1;
        } 
        if(outcome[i,t]>0){ // gain
          Q[choice] = Q[choice] + (alpha_g[i]*(1-day2[i,t])+alpha_g_day2[i]*day2[i,t]) * (outcome[i,t]-Q[choice]);
        } else { // loss
          Q[choice] = Q[choice] + (alpha_l[i]*(1-day2[i,t])+alpha_l_day2[i]*day2[i,t]) * (outcome[i,t]-Q[choice]);
        }

      } // end of t loop (T trials)
    } // end of i loop (N subjects)
  } // end of local section
}  
