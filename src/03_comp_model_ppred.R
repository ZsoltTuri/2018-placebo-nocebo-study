##================================================================================================
# posterior predictive

get.ppred <- function(mod, nrep, stan.data){
  mu_p=extract(mod, "mu_p")[[1]]
  sigma=extract(mod, "sigma")[[1]]
  nsamp=dim(mu_p)[1]
  alpha_gs = extract(mod, "alpha_g")[[1]]
  alpha_ls = extract(mod, "alpha_l")[[1]]
  betas  = extract(mod, "beta")[[1]]
  alpha_g_day2s=extract(mod, "alpha_g_day2")[[1]]
  alpha_l_day2s=extract(mod, "alpha_l_day2")[[1]]
  beta_day2s =extract(mod, "beta_day2")[[1]]
  
  contingencies=matrix(c(0.8,0.2,
                         0.7,0.3,
                         0.6,0.4), nrow=3,byrow=T)
  
  attach(stan.data)
  mat=cbind(group_p_cntrl,group_n_cntrl, group_p_cond, group_n_cond)
  groupfac=factor(mat%*%(1:4), labels = c("nhg",colnames(mat)) )
  
  drep=data.frame(
    nrep=0, ## nrep==0 is real data
    subj=rep(1:N, each=T),
    trial=rep(1:(T/2),2*N),
    group=rep(groupfac,each=T),
    pair=as.vector(t(pair)),
    day=as.vector(t(stan.data$day)),
    outcome=as.vector(t(outcome)),
    accuracy=as.vector(t(accuracy))
  )

  for( k in 1:nrep ){
    cat(".")
    acc=matrix(0, nrow=N, ncol=T)
    rew=matrix(0, nrow=N, ncol=T)
    ix=sample(1:nsamp, 1)
    alpha_g=alpha_gs[ix,]
    alpha_l=alpha_ls[ix,]
    beta =betas[ix,]
    alpha_g_day2=alpha_g_day2s[ix,]
    alpha_l_day2=alpha_l_day2s[ix,]
    beta_day2 =beta_day2s[ix,]
    
    
    for (i in 1:N) {
      Q=rep(0,6);

      for (t in 1:T)  {
        if(t==241){
          Q=rep(0,6); # reset before start of day2
        }
        pcorrect = boot::inv.logit( (Q[2*pair[i,t]-1]-Q[2*pair[i,t]])/( beta[i]*(1-day2[i,t])+beta_day2[i]*day2[i,t]) );
        acc[i,t] = sample(c(0,1), size=1, replace = T, prob = c( 1-pcorrect, pcorrect ) )
        
        # update action values
        choice=2*pair[i,t];
        if(acc[i,t]==1){
          choice=choice-1;
        } 
        if(acc[i,t]==1)
          rew[i,t]=sample( c(1,0), size=1, prob=contingencies[pair[i,t],])
        else
          rew[i,t]=sample( c(0,1), size=1, prob=contingencies[pair[i,t],])
        
        if(rew[i,t]>0){
          Q[choice] = Q[choice] + (alpha_g[i]*(1-day2[i,t])+alpha_g_day2[i]*day2[i,t]) * (rew[i,t]-Q[choice]);
        } else {
          Q[choice] = Q[choice] + (alpha_l[i]*(1-day2[i,t])+alpha_l_day2[i]*day2[i,t]) * (rew[i,t]-Q[choice]);
        }
      } # end of t loop (T trials)
    } # end of i loop (N subjects)
  
    dd<-data.frame(
      nrep=k,
      subj=rep(1:N, each=T),
      trial=rep(1:(T/2),2*N),
      group=rep(groupfac,each=T),
      pair=as.vector(t(pair)),
      day=as.vector(t(stan.data$day)),
      outcome=as.vector(t(rew)),
      accuracy=as.vector(t(acc))
    )
    drep = rbind(drep,dd)
  }
  detach("stan.data")
  
  return(drep)
}

## return probability for each given response under the model for many draws from the posterior distribution
get.presp <- function(mod, nrep, stan.data){
  nsamp=dim(as.matrix(mod))[1]
  alpha_gs = extract(mod, "alpha_g")[[1]]
  alpha_ls = extract(mod, "alpha_l")[[1]]
  betas  = extract(mod, "beta")[[1]]
  alpha_g_day2s=extract(mod, "alpha_g_day2")[[1]]
  alpha_l_day2s=extract(mod, "alpha_l_day2")[[1]]
  beta_day2s =extract(mod, "beta_day2")[[1]]
  
  attach(stan.data)
  mat=cbind(group_p_cntrl,group_n_cntrl, group_p_cond, group_n_cond)
  groupfac=factor(mat%*%(1:4), labels = c("nhg",colnames(mat)) )
  
  drep=NULL
  
  for( k in 1:nrep ){
    cat(".")
    pcorrect=matrix(0, nrow=N, ncol=T)
    presponse=matrix(0, nrow=N, ncol=T)
    ix=sample(1:nsamp, 1)
    alpha_g=alpha_gs[ix,]
    alpha_l=alpha_ls[ix,]
    beta =betas[ix,]
    alpha_g_day2=alpha_g_day2s[ix,]
    alpha_l_day2=alpha_l_day2s[ix,]
    beta_day2 =beta_day2s[ix,]
    
    
    for (i in 1:N) {
      Q=rep(0,6);
      
      for (t in 1:T)  {
        if(t==241){
          Q=rep(0,6); # reset before start of day2
        }
        pcorrect[i,t] = boot::inv.logit( (Q[2*pair[i,t]-1]-Q[2*pair[i,t]])/( beta[i]*(1-day2[i,t])+beta_day2[i]*day2[i,t]) );
        presponse[i,t] = ifelse(accuracy[i,t]>0, pcorrect[i,t], 1-pcorrect[i,t])

        # update action values
        choice=2*pair[i,t];
        if(accuracy[i,t]==1){
          choice=choice-1;
        } 
        if(outcome[i,t]>0){
          Q[choice] = Q[choice] + (alpha_g[i]*(1-day2[i,t])+alpha_g_day2[i]*day2[i,t]) * (outcome[i,t]-Q[choice]);
        } else {
          Q[choice] = Q[choice] + (alpha_l[i]*(1-day2[i,t])+alpha_l_day2[i]*day2[i,t]) * (outcome[i,t]-Q[choice]);
        }
      } # end of t loop (T trials)
    } # end of i loop (N subjects)
    
    dd<-data.frame(
      nrep=k,
      subj=rep(1:N, each=T),
      trial=rep(1:(T/2),2*N),
      group=rep(groupfac,each=T),
      pair=as.vector(t(pair)),
      day=as.vector(t(stan.data$day)),
      outcome=as.vector(t(outcome)),
      accuracy=as.vector(t(accuracy)),
      pcorrect=as.vector(t(pcorrect)),
      presponse=as.vector(t(presponse))
    )
    drep = rbind(drep,dd)
  }
  detach("stan.data")
  
  return(drep)
}