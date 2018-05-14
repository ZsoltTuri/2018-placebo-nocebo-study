
sim.and.eval <- function(alpha, beta, ntrials, npairs=3, reward.prob=c(.8, .7, .6)){
  Q=rep(0, npairs*2)
  ACC=0
  for(pair in 1:npairs){
    for(t in 1:ntrials){
      A=2*pair-1
      B=A+1
      pA=exp(Q[A]/beta-log.sum.exp( c(Q[A], Q[B])/beta ))
      a=sample(c(A,B), 1, prob=c(pA, 1-pA))
      reward=sample( c(a==A,a==B), 1, prob=c(reward.prob[pair], 1-reward.prob[pair]))
      Q[a] = Q[a] + alpha*(reward-Q[a])
      ACC = ACC+(a==A)
    }
  }
  return(ACC/(ntrials*npairs))
}

library(Rcpp)
cppFunction('double csimandeval_alphalg(double alphag, double alphal, double beta, int ntrials){
  int npairs=3;
  double reward_prob[3]={.8, .7, .6};
  double ACC=0.0;
  int pair, t;
  double Q[6]={0.0};
  int A,B,a,reward;
  double pA, offset;

  for(pair=0; pair<npairs; pair++){
    for(t=0; t<ntrials; t++){
      A=2*pair;
      B=A+1;
      offset=((Q[A]/beta)>(Q[B]/beta) ? (Q[A]/beta) : (Q[B]/beta));
      pA=exp(Q[A]/beta - (log( exp( Q[A]/beta - offset)+exp( Q[B]/beta - offset)) + offset));
      a=(((double)runif(1)[0])<pA) ? A : B;
      reward=(((double)runif(1)[0])<reward_prob[pair]) ? (a==A) : (a==B);
      //Rprintf("pA=%f, a=%i, reward=%i, ACC=%i\\n", pA, a, reward, a==A);

      if(reward>0){
        Q[a] = Q[a] + alphag*(reward-Q[a]);
      } else {
        Q[a] = Q[a] + alphal*(reward-Q[a]);
      }
      ACC = ACC+(a==A);
    }
  }
  return ACC/(ntrials*npairs);
}
')


cppFunction('double csimandeval(double alpha, double beta, int ntrials){
  int npairs=3;
  double reward_prob[3]={.8, .7, .6};
  double ACC=0.0;
  int pair, t;
  double Q[6]={0.0};
  int A,B,a,reward;
  double pA, offset;

  for(pair=0; pair<npairs; pair++){
    for(t=0; t<ntrials; t++){
      A=2*pair;
      B=A+1;
      offset=((Q[A]/beta)>(Q[B]/beta) ? (Q[A]/beta) : (Q[B]/beta));
      pA=exp(Q[A]/beta - (log( exp( Q[A]/beta - offset)+exp( Q[B]/beta - offset)) + offset));
      a=(((double)runif(1)[0])<pA) ? A : B;
      reward=(((double)runif(1)[0])<reward_prob[pair]) ? (a==A) : (a==B);
      //Rprintf("pA=%f, a=%i, reward=%i, ACC=%i\\n", pA, a, reward, a==A);

      Q[a] = Q[a] + alpha*(reward-Q[a]);
      ACC = ACC+(a==A);
    }
  }
  return ACC/(ntrials*npairs);
}
')


log.sum.exp<- function(x) {
  # Computes log(sum(exp(x))
  # Uses offset trick to avoid numeric overflow: http://jblevins.org/notes/log-sum-exp
  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset
}


## has trouble with numerical overflow in Pt calculation
llqlearning.org <- function( alpha, beta, ACC, pair, reward){
  npairs=max(pair)
  ntrials=length(pair)
  Q=rep(0, 2*npairs)
  
  LL=0
  for( t in 1:ntrials ) {
    if(ACC[t]<0) # catch missing values
      next
    
    # a is selected action in trial t 
    # correct action is always first of pair (A,C,E)
    A=2*pair[t]-1
    B=A+1
    a=A+(ACC[t]<1)
    #cat(sprintf("%i, %i, %i\n", A, B, a))
    Q[a] = Q[a] + alpha*(reward[t]-Q[a])
    Pt=exp(Q[a]/beta)/( exp(Q[A]/beta)+exp(Q[B]/beta) )
    #logPt=Q[a]-log.sum.exp( Q[c(A,B)]/beta )
    LL=LL+log(Pt)
  }
  LL
}

llqlearning <- function( alpha, beta, ACC, pair, reward){
  npairs=max(pair)
  ntrials=length(pair)
  Q=rep(0, 2*npairs)
  
  LL=0
  for( t in 1:ntrials ) {
    if(ACC[t]<0) # catch missing values
      next
    
    # a is selected action in trial t 
    # correct action is always first of pair (A,C,E)
    A=2*pair[t]-1
    B=A+1
    a=A+(ACC[t]<1)
    #cat(sprintf("%i, %i, %i\n", A, B, a))
    Q[a] = Q[a] + alpha*(reward[t]-Q[a])
    #Pt=exp(Q[a]/beta)/( exp(Q[A]/beta)+exp(Q[B]/beta) )
    logPt=Q[a]/beta-log.sum.exp( c(Q[A], Q[B])/beta )
    LL=LL+logPt#log(Pt)
  }
  LL
}

## return probability correct under the model for each trial
paqlearning <- function( alpha, beta, ACC, pair, reward){
  npairs=max(pair)
  ntrials=length(pair)
  Q=rep(0, 2*npairs)
  pa=vector(length=ntrials)
  for( t in 1:ntrials ) {
    if(ACC[t]<0){ # catch missing values
      pa[t]=NA
      next
    }
    
    # a is selected action in trial t 
    # correct action is always first of pair (A,C,E)
    A=2*pair[t]-1
    B=A+1
    a=A+(ACC[t]<1)
    #cat(sprintf("%i, %i, %i\n", A, B, a))
    Q[a] = Q[a] + alpha*(reward[t]-Q[a])
    #Pt=exp(Q[a]/beta)/( exp(Q[A]/beta)+exp(Q[B]/beta) )
    logPt=Q[A]/beta-log.sum.exp( c(Q[A], Q[B])/beta )
    pa[t]=exp(logPt)
  }
  pa
}
