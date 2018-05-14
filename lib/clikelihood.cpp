#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double csimandeval(double alpha, double beta, int ntrials){
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

// [[Rcpp::export]]
double csimandeval_posneg(double alphaG, double alphaL, double beta, int ntrials){
  int npairs=3;
  double reward_prob[3]={.8, .7, .6};
  double ACC=0.0;
  int pair, t;
  double Q[6]={0.0};
  int A,B,a,reward;
  double pA, offset, pe;
  
  for(pair=0; pair<npairs; pair++){
    for(t=0; t<ntrials; t++){
      A=2*pair;
      B=A+1;
      offset=((Q[A]/beta)>(Q[B]/beta) ? (Q[A]/beta) : (Q[B]/beta));
      pA=exp(Q[A]/beta - (log( exp( Q[A]/beta - offset)+exp( Q[B]/beta - offset)) + offset));
      a=(((double)runif(1)[0])<pA) ? A : B;
      reward=(((double)runif(1)[0])<reward_prob[pair]) ? (a==A) : (a==B);
      //Rprintf("pA=%f, a=%i, reward=%i, ACC=%i\\n", pA, a, reward, a==A);

      pe=(reward-Q[a]);
      Q[a] = Q[a] + ((pe>0)?(alphaG):(alphaL))*pe;      
      ACC = ACC+(a==A);
    }
  }
  return ACC/(ntrials*npairs);
}

/** model w learning rate for positive and negative feedback
*/ 
// [[Rcpp::export]]
double cllqlearning_posneg(double alphaG, double alphaL, double beta, NumericVector ACC, NumericVector pair, NumericVector reward){
  int npairs=max(pair);
  int ntrials=pair.size();
  
  NumericVector Q=NumericVector(2*npairs, 0.0);
  double LL=0.0;
  double logPt;
  
  int A,B,a;
  double pe; // prediction error
  double offset; // for logsumexp trick
  
  for( int t=0; t<ntrials; t++ ){
    if(ACC[t]<0) // missing responses
      continue;
    
    // a is selected action in trial t
    // correct action is always first of pair (A,C,E)
    A=2*(((int)pair[t])-1);
    B=A+1;
    a=A+(ACC[t]<1);
    
    pe=(reward[t]-Q[a]);
    Q[a] = Q[a] + ((pe>0)?(alphaG):(alphaL))*pe;
    // printf("t=%i, Q[%i]=%f\n",t, a, Q[a]);
    
    // this is what is computed but using the log.sum.exp trick
    //Pt=exp(Q[a]/beta)/( exp(Q[A]/beta)+exp(Q[B]/beta) );
    //LL+=log(Pt);
    offset=(((Q[A]/beta)>(Q[B]/beta))?(Q[A]/beta):(Q[B]/beta));
    logPt=Q[a]/beta - (log( exp( Q[A]/beta - offset) + exp( Q[B]/beta - offset)) + offset);
    LL+=logPt;
  }
  return(LL);
}
