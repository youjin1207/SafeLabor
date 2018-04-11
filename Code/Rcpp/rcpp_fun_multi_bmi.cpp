#include <Rcpp.h>
#include <cmath>
//[[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;

//[[Rcpp::export]]
double sumC(NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}


//[[Rcpp::export]]
double splines(NumericVector Theta){
	// 'ind' starts from 0(zero) 
	// global variables in R
	Environment env = Environment::global_env();
	// NumericMatrix Bsplines = env["Bsplines"]; 
	// NumericVector obsX = env["obs.X"];

	return(1.23);
}


//[[Rcpp::export]]
double joint(double frailty, NumericVector Theta, int ind){
	// 'ind' starts from 0(zero) 
	// global variables in R
	Environment env = Environment::global_env();
	NumericMatrix Bsplines = env["Bsplines"]; 
	NumericMatrix Isplines = env["Isplines"]; 
	IntegerVector obsDelta = env["obs.Delta"];
	IntegerVector obsU = env["obs.U"];
	NumericVector obsbmi2 = env["obs.bmi2"];
	NumericVector obsbmi3 = env["obs.bmi3"];
	
	IntegerVector idx1 = IntegerVector::create(6,7,8,9,10,11,12,13,14);
	IntegerVector idx2 = IntegerVector::create(15,16,17,18,19,20,21,22,23);
	IntegerVector idx3 = IntegerVector::create(24,25,26,27,28,29,30,31,32);
	
	double lambda1 = sumC(Bsplines.row(ind) * exp(Theta[idx1]));
	double Lambda1 = sumC(Isplines.row(ind) * exp(Theta[idx1]));

	double lambda2 = sumC(Bsplines.row(ind) * exp(Theta[idx2]));
	double Lambda2 = sumC(Isplines.row(ind) * exp(Theta[idx2]));

	double H = sumC(Isplines.row(ind) * exp(Theta[idx3]));
	
	double cause1 = pow(exp(-frailty)*lambda1*exp(obsbmi2[ind]*Theta[0] + obsbmi3[ind]*Theta[1]), obsDelta[ind] == 1); 
  	double cause2 = pow(exp(-frailty)*lambda2*exp(obsbmi2[ind]*Theta[2] + obsbmi3[ind]*Theta[3]), obsDelta[ind] == 2); 
  	
  	double cumcause1 =  exp(-frailty)*Lambda1*exp(obsbmi2[ind]*Theta[0] + obsbmi3[ind]*Theta[1]); 
  	double cumcause2 =  exp(-frailty)*Lambda2*exp(obsbmi2[ind]*Theta[2] + obsbmi3[ind]*Theta[3]); 
  	
  	double dtimel = cause1*cause2*exp( -(cumcause1 + cumcause2));

  	double with = pow( 1 - exp( - exp(-frailty)*H*exp(obsbmi2[ind]*Theta[4] + obsbmi3[ind]*Theta[5])), obsU[ind] == 1);
  	double without = pow( exp(- exp(-frailty)*H*exp(obsbmi2[ind]*Theta[4] + obsbmi3[ind]*Theta[5])), obsU[ind] == 0); 
       
  	// current-status data model
  	double morbl = with*without;
   
  	// individual contribution to the likelihood.
  	double marginal = exp(-pow(frailty,2) /(2*pow(exp(Theta[33]),2)) + pow(frailty,2))/sqrt(2*M_PI*pow(exp(Theta[33]),2));
  	double indcontri =  dtimel*morbl*marginal;

	return(indcontri);
}


// [[Rcpp::export]]
double loglikeli(NumericVector parameter){

	// global variables in R
	Environment env = Environment::global_env();
	NumericVector nodes = env["nodes"]; 
	NumericVector weights = env["weights"]; 
	const int n = env["n"];

	IntegerVector id = seq_len(n)-1;
	NumericVector indlike(n);  

	for(int i = 0; i < n; ++i){
		for(int w = 0; w < weights.size(); ++w) {
			indlike[i] += joint(nodes[w], parameter, id[i])*weights[w];
	  }
		indlike[i] = std::max(indlike[i], 1 / pow(10, 20));
	}
    
   	NumericVector loglike = log(indlike);
   	double negloglike = - sumC(loglike);

   	return(negloglike);
}