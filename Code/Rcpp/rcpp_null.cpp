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
	NumericVector obsage = env["obs.age"];
	NumericVector obsage2 = env["obs.age2"];
	
	IntegerVector idx1 = IntegerVector::create(12,13,14,15,16,17,18,19,20);
	IntegerVector idx2 = IntegerVector::create(21,22,23,24,25,26,27,28,29);
	IntegerVector idx3 = IntegerVector::create(30,31,32,33,34,35,36,37,38);
	IntegerVector idx4 = IntegerVector::create(43,44,45,46,47,48,49,50,51);
	
	double lambda1 = sumC(Bsplines.row(ind) * exp(Theta[idx1]));
	double Lambda1 = sumC(Isplines.row(ind) * exp(Theta[idx1]));

	double lambda2 = sumC(Bsplines.row(ind) * exp(Theta[idx2]));
	double Lambda2 = sumC(Isplines.row(ind) * exp(Theta[idx2]));

	double lambda3 = sumC(Bsplines.row(ind) * exp(Theta[idx3]));
	double Lambda3 = sumC(Isplines.row(ind) * exp(Theta[idx3]));

	double H = sumC(Isplines.row(ind) * exp(Theta[idx4]));
	
	double cause1 = pow(frailty*lambda1*exp(obsbmi2[ind]*Theta[0] + obsbmi3[ind]*Theta[1] + obsage[ind]*Theta[2] + obsage2[ind]*Theta[3]), obsDelta[ind] == 1); 
  	double cause2 = pow(frailty*lambda2*exp(obsbmi2[ind]*Theta[4] + obsbmi3[ind]*Theta[5] + obsage[ind]*Theta[6] + obsage2[ind]*Theta[7]), obsDelta[ind] == 2); 
  	double cause3 = pow(frailty*lambda3*exp(obsbmi2[ind]*Theta[8] + obsbmi3[ind]*Theta[9] + obsage[ind]*Theta[10] + obsage2[ind]*Theta[11]), obsDelta[ind] == 3); 


  	double cumcause1 =  frailty*Lambda1*exp(obsbmi2[ind]*Theta[0] + obsbmi3[ind]*Theta[1] + obsage[ind]*Theta[2] + obsage2[ind]*Theta[3]); 
  	double cumcause2 =  frailty*Lambda2*exp(obsbmi2[ind]*Theta[4] + obsbmi3[ind]*Theta[5] + obsage[ind]*Theta[6] + obsage2[ind]*Theta[7]); 
  	double cumcause3 =  frailty*Lambda3*exp(obsbmi2[ind]*Theta[8] + obsbmi3[ind]*Theta[9] + obsage[ind]*Theta[10] + obsage2[ind]*Theta[11]); 

  	double dtimel = cause1*cause2*cause3*exp( -(cumcause1 + cumcause2 + cumcause3));

  	double with = pow( 1 - exp( - frailty*H*exp(obsbmi2[ind]*Theta[39] + obsbmi3[ind]*Theta[40] + obsage[ind]*Theta[41] + obsage2[ind]*Theta[42])), obsU[ind] == 1);
  	double without = pow( exp(- frailty*H*exp(obsbmi2[ind]*Theta[39] + obsbmi3[ind]*Theta[40] + obsage[ind]*Theta[41] + obsage2[ind]*Theta[42])), obsU[ind] == 0); 
       
  	// current-status data model
  	double morbl = with*without;
   
  	// individual contribution to the likelihood.
  	double marginal = exp(Theta[52])*pow(frailty, exp(Theta[52])-1)*exp(-pow(frailty,exp(Theta[52])))*exp(frailty);
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
		indlike[i] = std::max(indlike[i], 1 / pow(10, 10));
	}
    
   	NumericVector loglike = log(indlike);
   	double negloglike = - sumC(loglike);

   	return(negloglike);
}

