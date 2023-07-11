#include <RcppDist.h>
#include <RcppArmadillo.h>
#include <string>
// [[Rcpp::depends(RcppArmadillo, RcppDist, Rcpp)]]


// [[Rcpp::export]]

arma::mat xdist(arma::mat X,
                arma::vec mean0,
                arma::mat GPsigma,
                arma::vec addmean,
                float sdobs,
                std::string maxmin = "both",
                int K = 1000){
  
  int M = mean0.n_elem;
  int m = X.n_cols;
  
  arma::vec newmean(M);
  arma::mat ypred(K,M);
  
  arma::mat L = chol(GPsigma, "lower");
  
  arma::mat norm(M, 1);
  for(int i = 0; i < K; i++){
    newmean = mean0 + addmean * (arma::randn(1)*sdobs);
    norm = arma::randn(M);
    ypred.row(i) = trans(newmean + L*norm);
  }
  
  
  if(maxmin == "both"){
    m = 2*m;
  }
  
  arma::mat counts(K, m);
  
  if(maxmin == "max"){
    
    
    arma::ucolvec maxind = index_max(ypred, 1);
    
    for(int i = 0; i < K; i++){
      counts.row(i) = X.row(maxind(i));
    }
  } else if(maxmin == "min"){
    
    
    arma::ucolvec minind = index_min(ypred, 1);
    
    for(int i = 0; i < K; i++){
      counts.row(i) = X.row(minind(i));
    }
  } else if(maxmin == "both"){
    
    arma::ucolvec maxind = index_max(ypred, 1);
    arma::ucolvec minind = index_min(ypred, 1);
    
    for(int i = 0; i < K; i++){
      counts.row(i) = join_rows(X.row(maxind(i)), X.row(minind(i)));
    }
  }
  
  
  return counts;
  
}


