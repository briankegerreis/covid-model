#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List blockListCpp(NumericMatrix B,NumericVector lam,IntegerVector typ) {
  int N = lam.size();
  std::vector<int> to;
  std::vector<int> from;
  for (int i = 0; i < N; i++){
    for (int j = i+1; j < N; j++){
      float p = pow(B(typ[i],typ[j]),1+lam[i]+lam[j]);
      if(R::rbinom(1,p) == 1){
        to.push_back(j+1);
        from.push_back(i+1);
      }
    }
  }
  return List::create(_["from"] = from,
                      _["to"] = to);
}
