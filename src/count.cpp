#include <Rcpp.h>
using namespace Rcpp;


NumericMatrix checkRow( NumericMatrix data, NumericMatrix out,
                        int i, int numObs, int numVars,
                        int n_categories, bool freq = true ) {  
                          
  int r = 0;
  int num_true = 0;
  double value;
  
  if (freq) {
    value = 1.0;
  } else {
    value = 1.0 / numObs;
  }
  
  while (num_true <= (numVars-1)) {
    num_true = 0;
    for (int j = 0; j < numVars; j++) {
      if (data(i, j) == out(r, j)) num_true++;
    }
    r++;
    
    if (r >= n_categories) break;
  }
  
  if (num_true == numVars)
    out(r-1, numVars) = out(r-1, numVars) + value;
  
  if (i+1 <= numObs) {
    return checkRow( data, out, i+1, numObs, numVars, n_categories, freq );
  } else {
    return out;
  }
}


std::map<int, NumericVector> list_unique( NumericMatrix data ) {
  
  int numVars = data.ncol();
  int numObs = data.nrow();
  std::map<int, NumericVector> out;

  for (int j = 0; j < numVars; j++) {
    for (int i = 0; i < numObs; i++) {
        
      bool not_in = true;
      
      for (int u = 0; u < out[j].length(); u++) {
        if ( out[j][u] == data(i, j) )
          not_in = false;
      }
      
      if ( not_in )
        out[j].push_back( data(i, j) );
    }
    std::sort(out[j].begin(), out[j].end());
  }
  
  return out;
}


std::map<int, NumericVector> values_range( NumericMatrix data ) {

  int numVars = data.ncol();
  int numObs = data.nrow();
  int max_val;
  int min_val;
  std::map<int, NumericVector> out;
  
  for (int j = 0; j < numVars; j++) {
    for (int i = 0; i < numObs; i++) {
      if (fmod(data(i, j), 1) != 0)
        Rcpp::stop("not an integer");
      
      if (i == 0 || data(i, j) > max_val)
        max_val = data(i, j);
        
      if (i == 0 || data(i, j) < min_val)
        min_val = data(i, j);
    }
        
    for (int i = min_val; i < max_val+1; i++) {
      out[j].push_back(i);
    }  
  }
  
  return out;
}


// [[Rcpp::export]]


NumericMatrix count_cpp(NumericMatrix data, bool unique = true, bool freq = true) {
  
  int numVars = data.ncol();
  int numObs = data.nrow();
  NumericVector max_val(numVars);
  NumericVector cum_max(numVars);
  std::map<int, NumericVector> levels;
  
  if (unique) {
    levels = list_unique(data);
  } else {
    levels = values_range(data);
  }
  
  for (int j = 0; j < numVars; j++) {
    max_val[j] = levels[j].length();
  }
                          
  // create matrix of possible values
  
  cum_max[0] = max_val[0];
  
  for (int j = 1; j < numVars; j++) {
    cum_max[j] = max_val[j] * cum_max[j-1];
  }
    
  NumericMatrix out(cum_max[numVars-1], numVars+1);
  
  for (int r = 0; r < cum_max[numVars-1]; r++) {
    int val = fmod(r, max_val[0]); 
    out(r, 0) = levels[0][val]; 
    for (int j = 1; j < numVars; j++) {
      int val = fmod(ceil((r+1) / cum_max[j-1])-1, max_val[j]);
      out(r, j) = levels[j][val];  
    }
  }
    
  // search for values in the dataset
  
  out = checkRow( data, out, 0, numObs, numVars, cum_max[numVars-1], freq );
  return out;
}

