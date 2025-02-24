#include <Rcpp.h>
#include <unordered_set>
#include <limits>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
List find_best_logical_region_cpp(const NumericMatrix& X, // binary in {0,1}
                                  const NumericVector& residuals,
                                  const IntegerVector& used_rows, 
                                  double min_obs_pct = 0.05) {
  
  int n = X.nrow();
  int p = X.ncol();
  int min_obs = std::max((int)std::floor(n * min_obs_pct), 5);
  
  // Convert used_rows to a set, meaning these rows are "already assigned".
  std::unordered_set<int> assigned;
  for(int i=0; i<used_rows.size(); i++){
    assigned.insert( used_rows[i] );
  }
  
  double best_metric = -std::numeric_limits<double>::infinity();
  std::string best_rule;
  double best_meanres = 0.0;
  std::vector<int> best_points;
  
  // Try all pairs of binary vars, each in {0 or 1} form, plus AND/OR
  // For simplicity: (var1Val == a) AND/OR (var2Val == b)
  for(int var1=0; var1<p-1; var1++){
    for(int var2=var1+1; var2<p; var2++){
      for(int val1=0; val1<=1; val1++){ // 0 or 1
        for(int val2=0; val2<=1; val2++){
          // AND
          {
            std::vector<int> in_idx;
            double sumr=0.0;
            for(int i=0; i<n; i++){
              if(assigned.find(i) != assigned.end()) continue; // skip already used
              int x1 = (X(i,var1)>0.5)?1:0;
              int x2 = (X(i,var2)>0.5)?1:0;
              if(x1==val1 && x2==val2){
                in_idx.push_back(i);
                sumr += residuals[i];
              }
            }
            if((int)in_idx.size()>=min_obs){
              double meanr = sumr/(double)in_idx.size();
              double m = std::fabs(meanr)*std::sqrt((double)in_idx.size());
              if(m>best_metric){
                best_metric=m;
                best_meanres=meanr;
                best_points=in_idx;
                best_rule = "X" + std::to_string(var1+1) + "==" + std::to_string(val1) +
                  " & " +
                  "X" + std::to_string(var2+1) + "==" + std::to_string(val2);
              }
            }
          }
          // OR
          {
            std::vector<int> in_idx;
            double sumr=0.0;
            for(int i=0; i<n; i++){
              if(assigned.find(i)!=assigned.end()) continue;
              int x1 = (X(i,var1)>0.5)?1:0;
              int x2 = (X(i,var2)>0.5)?1:0;
              // Condition: (x1==val1) OR (x2==val2)
              if(x1==val1 || x2==val2){
                in_idx.push_back(i);
                sumr += residuals[i];
              }
            }
            if((int)in_idx.size()>=min_obs){
              double meanr = sumr/(double)in_idx.size();
              double m = std::fabs(meanr)*std::sqrt((double)in_idx.size());
              if(m>best_metric){
                best_metric=m;
                best_meanres=meanr;
                best_points=in_idx;
                best_rule = "(X" + std::to_string(var1+1) + "==" + std::to_string(val1) +
                  ") | (X" + std::to_string(var2+1) + "==" + std::to_string(val2) +")";
              }
            }
          }
        }// val2
      }// val1
    }
  }
  
  if(best_metric<0){
    return List::create(_["found_region"]=false);
  }
  
  return List::create(
    _["found_region"] = true,
    _["rule"] = best_rule,
    _["mean_residual"] = best_meanres,
    _["metric"] = best_metric,
    _["points_inside"] = best_points,
    _["continuous_pair"] = false
  );
}
