// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <vector>
#include <string>
#include <algorithm>
#include <unordered_set>
#include <limits>
#include <cmath>
#include <sstream>

using namespace Rcpp;

struct CandSet {
  std::vector<int> coverageSet; 
  double sumRes;
  int coverage;
  double metric;     // coverage*(meanRes^2)
  std::string ruleText;
  std::vector<int> featUsed; 
};

struct BestRule {
  bool found;
  double metric;
  double meanResidual;
  int coverage;
  std::string ruleText;
  std::vector<int> inside;
  BestRule() : found(false), metric(-1e15), meanResidual(0.0), coverage(0) {}
};

inline std::vector<int> setIntersect(const std::vector<int> &A, const std::vector<int> &B) {
  std::vector<int> out; 
  out.reserve(std::min(A.size(), B.size()));
  size_t i=0, j=0;
  while(i<A.size() && j<B.size()){
    if(A[i]<B[j]) { i++; }
    else if(B[j]<A[i]){ j++; }
    else {
      out.push_back(A[i]);
      i++; j++;
    }
  }
  return out;
}

inline std::vector<int> setUnionV(const std::vector<int> &A, const std::vector<int> &B) {
  std::vector<int> out; 
  out.reserve(A.size()+B.size());
  size_t i=0, j=0;
  while(i<A.size() && j<B.size()){
    if(A[i]<B[j]) { out.push_back(A[i]); i++; }
    else if(B[j]<A[i]){ out.push_back(B[j]); j++; }
    else {
      out.push_back(A[i]);
      i++; j++;
    }
  }
  while(i<A.size()) { out.push_back(A[i]); i++; }
  while(j<B.size()) { out.push_back(B[j]); j++; }
  return out;
}

// sum of residuals for coverageSet
inline double sumResiduals(const std::vector<int> &cover, const NumericVector &resid) {
  double s=0;
  for(int idx : cover){ s += resid[idx]; }
  return s;
}

//-----------------------------------------------------
// buildSingleLiterals: same as before
//-----------------------------------------------------
std::vector<CandSet> buildSingleLiterals(const NumericMatrix &X,
                                         const NumericVector &resid) {
  int n = X.nrow();
  int p = X.ncol();
  std::vector<CandSet> out(2*p);
  
#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int j=0; j<p; j++){
    std::vector<int> zeroIdx; zeroIdx.reserve(n);
    std::vector<int> oneIdx;  oneIdx.reserve(n);
    for(int i=0; i<n; i++){
      int xij = (X(i,j)>0.5)?1:0;
      if(xij==0) zeroIdx.push_back(i);
      else        oneIdx.push_back(i);
    }
    double sumZ=0, sumO=0;
    for(auto &r: zeroIdx) sumZ+=resid[r];
    for(auto &r: oneIdx)  sumO+=resid[r];
    
    // (X_j==0)
    {
      CandSet c;
      c.coverageSet = std::move(zeroIdx);
      c.coverage    = c.coverageSet.size();
      double meanr  = (c.coverage>0)?(sumZ/double(c.coverage)):0.0;
      c.sumRes      = sumZ;
      c.metric      = double(c.coverage)*(meanr*meanr);
      std::ostringstream oss; 
      oss << "(X" << (j+1) << "==0)";
      c.ruleText    = oss.str();
      c.featUsed    = {j};
      out[2*j]      = std::move(c);
    }
    // (X_j==1)
    {
      CandSet c;
      c.coverageSet = std::move(oneIdx);
      c.coverage    = c.coverageSet.size();
      double meanr  = (c.coverage>0)?(sumO/double(c.coverage)):0.0;
      c.sumRes      = sumO;
      c.metric      = double(c.coverage)*(meanr*meanr);
      std::ostringstream oss; 
      oss << "(X" << (j+1) << "==1)";
      c.ruleText    = oss.str();
      c.featUsed    = {j};
      out[2*j+1]    = std::move(c);
    }
  }
  
  return out;
}

//-----------------------------------------------------
// keepTopK: Given a vector of CandSet, keep only top K by c.metric
//-----------------------------------------------------
std::vector<CandSet> keepTopK(const std::vector<CandSet> &all, int K) {
  if((int)all.size() <= K) {
    return all;  // already small
  }
  // partial sort
  std::vector<CandSet> out = all;  // copy
  std::nth_element(out.begin(), out.begin()+K, out.end(),
                   [](const CandSet &a, const CandSet &b){
                     return a.metric > b.metric;
                   });
  out.resize(K);
  // we might want to fully sort the top K
  std::sort(out.begin(), out.end(), [](auto &a, auto &b){
    return a.metric > b.metric;
  });
  return out;
}

//-----------------------------------------------------
// buildTwoLiteralSets_topK: build 2-literal from singles but only keep topK
//-----------------------------------------------------
std::vector<CandSet> buildTwoLiteralSets_topK(
    const std::vector<CandSet> &singles,
    const NumericVector &resid,
    int n,
    int K
) {
  // we'll collect *all* 2-literal, then keep topK. Or we can do chunk-based partial merges.
  // For demonstration, let's do full, then keep topK.
  
  std::vector<CandSet> results;
  results.reserve(singles.size()*singles.size()/4); // guess
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<CandSet> localVec;
  localVec.reserve(1000);
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(int i=0; i<(int)singles.size()-1; i++){
    const CandSet &A = singles[i];
    // We skip expansions if coverage is 0
    if(A.coverage==0) continue;
    for(int j=i+1; j<(int)singles.size(); j++){
      const CandSet &B = singles[j];
      if(B.coverage==0) continue;
      // distinct features?
      if(A.featUsed[0]==B.featUsed[0]) continue;
      
      // AND
      {
        std::vector<int> cov = setIntersect(A.coverageSet,B.coverageSet);
        int c = cov.size();
        if(c>0){
          double sumr=0.0;
          for(int idx: cov) sumr += resid[idx];
          double meanr = sumr/double(c);
          CandSet newC;
          newC.coverageSet = std::move(cov);
          newC.coverage    = c;
          newC.sumRes      = sumr;
          newC.metric      = double(c)*(meanr*meanr);
          std::ostringstream oss;
          oss<<"("<<A.ruleText<<" & "<<B.ruleText<<")";
          newC.ruleText = oss.str();
          newC.featUsed = { A.featUsed[0], B.featUsed[0] };
          localVec.push_back(std::move(newC));
        }
      }
      // OR
      {
        std::vector<int> cov = setUnionV(A.coverageSet,B.coverageSet);
        int c = cov.size();
        if(c>0){
          double sumr=0.0;
          for(int idx: cov) sumr += resid[idx];
          double meanr = sumr/double(c);
          CandSet newC;
          newC.coverageSet = std::move(cov);
          newC.coverage    = c;
          newC.sumRes      = sumr;
          newC.metric      = double(c)*(meanr*meanr);
          std::ostringstream oss;
          oss<<"("<<A.ruleText<<" | "<<B.ruleText<<")";
          newC.ruleText = oss.str();
          newC.featUsed = { A.featUsed[0], B.featUsed[0] };
          localVec.push_back(std::move(newC));
        }
      }
      
    } // j
  }   // i
  
#ifdef _OPENMP
#pragma omp critical
#endif
{
  results.insert(results.end(), localVec.begin(), localVec.end());
}
} // parallel

// keep topK
return keepTopK(results, K);
}

//-----------------------------------------------------
// buildThreeLiteralSets_topK: from topK 2-literal, expand with single-literal
//-----------------------------------------------------
std::vector<CandSet> buildThreeLiteralSets_topK(
    const std::vector<CandSet> &twos, 
    const std::vector<CandSet> &singles,
    const NumericVector &resid,
    int n,
    int K
){
  std::vector<CandSet> results;
  results.reserve(twos.size()*singles.size()/4);
  
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<CandSet> localVec;
  localVec.reserve(1000);
  
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(int i=0; i<(int)twos.size(); i++){
    const CandSet &A = twos[i];
    if(A.coverage==0) continue;
    int f1 = A.featUsed[0];
    int f2 = A.featUsed[1];
    for(int j=0; j<(int)singles.size(); j++){
      const CandSet &B = singles[j];
      if(B.coverage==0) continue;
      int fB = B.featUsed[0];
      // distinct features => skip if overlap
      if(fB==f1 || fB==f2) continue;
      
      // AND
      {
        std::vector<int> cov = setIntersect(A.coverageSet,B.coverageSet);
        int c = cov.size();
        if(c>0){
          double sumr=0;
          for(int idx: cov){ sumr+=resid[idx]; }
          double meanr = sumr/double(c);
          CandSet newC;
          newC.coverageSet = std::move(cov);
          newC.coverage    = c;
          newC.sumRes      = sumr;
          newC.metric      = double(c)*(meanr*meanr);
          std::ostringstream oss;
          oss<<"("<<A.ruleText<<" & "<<B.ruleText<<")";
          newC.ruleText = oss.str();
          newC.featUsed = {f1,f2,fB};
          localVec.push_back(std::move(newC));
        }
      }
      // OR
      {
        std::vector<int> cov = setUnionV(A.coverageSet,B.coverageSet);
        int c = cov.size();
        if(c>0){
          double sumr=0;
          for(int idx: cov){ sumr+=resid[idx]; }
          double meanr = sumr/double(c);
          CandSet newC;
          newC.coverageSet = std::move(cov);
          newC.coverage    = c;
          newC.sumRes      = sumr;
          newC.metric      = double(c)*(meanr*meanr);
          std::ostringstream oss;
          oss<<"("<<A.ruleText<<" | "<<B.ruleText<<")";
          newC.ruleText = oss.str();
          newC.featUsed = {f1,f2,fB};
          localVec.push_back(std::move(newC));
        }
      }
    }
  }
  
#ifdef _OPENMP
#pragma omp critical
#endif
{
  results.insert(results.end(), localVec.begin(), localVec.end());
}
} // parallel

return keepTopK(results, K);
}

//-----------------------------------------------------
// pickBest: coverage among unassigned rows
//-----------------------------------------------------
BestRule pickBest(const std::vector<CandSet> &cands,
                  const NumericVector &resid,
                  const IntegerVector &unassigned) {
  // build hash set of unassigned rows
  std::unordered_set<int> unassignedSet;
  unassignedSet.reserve(unassigned.size());
  for(int i=0; i<unassigned.size(); i++){
    unassignedSet.insert(unassigned[i]);
  }
  BestRule best;
  for(const auto &cs: cands){
    double sumr=0;
    std::vector<int> inside;
    inside.reserve(cs.coverageSet.size());
    for(int idx: cs.coverageSet){
      if(unassignedSet.find(idx)!=unassignedSet.end()){
        inside.push_back(idx);
        sumr+=resid[idx];
      }
    }
    int c = inside.size();
    if(c>0){
      double meanr = sumr/double(c);
      double metric = double(c)*(meanr*meanr);
      if(metric>best.metric){
        best.found=true;
        best.metric=metric;
        best.meanResidual=meanr;
        best.coverage=c;
        best.ruleText=cs.ruleText;
        best.inside=std::move(inside);
      }
    }
  }
  return best;
}

//-----------------------------------------------------
// main function
// [[Rcpp::export]]
List find_best_binary_rule_3way_topK(
    const NumericMatrix &X,
    const NumericVector &residuals,
    const IntegerVector &unassigned,
    int K_twoWay = 100,
    int K_threeWay = 100
) {
  int n = X.nrow();
  int p = X.ncol();
  
  // 1) single-literal
  std::vector<CandSet> singles = buildSingleLiterals(X, residuals);
  
  // 2) two-literal => keep top K
  std::vector<CandSet> twos = buildTwoLiteralSets_topK(singles, residuals, n, K_twoWay);
  
  // 3) three-literal => keep top K
  std::vector<CandSet> threes = buildThreeLiteralSets_topK(twos, singles, residuals, n, K_threeWay);
  
  // combine singles, twos, threes
  // (optionally keep top K from this union, but let's just pick the best among them)
  std::vector<CandSet> all;
  all.reserve(singles.size()+twos.size()+threes.size());
  all.insert(all.end(), singles.begin(), singles.end());
  all.insert(all.end(), twos.begin(),   twos.end());
  all.insert(all.end(), threes.begin(), threes.end());
  
  // pick best w.r.t. coverage among unassigned
  BestRule best = pickBest(all, residuals, unassigned);
  if(!best.found){
    return List::create(_["found_rule"]=false);
  }
  IntegerVector insideR(best.inside.size());
  for(int i=0; i<(int)best.inside.size(); i++){
    insideR[i] = best.inside[i];
  }
  return List::create(
    _["found_rule"]=true,
    _["rule_string"]=best.ruleText,
    _["coverage"]=best.coverage,
    _["mean_residual"]=best.meanResidual,
    _["metric"]=best.metric,
    _["points_inside"]=insideR
  );
}
