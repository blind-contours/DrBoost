// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>

using namespace Rcpp;

// A bounding box in p-dimensional space
struct RegionBox {
  std::vector<double> lower, upper;
  RegionBox() {}
  RegionBox(int p){
    lower.resize(p, -std::numeric_limits<double>::infinity());
    upper.resize(p,  std::numeric_limits<double>::infinity());
  }
};

// A candidate "beam" item
struct CandBox {
  RegionBox box;
  double sumRes   = 0.0;
  int coverage    = 0;
  double metric   = -1e15; // coverage * (meanRes^2)
  std::vector<int> inside; // row indices (0-based)
};

static void evaluateCandidate(
    const RegionBox &cand,
    const NumericMatrix &X,
    const NumericVector &residuals,
    const IntegerVector &coverageCount,
    int max_overlap,
    int minCount,
    int maxCount,
    std::vector<CandBox> &collector
) {
  int n = X.nrow();
  int p = X.ncol();
  std::vector<int> inside;
  inside.reserve(n);
  
  double sumr = 0.0;
  for(int i=0; i<n; i++){
    if(coverageCount[i] >= max_overlap) {
      continue;  // skip if we exceed overlap
    }
    bool inBox = true;
    for(int d=0; d<p; d++){
      double xv = X(i,d);
      if(xv < cand.lower[d] || xv > cand.upper[d]) {
        inBox = false;
        break;
      }
    }
    if(inBox){
      inside.push_back(i);
      sumr += residuals[i];
    }
  }
  
  int covSize = (int)inside.size();
  if(covSize < minCount) return;
  if(covSize > maxCount) return;
  
  // coverage * (mean residual)^2
  double meanr  = sumr / (double)covSize;
  double metric = (double)covSize * (meanr * meanr);
  if(metric <= 0.0) return;
  
  // Build the candidate
  CandBox cb;
  cb.box      = cand;
  cb.sumRes   = sumr;
  cb.coverage = covSize;
  cb.metric   = metric;
  cb.inside   = std::move(inside);
  collector.push_back(std::move(cb));
}

// Simple helper to keep top K by 'metric'
static void keepTopK(std::vector<CandBox> &vec, int K) {
  if((int)vec.size() <= K) return;
  std::nth_element(
    vec.begin(),
    vec.begin() + K,
    vec.end(),
    [](const CandBox &A, const CandBox &B){
      return A.metric > B.metric; // descending
    }
  );
  vec.resize(K);
  std::sort(
    vec.begin(),
    vec.end(),
    [](const CandBox &A, const CandBox &B){
      return A.metric > B.metric;
    }
  );
}

static std::vector<double> collectThresholds(
    const NumericMatrix &thresh,
    int jcol
) {
  std::vector<double> vals;
  vals.reserve(thresh.nrow());
  for(int r=0; r<thresh.nrow(); r++){
    double tv = thresh(r,jcol);
    if(!NumericVector::is_na(tv)){
      vals.push_back(tv);
    }
  }
  if(!vals.empty()) {
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
  }
  return vals;
}

// [[Rcpp::export]]
List find_topk_candidates_beam_cpp(
    const NumericMatrix &X,
    const NumericVector &residuals,
    const IntegerVector &coverageCount,
    const int max_overlap,
    const NumericMatrix &thresholds_list,
    double min_obs_pct=0.05,
    double max_obs_frac=1.0,
    int K1=200, int K2=200, int K3=200, int K4=200,
    int topK=50
){
  // 1) basic checks
  int n = X.nrow();
  int p = X.ncol();
  if(n < 1 || p < 1){
    return List::create(_["found_any"] = false);
  }
  
  int usableCount = 0;
  for(int i=0; i<n; i++){
    if(coverageCount[i] < max_overlap){
      usableCount++;
    }
  }
  if(usableCount == 0){
    return List::create(_["found_any"] = false);
  }
  int minCount = std::max(1, (int)std::floor(usableCount * min_obs_pct));
  int maxCount = (int)std::floor(usableCount * max_obs_frac);
  if(maxCount < minCount){
    return List::create(_["found_any"] = false);
  }
  
  // 2) gather thresholds
  std::vector< std::vector<double> > allThresh(p);
  for(int j=0; j<p; j++){
    allThresh[j] = collectThresholds(thresholds_list, j);
  }
  
  // We'll do 1D,2D,3D,4D expansions using beams
  std::vector<CandBox> beam1, beam2, beam3, beam4;
  beam1.reserve(10000);
  beam2.reserve(10000);
  beam3.reserve(10000);
  beam4.reserve(10000);
  
  // A) 1D expansions
  for(int j=0; j<p; j++){
    const auto &tj = allThresh[j];
    if(tj.empty()) continue;
    int sz = (int)tj.size();
    for(int a=0; a<sz; a++){
      for(int b=a; b<sz; b++){
        RegionBox cand(p);
        for(int d=0; d<p; d++){
          cand.lower[d] = -std::numeric_limits<double>::infinity();
          cand.upper[d] =  std::numeric_limits<double>::infinity();
        }
        cand.lower[j] = tj[a];
        cand.upper[j] = tj[b];
        evaluateCandidate(cand, X, residuals, coverageCount,
                          max_overlap, minCount, maxCount, beam1);
      }
    }
  }
  keepTopK(beam1, K1);
  
  // B) 2D expansions
  for(size_t i=0; i<beam1.size(); i++){
    const RegionBox &b1 = beam1[i].box;
    for(int j2=0; j2<p; j2++){
      const auto &t2 = allThresh[j2];
      if(t2.empty()) continue;
      int sz2 = (int)t2.size();
      for(int a=0; a<sz2; a++){
        for(int b=a; b<sz2; b++){
          RegionBox cand = b1;
          cand.lower[j2] = t2[a];
          cand.upper[j2] = t2[b];
          evaluateCandidate(cand, X, residuals, coverageCount,
                            max_overlap, minCount, maxCount, beam2);
        }
      }
    }
  }
  keepTopK(beam2, K2);
  
  // C) 3D expansions
  if(K3>0){
    for(size_t i=0; i<beam2.size(); i++){
      const RegionBox &b2 = beam2[i].box;
      for(int j3=0; j3<p; j3++){
        const auto &t3 = allThresh[j3];
        if(t3.empty()) continue;
        int sz3 = (int)t3.size();
        for(int a=0; a<sz3; a++){
          for(int b=a; b<sz3; b++){
            RegionBox cand = b2;
            cand.lower[j3] = t3[a];
            cand.upper[j3] = t3[b];
            evaluateCandidate(cand, X, residuals, coverageCount,
                              max_overlap, minCount, maxCount, beam3);
          }
        }
      }
    }
    keepTopK(beam3, K3);
  }
  
  // D) 4D expansions
  if(K4>0){
    for(size_t i=0; i<beam3.size(); i++){
      const RegionBox &b3 = beam3[i].box;
      for(int j4=0; j4<p; j4++){
        const auto &t4 = allThresh[j4];
        if(t4.empty()) continue;
        int sz4 = (int)t4.size();
        for(int a=0; a<sz4; a++){
          for(int b=a; b<sz4; b++){
            RegionBox cand = b3;
            cand.lower[j4] = t4[a];
            cand.upper[j4] = t4[b];
            evaluateCandidate(cand, X, residuals, coverageCount,
                              max_overlap, minCount, maxCount, beam4);
          }
        }
      }
    }
    keepTopK(beam4, K4);
  }
  
  // E) merge all
  std::vector<CandBox> allCands;
  allCands.reserve(beam1.size() + beam2.size()
                     + beam3.size() + beam4.size());
  allCands.insert(allCands.end(), beam1.begin(), beam1.end());
  allCands.insert(allCands.end(), beam2.begin(), beam2.end());
  if(K3>0) allCands.insert(allCands.end(), beam3.begin(), beam3.end());
  if(K4>0) allCands.insert(allCands.end(), beam4.begin(), beam4.end());
  
  // keep topK overall
  keepTopK(allCands, topK);
  if(allCands.empty()){
    return List::create(_["found_any"] = false);
  }
  
  // Build return
  // We'll return a list of boxes: each has "lb","ub","coverage","sumRes","insideIdx"
  // in ascending order of rank
  List out;
  out["found_any"] = true;
  
  // We'll store them in a structured manner
  int outSize = (int)allCands.size();
  NumericMatrix LBs(outSize, p), UBs(outSize, p);
  IntegerVector coverVec(outSize);
  NumericVector sumResVec(outSize), metricVec(outSize);
  List insideIdxList(outSize);
  
  for(int i=0; i<outSize; i++){
    const auto &cb = allCands[i];
    for(int d=0; d<p; d++){
      LBs(i,d) = cb.box.lower[d];
      UBs(i,d) = cb.box.upper[d];
    }
    coverVec[i]   = cb.coverage;
    sumResVec[i]  = cb.sumRes;
    metricVec[i]  = cb.metric;
    
    // convert inside to an integer vector
    IntegerVector indices(cb.inside.size());
    for(int k=0; k<(int)cb.inside.size(); k++){
      indices[k] = cb.inside[k];
    }
    insideIdxList[i] = indices;
  }
  
  out["lowerBounds"] = LBs;
  out["upperBounds"] = UBs;
  out["coverage"]    = coverVec;
  out["sumRes"]      = sumResVec;
  out["metric"]      = metricVec;
  out["inside"]      = insideIdxList;
  
  return out;
}
