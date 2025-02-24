// [[Rcpp::plugins(openmp)]]
#include <Rcpp.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <vector>
#include <string>
#include <sstream>
#include <limits>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

/* --------------------------------------------------------------------------
 1) Data Structures
 -------------------------------------------------------------------------- */

// A bounding box in p-dimensional space
// lower[d] = -∞ means "no lower bound" on dimension d.
// upper[d] = +∞ means "no upper bound" on dimension d.
struct RegionBox {
  std::vector<double> lower, upper;
  RegionBox() {}
  RegionBox(int p) {
    lower.resize(p, -std::numeric_limits<double>::infinity());
    upper.resize(p,  std::numeric_limits<double>::infinity());
  }
};

struct CandBox {
  RegionBox box;
  double sumRes;
  int coverage;
  double metric;  // coverage * (meanResid)^2
  std::vector<int> inside;
  
  CandBox() : sumRes(0.0), coverage(0), metric(-1e15) {}
};


/* --------------------------------------------------------------------------
 2) Evaluate coverage (only among unassigned points) and compute metric
 -------------------------------------------------------------------------- */
inline void evaluateCandidate(const RegionBox &candBox,
                              const NumericMatrix &X,
                              const NumericVector &residuals,
                              const LogicalVector &assignedMask,
                              int minCount,
                              int maxCount,  // coverage upper bound
                              CandBox &bestSoFar)
{
  int n = X.nrow();
  int p = X.ncol();
  
  std::vector<int> inside;
  inside.reserve(n);
  
  double sumr = 0.0;
  // For each row i, skip if assigned or if outside candBox.
  for(int i=0; i<n; i++){
    if(assignedMask[i]) {
      continue;  // already used in an earlier region
    }
    bool inBox = true;
    for(int d=0; d<p; d++){
      double xval = X(i,d);
      if(xval < candBox.lower[d] || xval > candBox.upper[d]){
        inBox = false;
        break;
      }
    }
    if(inBox){
      inside.push_back(i);
      sumr += residuals[i];
    }
  }
  
  // Check minCount
  int covSize = (int)inside.size();
  if(covSize < minCount) return;
  
  // Check maxCount
  if(covSize > maxCount) return;
  
  double coverageD = (double)covSize;
  double meanr     = sumr / coverageD;
  double metric    = coverageD * (meanr * meanr);
  
  // Update if better
  if(metric > bestSoFar.metric){
    bestSoFar.metric   = metric;
    bestSoFar.sumRes   = sumr;
    bestSoFar.coverage = covSize;
    bestSoFar.box      = candBox;
    bestSoFar.inside   = std::move(inside);
  }
}


/* --------------------------------------------------------------------------
 3) Single-/Two-/Three-/Four-D expansions + beam search
 -------------------------------------------------------------------------- */

// Sort & unique a single column of thresholds
static std::vector<double> collectThresholds(const NumericMatrix &thresh,
                                             int jcol)
{
  std::vector<double> vals;
  vals.reserve(thresh.nrow());
  for(int r=0; r<(int)thresh.nrow(); r++){
    double tv = thresh(r,jcol);
    if(!NumericVector::is_na(tv)){
      vals.push_back(tv);
    }
  }
  if(!vals.empty()){
    std::sort(vals.begin(), vals.end());
    vals.erase(std::unique(vals.begin(), vals.end()), vals.end());
  }
  return vals;
}

// Keep top-K by metric
std::vector<CandBox> keepTopK(const std::vector<CandBox> &all, int K)
{
  if((int)all.size() <= K) return all;
  
  std::vector<CandBox> out = all;
  // partial sort
  std::nth_element(out.begin(), out.begin()+K, out.end(),
                   [](const CandBox &a, const CandBox &b){
                     return a.metric > b.metric;
                   });
  out.resize(K);
  // fully sort top-K
  std::sort(out.begin(), out.end(),
            [](const CandBox &a, const CandBox &b){
              return a.metric > b.metric;
            });
  return out;
}


/* --------------------------------------------------------------------------
 4) Helpers for excluding boxes
 -------------------------------------------------------------------------- */

// Compare two RegionBoxes for equality, dimension by dimension
bool sameBox(const RegionBox &a, const RegionBox &b, double eps = 1e-14)
{
  int p = (int)a.lower.size();
  for(int d=0; d<p; d++){
    // We treat them as "same" if lower and upper match exactly
    // (or within a tiny numerical tolerance).
    double diff1 = std::fabs(a.lower[d] - b.lower[d]);
    double diff2 = std::fabs(a.upper[d] - b.upper[d]);
    if(diff1 > eps || diff2 > eps) {
      return false;
    }
  }
  return true;
}

// Check if a candidate box is in the excluded list
bool isExcluded(const RegionBox &cand, const std::vector<RegionBox> &excludedVec)
{
  for(const auto &badBox : excludedVec){
    if(sameBox(cand, badBox)) {
      return true;
    }
  }
  return false;
}


/* --------------------------------------------------------------------------
 5) Converting box -> textual rule
 -------------------------------------------------------------------------- */

// Convert double => string, with +/-Inf
static std::string doubleToStr(double val){
  if(std::isinf(val)){
    return (val > 0.0) ? "+Inf" : "-Inf";
  }
  std::ostringstream ss;
  ss << val;
  return ss.str();
}

// Make a dimension constraint
static std::string makeOneDimConstraint(const std::string &featName,
                                        double lowVal,
                                        double highVal)
{
  bool hasLower = (lowVal != -std::numeric_limits<double>::infinity());
  bool hasUpper = (highVal!=  std::numeric_limits<double>::infinity());
  
  if(hasLower && hasUpper){
    std::ostringstream ss;
    ss << "(" << doubleToStr(lowVal) << " <= " << featName
       << " & " << featName << " <= " << doubleToStr(highVal) << ")";
    return ss.str();
  } else if(hasLower){
    std::ostringstream ss;
    ss << "(" << featName << " >= " << doubleToStr(lowVal) << ")";
    return ss.str();
  } else if(hasUpper){
    std::ostringstream ss;
    ss << "(" << featName << " <= " << doubleToStr(highVal) << ")";
    return ss.str();
  } else {
    // no constraints => always true in that dimension
    return "";
  }
}

// Combine dimension constraints into final rule string
static std::string boxToRuleString(const RegionBox &box,
                                   const CharacterVector &featNames)
{
  int p = featNames.size();
  std::ostringstream oss;
  bool first = false;
  for(int d=0; d<p; d++){
    std::string fn = Rcpp::as<std::string>(featNames[d]);
    std::string expr = makeOneDimConstraint(fn, box.lower[d], box.upper[d]);
    if(!expr.empty()){
      if(first) oss << " & ";
      oss << expr;
      first = true;
    }
  }
  // If no constraints, return "TRUE"
  if(!first) return std::string("TRUE");
  return oss.str();
}


/* --------------------------------------------------------------------------
 6) Main function with beam search for up to 4D bounding boxes,
 now with excludedBoxes
 -------------------------------------------------------------------------- */

/**
 * @param excludedBoxes R list, each element must be a list with:
 *        - lowerBound (NumericVector)
 *        - upperBound (NumericVector)
 *   We will parse them into std::vector<RegionBox> and skip them.
 */
// [[Rcpp::export]]
List find_best_region_beam_cpp(const NumericMatrix &X,
                               const NumericVector &residuals,
                               const LogicalVector &assignedMask,
                               const NumericMatrix &thresholds_list,
                               double min_obs_pct = 0.05,
                               double max_obs_frac = 1.0,
                               int K1 = 200,   // top-K for 1D
                               int K2 = 200,   // top-K for 2D
                               int K3 = 200,   // top-K for 3D
                               int K4 = 200,   // top-K for 4D
                               Rcpp::CharacterVector featureNames = CharacterVector(),
                               SEXP excludedBoxes = R_NilValue  // optional argument
)
{
  int n = X.nrow();
  int p = X.ncol();
  
  //---------------------------------------------
  // 1) Convert excludedBoxes from R into C++ vector
  //---------------------------------------------
  std::vector<RegionBox> excludedVec;  
  if(!Rf_isNull(excludedBoxes)) {
    // Expecting: a list of region-lists, each has lowerBound, upperBound
    List eboxes(excludedBoxes);
    for(int i=0; i < eboxes.size(); i++){
      List one = eboxes[i];
      NumericVector lb = one["lowerBound"];
      NumericVector ub = one["upperBound"];
      if(lb.size() != p || ub.size() != p){
        Rcpp::stop("excludedBoxes: dimension mismatch");
      }
      RegionBox rb(p);
      for(int d=0; d<p; d++){
        rb.lower[d] = lb[d];
        rb.upper[d] = ub[d];
      }
      excludedVec.push_back(rb);
    }
  }
  
  // Count how many points are unassigned
  int unassignedCount = 0;
  for(int i=0; i<n; i++){
    if(!assignedMask[i]) {
      unassignedCount++;
    }
  }
  if(unassignedCount == 0){
    // No unassigned points => no valid region
    return List::create(_["found_region"] = false);
  }
  
  // coverage thresholds
  int minCount = std::max((int)std::floor(unassignedCount * min_obs_pct), 1);
  int maxCount = (int)std::floor((double)unassignedCount * max_obs_frac);
  if(maxCount < minCount){
    return List::create(_["found_region"] = false);
  }
  
  // Precompute thresholds for each dimension
  std::vector< std::vector<double> > allThresh(p);
  for(int j=0; j<p; j++){
    allThresh[j] = collectThresholds(thresholds_list, j);
  }
  
  //---------------------------------------------
  // 2) Build candidate sets for 1D/2D/3D/4D
  //---------------------------------------------
  
  // Step A: single-dimensional expansions
  std::vector<CandBox> beam1;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<CandBox> localVec;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(int j=0; j<p; j++){
    const auto &tj = allThresh[j];
    if(tj.empty()) continue;
    int szj = (int)tj.size();
    for(int a=0; a<szj; a++){
      for(int b=a; b<szj; b++){
        RegionBox cand(p);
        cand.lower[j] = tj[a];
        cand.upper[j] = tj[b];
        
        // skip if in excludedVec
        if(isExcluded(cand, excludedVec)) {
          continue;
        }
        
        CandBox cbox;
        evaluateCandidate(cand, X, residuals, assignedMask,
                          minCount, maxCount, cbox);
        
        if(cbox.coverage >= minCount && cbox.coverage <= maxCount && cbox.metric > 0){
          localVec.push_back(std::move(cbox));
        }
      }
    }
  }
#ifdef _OPENMP
#pragma omp critical
#endif
{
  beam1.insert(beam1.end(), localVec.begin(), localVec.end());
}
} // end parallel
beam1 = keepTopK(beam1, K1);

// Step B: Expand 1D -> 2D
std::vector<CandBox> beam2;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<CandBox> localVec;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(int i=0; i<(int)beam1.size(); i++){
    const RegionBox &box1 = beam1[i].box;
    // Try adding a second dimension
    for(int j2=0; j2<p; j2++){
      // If box1 already restricts j2, skip
      if(box1.lower[j2] != -std::numeric_limits<double>::infinity() ||
         box1.upper[j2] !=  std::numeric_limits<double>::infinity())
      {
        continue;
      }
      const auto &t2 = allThresh[j2];
      if(t2.empty()) continue;
      int sz2 = (int)t2.size();
      for(int c=0; c<sz2; c++){
        for(int d=c; d<sz2; d++){
          RegionBox cand = box1;
          cand.lower[j2] = t2[c];
          cand.upper[j2] = t2[d];
          
          if(isExcluded(cand, excludedVec)) {
            continue;
          }
          
          CandBox cbox;
          evaluateCandidate(cand, X, residuals, assignedMask,
                            minCount, maxCount, cbox);
          
          if(cbox.coverage >= minCount && cbox.coverage <= maxCount && cbox.metric > 0){
            cbox.box = cand;
            localVec.push_back(std::move(cbox));
          }
        }
      }
    }
  }
#ifdef _OPENMP
#pragma omp critical
#endif
{
  beam2.insert(beam2.end(), localVec.begin(), localVec.end());
}
} // end parallel
beam2 = keepTopK(beam2, K2);

// Step C: Expand 2D -> 3D
std::vector<CandBox> beam3;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<CandBox> localVec;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(int i=0; i<(int)beam2.size(); i++){
    const RegionBox &box2 = beam2[i].box;
    // find a third dimension j3 not used
    for(int j3=0; j3<p; j3++){
      if(box2.lower[j3] != -std::numeric_limits<double>::infinity() ||
         box2.upper[j3] !=  std::numeric_limits<double>::infinity())
      {
        continue;
      }
      const auto &t3 = allThresh[j3];
      if(t3.empty()) continue;
      int sz3 = (int)t3.size();
      for(int e=0; e<sz3; e++){
        for(int f=e; f<sz3; f++){
          RegionBox cand = box2;
          cand.lower[j3] = t3[e];
          cand.upper[j3] = t3[f];
          
          if(isExcluded(cand, excludedVec)) {
            continue;
          }
          
          CandBox cbox;
          evaluateCandidate(cand, X, residuals, assignedMask,
                            minCount, maxCount, cbox);
          
          if(cbox.coverage >= minCount && cbox.coverage <= maxCount && cbox.metric > 0){
            cbox.box = cand;
            localVec.push_back(std::move(cbox));
          }
        }
      }
    }
  }
#ifdef _OPENMP
#pragma omp critical
#endif
{
  beam3.insert(beam3.end(), localVec.begin(), localVec.end());
}
} // end parallel
beam3 = keepTopK(beam3, K3);

// Step D: Expand 3D -> 4D
std::vector<CandBox> beam4;
#ifdef _OPENMP
#pragma omp parallel
#endif
{
  std::vector<CandBox> localVec;
#ifdef _OPENMP
#pragma omp for schedule(dynamic)
#endif
  for(int i=0; i<(int)beam3.size(); i++){
    const RegionBox &box3 = beam3[i].box;
    // find a fourth dimension j4 not used
    for(int j4=0; j4<p; j4++){
      if(box3.lower[j4] != -std::numeric_limits<double>::infinity() ||
         box3.upper[j4] !=  std::numeric_limits<double>::infinity())
      {
        continue;
      }
      const auto &t4 = allThresh[j4];
      if(t4.empty()) continue;
      int sz4 = (int)t4.size();
      
      for(int e=0; e<sz4; e++){
        for(int f=e; f<sz4; f++){
          RegionBox cand = box3;  // copy
          cand.lower[j4] = t4[e];
          cand.upper[j4] = t4[f];
          
          if(isExcluded(cand, excludedVec)) {
            continue;
          }
          
          CandBox cbox;
          evaluateCandidate(cand, X, residuals, assignedMask,
                            minCount, maxCount, cbox);
          
          if(cbox.coverage >= minCount && cbox.coverage <= maxCount && cbox.metric > 0){
            cbox.box = cand;
            localVec.push_back(std::move(cbox));
          }
        }
      }
    }
  }
#ifdef _OPENMP
#pragma omp critical
#endif
{
  beam4.insert(beam4.end(), localVec.begin(), localVec.end());
}
} // end parallel
beam4 = keepTopK(beam4, K4);

// Merge all candidate lists
std::vector<CandBox> allCands;
allCands.reserve(beam1.size() + beam2.size() + beam3.size() + beam4.size());
allCands.insert(allCands.end(), beam1.begin(), beam1.end());
allCands.insert(allCands.end(), beam2.begin(), beam2.end());
allCands.insert(allCands.end(), beam3.begin(), beam3.end());
allCands.insert(allCands.end(), beam4.begin(), beam4.end());

// Pick the single best CandBox across all
CandBox globalBest;
for(const auto &cb : allCands){
  if(cb.metric > globalBest.metric){
    globalBest = cb;
  }
}
if(globalBest.metric < 0.0){
  // no valid region found
  return List::create(_["found_region"] = false);
}

// Build final output
NumericVector lb(p), ub(p);
for(int d=0; d<p; d++){
  lb[d] = globalBest.box.lower[d];
  ub[d] = globalBest.box.upper[d];
}
IntegerVector ptsInside(globalBest.inside.size());
for(int i=0; i<(int)globalBest.inside.size(); i++){
  ptsInside[i] = globalBest.inside[i];
}
std::string ruleStr = boxToRuleString(globalBest.box, featureNames);

return List::create(
  _["found_region"]  = true,
  _["lowerBound"]    = lb,
  _["upperBound"]    = ub,
  _["points_inside"] = ptsInside,
  _["metric"]        = globalBest.metric,
  _["coverage"]      = globalBest.coverage,
  _["rule_string"]   = ruleStr
);
}
