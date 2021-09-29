#include <cassert>
#include <cctype>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <limits>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "htslib/faidx.h"
#include "htslib/hts.h"
#include "htslib/sam.h"

#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/binomial.hpp>

#include "../include/hmcnc.h"
#include "../include/hmcnc_defs.h"
#include "../include/hmcnc_io.h"

using boost::math::binomial;

using boost::math::poisson;
using boost::math::pdf;
using boost::math::negative_binomial_distribution;
using std::vector;
using std::cout;
using std::string;
using std::log;
using namespace std;

constexpr std::array<int8_t, 256> NucMap{
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 15
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 31
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 47
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 63

//  A   C        G
  4,0,4,1, 4,4,4,2, 4,4,4,4, 4,4,4,4,  // 79
//         T
  4,4,4,4, 3,4,4,4, 4,4,4,4, 4,4,4,4,  // 95
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 111
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 127

  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 143
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 159
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 175
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 191

  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 207
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 223
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 239
  4,4,4,4, 4,4,4,4, 4,4,4,4, 4,4,4,4,  // 255
};

struct BamHeaderDeleter {
  void operator()(bam_hdr_t* hdr) const noexcept {
    bam_hdr_destroy(hdr);
    hdr = nullptr;
  }
};

struct BamRecordDeleter {
  void operator()(bam1_t* b) const noexcept {
    bam_destroy1(b);
    b = nullptr;
  }
};

struct FastaIndexDeleter {
  void operator()(faidx_t* fai) const noexcept {
    fai_destroy(fai);
    fai = nullptr;
  }
};

struct HtslibFileDeleter {
  void operator()(htsFile* file) const noexcept {
    if (file) {
      hts_close(file);
    }
    file = nullptr;
  }
};

struct HtslibIndexDeleter {
  void operator()(hts_idx_t* index) const noexcept {
    hts_idx_destroy(index);
    index = nullptr;
  }
};

struct HtslibIteratorDeleter {
  void operator()(hts_itr_t* iter) const noexcept {
    hts_itr_destroy(iter);
    iter = nullptr;
  }
};

void Moments(const vector<double> &v, double &ex, double &var) {
  const int numElements = static_cast<int>(v.size());

  ex=0;
  for (int i=0; i < numElements; i++) {
    ex+=v[i]*i;
  }
  var=0;
  for (int i=0; i < numElements; i++) {
    var+= v[i]*(i-ex)*(i-ex);
  }
}

void Reset(vector<vector<double> > &v) {
  for (auto &e : v) {
    fill(e.begin(), e.end(), 0);
  }
}

const double oned = 1.0;
const double ONE = log(oned - 1e-10);
const double epsilon = 0.00000000001;

double PairSumOfLogP(double a, double b) {
  double res = b;
  if (a != 0) {
    const double m = max(a, b);

    assert(a <= 0);
    const double diff = min(a,b) - m;
    const double e = exp(diff);
    const double lg=log(1+e);
    res = m + lg;
    assert(res < epsilon);
    if (res > 0 and res < epsilon ) {
      return 0;
    }
  }
  if (res > minNeg) {
    res=minNeg;
  }

  return res;
}

double SumOfLogP(const vector<double> &vals) {
  if (vals.empty()) {
    // Return 0 for error.
    return 0;
  }
  double maxVal = vals[0];
  for (size_t i=1; i < vals.size(); ++i) {
    maxVal=max(maxVal,vals[i]);
  }
  double expSum=0;
  for (const auto& v : vals) {
    expSum += exp(v - maxVal);
  }
  return maxVal+log(expSum);
}

using namespace std;
class ThreadInfo {
public:
  std::unique_ptr<htsFile, HtslibFileDeleter> htsfp;
  std::shared_ptr<hts_idx_t> bamidx;
  std::shared_ptr<bam_hdr_t> samHeader;
  std::unique_ptr<faidx_t, FastaIndexDeleter> fai;
  int *lastSeq;
  string refFileName;
  pthread_mutex_t *semaphore;
  vector<string> *contigNames;
  vector<int>    *contigLengths;
  vector<int> procChroms;
  vector<vector<int>> *covBins;
  vector<vector<int>> *clipBins;
  vector<vector<SNV>> *snvs;
  vector<vector<int>> *copyNumber;
  vector<vector<double>> *transP, *emisP, *expTransP, *expEmisP;
  vector<double> *startP;
  vector<vector<Interval>> *copyIntervals;
  int maxCov, maxState;
  bool exit;
  double mean;
  double var;
  double lepsi;
  double scale;
  double *pModel;
  vector<int> *totalReads;
  vector<long> *totalBases;
};

static void printModel(const vector<vector<double>> &transP, ostream *strm)
{
  *strm << "\nTRANS: \n";
  for (size_t r=0; r<transP.size(); r++) {
    *strm << r <<":";
    for (size_t c=0;c<transP[r].size();c++) {
      *strm << '\t' << transP[r][c];
    }
    *strm << '\n';
  }
  *strm << '\n';
}

static void printEmissions(const vector<vector<double> > &emisP, ostream *strm) {
  *strm << "\nEMIS: \n";
  assert(!emisP.empty());
  for (size_t i=0; i < emisP[0].size(); i++) {
    *strm << std::setw(7) << i;
  }
  *strm << '\n';

  for (const auto &e : emisP) {
    for (size_t c=0;c<e.size(); c++ ) {
      *strm << std::setw(7) << std::setprecision(2) << e[c];
      if (c+1 < e.size()) {
        *strm << '\t';
      }
    }
    *strm << '\n';
  }
}

double max_over_row(const vector<vector<double>> &v, size_t col, size_t nStates) {

  double maxi=-1 * (std::numeric_limits<double>::max());
  for (const auto &val : v) {
    assert(col < val.size());
    maxi=std::max(val[col], maxi);
  }
  return maxi;
}

double max_over_rows(const vector<vector<double>> &v, size_t col,
                     const vector<vector<double>> &v2, size_t nextState,
                     size_t nStates ) {
  double maxi2=-1 * (std::numeric_limits<double>::max());
  for(size_t i=0;i< nStates;i++) {
    maxi2=std::max(v[i][col] + v2[i][nextState], maxi2);
  }
  return maxi2;
}

double LgNegBinom(int cn, int cov, float Hmean, float Hvar) {
  double result=0;
  float r, p;

  if (Hmean == 0 or Hvar == 0) {
    return 0;
  }
  //
  // Since actual variance is unknown at high copy number states, assume linear increase.
  //

  if (Hmean==0) { //no alignment in contig
    if (cn!= 0) {
      result=lepsi;
    }
    else {
      result=0;
    }
  }
  else if (cov >= MAX_CN*Hmean) {//max_obs filtered previously
    if (cn!= MAX_CN) {
      result=lepsi;
    }
    else {
      result=0;
    }
  }
  else if(cn==0) {//del_states
    //
    // Use poisson for 0-state assuming it's a random mismap process.
    const poisson distribution(MISMAP_RATE*Hmean);
    const double prob=pdf(distribution, cov);
    if (prob == 0) {
      result=lepsi;
    }
    else {
      result=max(lepsi, log(prob));
    }
  }
  else {
    Hmean*=cn;
    Hvar*=cn;
    p=Hmean/Hvar;
    r=Hmean*(Hmean/Hvar)/(1-Hmean/Hvar);

    const negative_binomial_distribution<double> distribution(r,p);
    const double prob=pdf(distribution, cov);
    if (prob == 0) {
      result=lepsi;
    }
    else {
      result=max(lepsi, log(prob));
    }
  }
  return result;
}

double LgBinom(double p, int s, int n) {
  const binomial snv(n, p);
  const double pVal=pdf(snv,s);
  if (pVal == 0) {
    return lepsi;
  }
  else {
    const double lP=log(pVal);
    const double retVal=max(lP, lepsi);
    assert(isnan(lP) == 0);
    return retVal;
  }
}

double LgPrpoiss(int cn,  int cov, int Hmean) {
  double result=0;

  // Boundary case
  if (Hmean==0) {//no alignment in contig
    if(cn!= 0) {
      // Cannot emit 0 reads from non-zero states.
      result=lepsi;
    }
    else {
      // Can only emit 0
      result=0;
    }
  }
  else if (cov >= MAX_CN*Hmean) {//max_obs filtered previously
    //
    // Have state that catches ultra-high probability states.
    if(cn != MAX_CN) {
      result=lepsi;
    }
    else {
      result=0;
    }
  }
  else if(cn==0) {//del_states
    const poisson distribution(MISMAP_RATE*Hmean);
    const double prob=pdf(distribution, cov);
    if (prob == 0) {
      result = lepsi;
    }
    else {
      result = max(lepsi, log(prob));
    }
  }
  else {
    const poisson distribution(cn*Hmean);
    const double prob=pdf(distribution, cov);
    if (prob == 0) {
      result=lepsi;
    }
    else {
      result=max(lepsi, log(prob));
    }
  }
  return result;
}

static void correctModel(vector<vector<double>> &transP,
                         int nStates)
{
  double sum;
  assert(nStates < transP.size());
  for (int i=0;i<nStates;i++) {
    sum = 0;
    assert(nStates < transP[i].size());
    for (int j=0;j<nStates;j++) {
      sum+= std::exp(transP[i][j]);
    }
    for (int j=0;j<nStates;j++) {
      transP[i][j]= log(std::exp(transP[i][j])/sum);
    }
  }
}//correctModel

void CombineEmissions(const vector<int> &obs,
                      const vector<SNV> &snvs,
                      vector<uint8_t> &isCov,
                      vector<int> &obsIndex) {
  const int totObs=obs.size() + snvs.size();
  isCov.resize(totObs);
  fill(isCov.begin(), isCov.end(), false);
  obsIndex.resize(totObs);
  fill(obsIndex.begin(), obsIndex.end(), -1);
  size_t c=0,s=0,p=0, obsi=0;
  while (c < obs.size() or s < snvs.size()) {
    const int curPos=c*100;
    while (s < snvs.size() and snvs[s].pos < curPos) {
      isCov[obsi] = false;
      obsIndex[obsi] = s;
      obsi++;
      s++;
    }
    isCov[obsi] = true;
    obsIndex[obsi] = c;
    obsi++;
    c++;
  }
  if (c != obs.size() or s != snvs.size()) {
    cerr << "ERROR computing obs incex." << '\n';
    exit(1);
  }
}

double CSEmisP(int state,
               int pos,
               const vector<int> &obs,
               const vector<SNV> &snvs,
               const vector<uint8_t> &isCov,
               const vector<int> &obsIndex,
               const vector<vector<double>> &emisP,
               const vector<vector<vector<double>>> &binoP ) {
  if (isCov[pos]) {
    const int covIndex=obsIndex[pos];
    return emisP[state][obs[covIndex]];
  }
  else {
    const int snvIndex = obsIndex[pos];
    const int maxCov = emisP[0].size();
    int ref=snvs[snvIndex].ref;
    int alt=snvs[snvIndex].alt;
    int nucCov = ref+alt;
    if (nucCov > maxCov) {
      const double f=((double)maxCov)/nucCov;
      ref=ref*f;
      alt=alt*f;
      nucCov=ref+alt;
    }
    assert(nucCov < binoP[state].size());
    const int m=min(ref,alt);
    assert(m < binoP[state][nucCov].size());
    return binoP[state][nucCov][m];
  }
}

double ForwardBackwards(const vector<double> &startP,
                        const vector<vector<double>> &covCovTransP,
                        const vector<vector<double>> &emisP,
                        const vector<int> &obs,
                        vector<vector<double>> &f,
                        vector<vector<double>> &b) {
  const int totObs    = static_cast<int>(obs.size());
  const int nCovObs   = static_cast<int>(obs.size());
  const int nCovStates = static_cast<int>(startP.size());
  assert(nCovStates > 0);
  //
  // Eventually this can be done with logic in the code, but for now just store a
  // flag at each observation if it is an SNV or cov emission.
  //

  //
  // Initialize first col from standing.
  //

  // f and b are of length nObs+1
  // The first/last prob from an observation are at:
  //  f[1], f[nObs], from obs[0...nObs-1]
  //  and b[nObs], b[1] from obs[nObs-1 ... 0]
  //
  f.resize(nCovStates);
  b.resize(nCovStates);

  for (int i = 0; i < nCovStates; i++) {
    f[i].resize(totObs+1);
    fill(f[i].begin(), f[i].end(), 0);
    b[i].resize(totObs+1);
    fill(b[i].begin(), b[i].end(), 0);
  }

  for (int j=0; j < nCovStates; j++) {
    f[j][0] = log(1./nCovStates);
  }

  const double lgthird=log(1/3.);
  //
  // Shouldn't have 0 states. Even if the coverage is empty for
  // an entire chrom, should have 0 state.
  //

  //
  // If just one state (e.g. all zeros), same prob. Set to -1 to flag for now.
  //
  if (nCovStates == 1) {
    fill(f[0].begin(), f[0].end(), -1);
    fill(b[0].begin(), b[0].end(), -1);
    return 0;
  }
  int prevCovIdx=0, curCovIdx=0;
  int prevSNVIdx=0, curSNVIdx=0;

  for (int k=0; k < totObs; k++) {
    for (int i=0; i < nCovStates; i++) {
      double colSum=0;
      for (int j=0; j < nCovStates; j++) {
        assert(j== 0 or colSum != 0);
        assert(j < f.size());
        assert(k < f[j].size());
        assert(j < covCovTransP.size());
        assert(i < covCovTransP[j].size());
        colSum = PairSumOfLogP(colSum, f[j][k] + covCovTransP[j][i]);
      }
      assert(obs[k] < emisP[i].size());
      f[i][k+1] = colSum + emisP[i][obs[k]];
    }
  }

  // back
  for (int j=0; j < nCovStates; j++) {
    b[j][totObs] = startP[j];
  }

  for (int k = totObs-1; k > 0; k--) {
    for (int i=0; i < nCovStates; i++) {
      double colSum=0;
      for (int j=0; j < nCovStates; j++) {
        assert(j== 0 or colSum != 0);
        assert(prevCovIdx < obs.size()+1);
        colSum = PairSumOfLogP(colSum, b[j][k+1] + covCovTransP[j][i] + emisP[j][obs[k]]);
      }
      b[i][k] = colSum;
    }
  }

  double finalCol=0;
  for (int j=0; j < nCovStates; j++) {
    finalCol = PairSumOfLogP(finalCol, f[j][totObs]+log(1.0/nCovStates));
  }
  return finalCol;
}

double BaumWelchEOnChrom(const vector<double> &startP,
			 vector<vector<double>> &covCovTransP,
			 vector<vector<double>> &emisP,
			 vector<int> &obs,
			 vector<vector<double>> &f,
			 vector<vector<double>> &b,
			 vector<vector<double>> &expCovCovTransP,
			 vector<vector<double>> &expEmisP) {

  const int nStates = static_cast<int>(startP.size());
  const int nObs = obs.size();
  const double px = ForwardBackwards(startP, covCovTransP, emisP, obs, f, b);
  /*
  ofstream fb("fb.tsv");
  for (int k=0; k < f[0].size()-1; k++) {
    fb << k << '\t';
    double maxfb=f[0][k]+b[0][k+1];
    int maxi=0;
    double fbSum=0;
    for (int j=0; j < nStates; j++) {
      double p=f[j][k] + b[j][k+1];
      fbSum=PairSumOfLogP(fbSum, p);
    }
    for (int j=0; j < nStates; j++) {
      double p=f[j][k] + b[j][k+1];
      fb << std::setprecision(8) << f[j][k] << ", " << b[j][k+1] << ", " << f[j][k] + b[j][k+1];
	//fCov[j][k] << ", " << bCov[j][k+1] << ", " << p-fbSum;
      if (j +1 < nStates) { fb << '\t';}
      if (maxfb < p) {
	maxfb=p;
	maxi=j;
      }
    }
    fb << '\t' << maxi << '\t' << maxfb/fbSum << '\t' << obs[k] << '\n';

  }
  fb.close();
*/

  for (int k=1; k< nObs-1; k++) {
    double logSum=0;
    //
    // Calculate total probability of all transitions at this step.
    //
    for (size_t i=0; i < covCovTransP.size(); i++) {
      covCovTransP[i].resize(covCovTransP[i].size());
      for (size_t j=0; j < covCovTransP[i].size(); j++) {
        logSum = PairSumOfLogP(logSum, f[i][k] + covCovTransP[i][j] + emisP[j][obs[k+1]] + b[j][k+1]);
      }
    }
    for (size_t i=0; i < covCovTransP.size(); i++) {
      for (size_t j=0; j < covCovTransP[i].size(); j++) {
        const double pEdge=f[i][k] + covCovTransP[i][j] + emisP[j][obs[k+1]] + b[j][k+1];
        assert(isnan(exp(pEdge-logSum)) == false);
        expCovCovTransP[i][j] += exp(pEdge - logSum);
      }
    }
  }
  for (size_t ei=0; ei < expEmisP.size(); ei++) {
    fill(expEmisP[ei].begin(), expEmisP[ei].end(), 0);
  }
  for (int k=1; k< nObs-1; k++) {
    double colSum=0;

    for (int j=0; j < nStates; j++) {
      colSum=PairSumOfLogP(colSum, f[j][k] + b[j][k]);
    }
    for (int j=0; j < nStates; j++) {
      expEmisP[j][obs[k+1]] += exp(f[j][k]+b[j][k] - colSum);
    }
  }
  for (size_t ei=0; ei < expEmisP.size(); ei++) {
    double rowSum=0;
    rowSum=std::accumulate(expEmisP[ei].begin(),expEmisP[ei].end(), rowSum);
    for (int pi=0; pi < expEmisP[ei].size(); pi++) {
      expEmisP[ei][pi] /= rowSum;
    }
  }

  return px;
}

void PrintIntervals(const string chrom, const vector<Interval> &intv, ostream &out) {
  for (const auto &i : intv) {
    out << chrom << '\t'
        << i.start << '\t'
        << i.end << '\t'
        << i.copyNumber << '\t'
        << i.averageCoverage << '\t'
        << i.pVal << '\n';
  }
}

void StorePosteriorMaxIntervals(const vector<int> &cov,
				const vector<vector<double>> &f,
				const vector<vector<double>> &b,
				vector<Interval> &intervals) {
  intervals.resize(0);
  int i=0;
  if (f.size() == 0) {
    return;
  }
  int prevCN=-1;
  int prevStart=0;
  int totCov=0;
  int curCov=0;
  int nSNV=0;
  double colSum;
  double avgPVal=0.0;
  for (i=0; i < static_cast<int>(f[0].size()-1); i++) {
    double colMax=f[0][i] + b[0][i];
    int colMaxCN=0;
    assert(i < cov.size());
    totCov+=cov[i];
    colSum=f[0][i] + b[0][i];
    for (size_t j = 1; j < f.size(); j++) {
      const double v=f[j][i] + b[j][i];
      if (v > colMax) {
	      colMax=v;
	      colMaxCN=j;
      }
      colSum=PairSumOfLogP(colSum, f[j][i] + b[j][i]);
    }
    avgPVal+=colMax-colSum;
    if (prevCN==-1) {
      prevCN=colMaxCN;
      prevStart=i;
      continue;
    }
    if (colMaxCN != prevCN) {
      intervals.push_back(Interval(prevStart*100, i*100, prevCN, totCov/(i-prevStart), avgPVal/(i-prevStart)));
      //      out << chrom << '\t' << prevStart*100 << '\t' << i*100 << '\t' << prevCN << '\t' << i-prevStart << '\t' << totCov/(i-prevStart) << '\t' << nSNV << '\n';
      prevCN=colMaxCN;
      prevStart=i;
      totCov=0;
      nSNV=0;
    }
  }
  if (i - prevStart > 0) {
    intervals.push_back(Interval(prevStart*100, i*100, prevCN, totCov/(i-prevStart), avgPVal/(i-prevStart)));
  }
}

void UpdateEmisP(vector<vector<double>> &emisP,
                 vector<vector<double>> &expEmisP,
                 int model) {
  const int nCovStates=emisP.size();
  assert(nCovStates <= expEmisP.size());
  for (int i=1;i<nCovStates;i++) {
    double mean, var;
    Moments(expEmisP[i], mean, var);
    const int maxCov=expEmisP[i].size()-1;
    double stateSum=0;
    emisP[i].resize(expEmisP[i].size());
    for (int j=0;j<=maxCov;j++) {
      if (model == POIS or mean > var) {
        emisP[i][j]=LgPrpoiss( (int) i , j , (int) mean );
        stateSum+=exp(LgPrpoiss( (int) i , j , (int) mean ));
      }
      else {
        emisP[i][j]=LgNegBinom((int)i, (int) j, mean, var);
        stateSum+=exp(LgNegBinom((int)i, (int) j, mean, var));
      }
    }
  }
}

void BaumWelchM(const vector<double> &startP,
		const vector<vector<double>> &transP,
		const vector<vector<double>> &emisP,
		const vector<vector<vector<double>>> &binoP,
		int model,
		const vector<long> &stateTotCov,
		const vector<long> &stateNCov,
		const vector<vector<double>> &expTransP,
		vector<vector<double>> &expEmisP,
		vector<vector<double>> &covCovPrior,
		vector<vector<double>> &updateTransP,
		vector<vector<double>> &updateEmisP) {

  //
  // M step.
  //
  const int nStates=static_cast<int>(startP.size());
  updateTransP.resize(nStates);
  vector<double> colSums;

  cerr << "Update trans: " << '\n';
  cerr << "p\t";
  for (int j=0; j < nStates; j++) {
    cerr << std::setw(8) << j << '\t';
  }
  cerr << '\n';
  for (int j=0; j < nStates; j++) {
    double colSum=0;
    updateTransP[j].resize(nStates);
    assert(nStates <= expTransP[j].size());
    for (int k=0; k< nStates; k++) {
      colSum +=expTransP[j][k]; //PairSumOfLogP(colSum, expTransP[j][k]);
    }
    colSums.push_back(colSum);
    cerr << j;
    for (int k=0; k < nStates; k++) {
      updateTransP[j][k] = log(expTransP[j][k]/colSum); //min(ONE, expTransP[j][k] - colSum);
      cerr << '\t' << std::setw(8) << updateTransP[j][k];
    }
    cerr << '\n';
  }
  //
  // M step for emissions -- use summary statistics to recompute emission values
  // Eventually to really be BW, this should fit a negative binomial dist to the data.
  //
  updateEmisP.resize(nStates);

  vector<double> stateMean(nStates);
  assert(nStates <= stateNCov.size());
  assert(nStates <= stateTotCov.size());
  for (int i=0; i < nStates; i++) {
    if (stateNCov[i] > 0) {
      stateMean[i]=stateTotCov[i]/stateNCov[i];
    }
    else {
      stateMean[i] = 0;
    }
  }
  assert(!emisP.empty());
  const int maxCov=static_cast<int>(emisP[0].size());
  UpdateEmisP(updateEmisP, expEmisP, model);
}

void viterbi(const vector<double> &startP,
             const vector<vector<double>> &covCovTransP,
             const vector<vector<double>> &covSnvTransP,
             const vector<vector<double>> &snvSnvTransP,
             const vector<vector<double>> &emisP,
             const vector<vector<vector<double>>> &binoP,
             const vector<int> &cov,
             const vector<SNV> &snvs,
             const vector<uint8_t> &isCov,
             const vector<int> &obsIndex,
             vector<int> &viterbiPath) {
  /*
  //size_t  nObservations  = observations.size();
  size_t nObs=obsIndex.size();

  int nCovStates=covCovTransP.size();
  int nSNVStates=snvSnvTransP.size() + 1;
  vector<vector<double> > v(nStates, vector<double>(nObs) );
  //  vector<vector<double> > f(nStates, vector<double>(observations.size()) );
  //  vector<vector<double> > b(nStates, vector<double>(observations.size()) );
  vector<vector<double> > opt(nStates, vector<double>(nObs)  );
  // Init


  for(size_t i=0;i<nCovStates;i++) {
    v[i][0] = startP[i] + emisP[i][cov[0]];
  }
  // Iteration
  vector<double> colProb(nCovStates);
  for(size_t k=1 ; k<nObs ; k++) {
    if (isCov[k] and isCov[k-1]) {
      int curCovIdx=obsIndex[k];
      int prevCovIdx=obsIndex[k-1];
      for(size_t i=0;i<nCovStates;i++) {
	double maxProb = v[0][k-1] + covCovTransP[0][i];
	int maxState=0;
	for(size_t j=1; j < nCovStates; j++) {
	  double rowProb = v[j][k-1] + covCovTransP[j][i];
	  if (rowProb > maxProb) {
	    maxState=j;
	    maxProb=rowProb;
	  }
	}
	v[i][k] = maxProb + emisP[i][obs[curCovIdx]];
	opt[i][k] = maxState;
      }
    }
    else if (isCov[k] and isCov[k-1] == false) {
      //
      // Transition from snv to cov state. This one is pretty complicated.
      // For states that are between 1 and 3, the copy number does not change.
      //
      int prevSnvIdx    = obsIndex[k-1];
      int curCovIdx     = obsIndex[k];
      double maxSNVProb = v[1][k-1];
      int maxSNVIndex   = 1;
      for (int i=2; i <= 3; i++ ){
	if (maxSNVProb < v[i][k-1]) {
	  maxSNVProb = v[i][k-1];
	  maxSNVState = i;
	}
      }
      for (int i=0; i < nCovStates; i++ ) {
	double maxProb = v[0][k-1] + covSnvTransP[0][i];
	maxState = 0;
	if (i >= 1 and i <= 3) {
	  v[i][k] = v[i][k-1];
	  opt[i][k] = i;
	}
	else {
	  v[i][k] = maxSNVProb;
	  opt[i][k] = maxSNVState;
	}
      }
      else if (isCov[k] == false and isCov[k-1] == true) {
	//
	// COV->SNV state.
	//
	// Not done, assert here.
	assert(0);


      }
    }
  }


  // Traceback
  for(size_t i=0;i<nStates;i++)
    {
      if( max_over_row(v,nObs-1,nStates) == v[i][nObs-1] )
        {
	  viterbiPath[nObs-1] = i;
	  break;
        }
    }
  size_t lastObservation=nObs-2;
  //
  // Find highest-scoring final entry.
  //
  size_t rowIndex=nObs-2;
  double maxScore=v[0][rowIndex];
  int    maxState=0;
  for (size_t i=1;i < nStates; i++) {
    if (maxScore < v[i][rowIndex]) {
      maxState=i;
      maxScore=v[i][rowIndex];
    }
  }
  viterbiPath[lastObservation] = maxState;


  for( size_t f=lastObservation; f > 0; f--) {
      viterbiPath[f] = maxState;
      maxState=opt[maxState][f];
    }
  */
}//viterbi

class CountNuc {
public:
  int index;
  int count;

  int operator<(const CountNuc &rhs) const {
    return count < rhs.count;
  }
};

int StoreSNVs(char *contigSeq, int contigLength, float mean,
              vector<int> &nA, vector<int> &nC,
              vector<int> &nG, vector<int> &nT,
              vector<int> &nDel, vector<SNV> &snvs) {

  //
  // Easier for typing
  //
  const char *nucs="ACGTd";

  assert(contigLength >= 0);
  assert(nA.size() <= contigLength);
  assert(nC.size() <= contigLength);
  assert(nG.size() <= contigLength);
  assert(nT.size() <= contigLength);
  assert(nDel.size() <= contigLength);

  vector<int*> fPtr(5);
  fPtr[0] = &nA[0];
  fPtr[1] = &nC[0];
  fPtr[2] = &nG[0];
  fPtr[3] = &nT[0];
  fPtr[4] = &nDel[0];

  vector<CountNuc> counts;
  counts.resize(5);
  long total=0;
  if (mean==0) {
    for (int i=0; i < contigLength; i++) {
      for (size_t j=0; j < fPtr.size(); j++) {
        total+=fPtr[j][i];
      }
    }
    mean=((double)total)/contigLength;
  }
  for (int i =0; i < contigLength; i++) {
    int tot=0;
    for (int n=0; n < 5; n++) {
      counts[n].index=n;
      counts[n].count = fPtr[n][i];
      tot+=fPtr[n][i];
    }
    if (tot == 0) {
      continue;
    }
    sort(counts.begin(), counts.end());
    const char refNuc=toupper(contigSeq[i]);
    if (counts[4].index != 4 and
        counts[3].index != 4 and
        refNuc != 'N' and
        counts[3].count > 0.25*mean and
        counts[4].count > 0.25*mean and
        counts[3].count < 2*mean ) {
      //
      // Top most frequent are not a deletion.
      //
      if (nucs[counts[4].index] == refNuc) {
        snvs.push_back(SNV(i, refNuc, nucs[counts[3].index], counts[4].count, counts[3].count));
      }
      else if (nucs[counts[3].index] == refNuc) {
        snvs.push_back(SNV(i, refNuc, nucs[counts[4].index], counts[4].count, counts[3].count));
      }
      if (snvs.size() % 100000 == 0) {
        cerr << "Stored " << snvs.size() << " at " << i << '\n';
      }
    }
  }
  return 1;
}

int IncrementCounts(bam1_t *b, int contigLength,
                    vector<int> &nA,
                    vector<int> &nC,
                    vector<int> &nG,
                    vector<int> &nT,
                    vector<int> &nDel) {
  const int readLength = b->core.l_qseq;
  if (readLength < BIN_LENGTH or b->core.qual < 10 or b->core.flag & 0x800) {
    return 0;
  }

  vector<char> seq(readLength);
  uint8_t *q = bam_get_seq(b);
  for (int i=0; i < readLength; i++) {
    seq[i]=seq_nt16_str[bam_seqi(q,i)];
  }
  uint32_t* cigar = bam_get_cigar(b);
  const int refLen = bam_cigar2rlen(b->core.n_cigar, cigar);
  int qPos=0;
  int refPos = b->core.pos;
  int regionOffset=refPos;
  bool first=true;
  for (size_t ci=0; ci < b->core.n_cigar; ci++) {
    const int opLen=bam_cigar_oplen(cigar[ci]);
    const int op=bam_cigar_op(cigar[ci]);

    if (op == BAM_CSOFT_CLIP) {
      qPos += opLen;
      continue;
    }
    if (op == BAM_CINS) {
      qPos += opLen;
      continue;
    }
    if (op == BAM_CDEL) {
      const int stop=refPos+opLen;
      for (; refPos < stop and refPos < contigLength; refPos++) {
        nDel[regionOffset]+=1;
        regionOffset++;
      }
      continue;
    }
    if (op == BAM_CMATCH or op == BAM_CEQUAL or op == BAM_CDIFF) {
      if (refPos + opLen <= 0) {
        refPos += opLen;
        qPos += opLen;
        continue;
      }
      else {
        for (int p=0; p < opLen; p++) {
          if (refPos >= contigLength) {
            break;
          }
          if (refPos >= 1) {
            first=false;
            const char nuc=toupper(seq[qPos]);
            assert(regionOffset < nA.size());
            /*
            if (regionOffset == 20504515 or refPos == 20504515 ) {
              cout << regionOffset << '\t' << nuc << '\t' << nC[regionOffset] << '\t' << nT[regionOffset] << '\t' << bam_get_qname(b) << '\n';
              }*/
            if (nuc == 'A') { nA[regionOffset]++; }
            if (nuc == 'C') { nC[regionOffset]++; }
            if (nuc == 'G') { nG[regionOffset]++; }
            if (nuc == 'T') { nT[regionOffset]++; }
            regionOffset++;
          }
          refPos++;
          qPos++;
        }
      }
    }
  }
  return 1;
}

void ThreadedBWE(ThreadInfo *threadInfo) {
  double pChrom;
  while (*(threadInfo->lastSeq) < (*(*threadInfo).contigNames).size()) {

    pthread_mutex_lock(threadInfo->semaphore);

    const int curSeq = *((*threadInfo).lastSeq);
    *(threadInfo->lastSeq) = *(threadInfo->lastSeq) + 1;
    pthread_mutex_unlock(threadInfo->semaphore);

    if (curSeq >= threadInfo->contigNames->size()) {
      break;
    }
    vector<vector<double>> f, b, expCovCovTransP, expEmisP;
    expCovCovTransP.resize(threadInfo->transP->size());
    for (int i=0; i < expCovCovTransP.size(); i++) {
      expCovCovTransP[i].resize(expCovCovTransP.size(), 0);
    }
    expEmisP.resize(threadInfo->emisP->size());
    for (int i=0; i < expEmisP.size(); i++){
      expEmisP[i].resize((*(threadInfo->emisP))[i].size());
    }

    pChrom = BaumWelchEOnChrom(*threadInfo->startP,
			       *threadInfo->transP,
			       *threadInfo->emisP,
			       (*threadInfo->covBins)[curSeq],
			       f, b,
			       expCovCovTransP,
			       expEmisP);

    //
    // Update expected transitions
    //
    pthread_mutex_lock(threadInfo->semaphore);
    for (size_t i=0; i < threadInfo->transP->size(); i++) {
      for (size_t j=0; j < (*threadInfo->transP)[i].size(); j++) {
        (*threadInfo->expTransP)[i][j] += expCovCovTransP[i][j];
      }
    }
    for (size_t i=0; i < threadInfo->emisP->size(); i++) {
      for (size_t j=0; j < (*threadInfo->emisP)[i].size(); j++) {
        (*threadInfo->expEmisP)[i][j] += expEmisP[i][j];
      }
    }

    *threadInfo->pModel  += pChrom;
    pthread_mutex_unlock(threadInfo->semaphore);
    StorePosteriorMaxIntervals((*threadInfo->covBins)[curSeq],
			       f, b,
			       (*threadInfo->copyIntervals)[curSeq]);
  }
}

void ParseChrom(ThreadInfo *threadInfo) {

  while (*(threadInfo->lastSeq) < (*(*threadInfo).contigNames).size()) {
    //
    // Grab current chrom to process
    //
    pthread_mutex_lock(threadInfo->semaphore);

    const int curSeq = *((*threadInfo).lastSeq);
    *(threadInfo->lastSeq) = *(threadInfo->lastSeq) + 1;
    (*(*threadInfo).totalReads)[curSeq] = 0;
    (*(*threadInfo).totalBases)[curSeq] = 0;

    //
    // Deal with race condition by double checking curSeq;
    //
    if (curSeq >= threadInfo->contigNames->size()) {
      pthread_mutex_unlock(threadInfo->semaphore);
      break;
    }

    (*threadInfo).procChroms.push_back(curSeq);
    //    (*threadInfo).covBins->push_back(vector<int>());
    //    (*threadInfo).snvs->push_back(vector<SNV>());
    pthread_mutex_unlock(threadInfo->semaphore);

    const int contigLength=threadInfo->contigLengths->at(curSeq);

    vector<int> nA(contigLength, 0), nC(contigLength, 0), nT(contigLength, 0), nG(contigLength,0), nDel(contigLength, 0);

    stringstream regionStrm;
    regionStrm << (*(*threadInfo).contigNames)[curSeq];// << ":1-" << contigLength;

    const string region=regionStrm.str();

    std::unique_ptr<hts_itr_t, HtslibIteratorDeleter> regionIter(
      sam_itr_querys(threadInfo->bamidx.get(), threadInfo->samHeader.get(), region.c_str()));

    int chromLen;
    char *chromSeq = fai_fetch(threadInfo->fai.get(), region.c_str(), &chromLen);

    bool continueParsing=true;
    vector<std::unique_ptr<bam1_t, BamRecordDeleter>> reads; //(bam_init1());
    long totalSize=0;
    int chunkNumber=0;
    int totalReads=0;
    int endpos=0;
    int startpos=0;
    while (continueParsing) {
      int totalSize=0;
      reads.resize(0);
      pthread_mutex_lock(threadInfo->semaphore);
      int bufSize=0;
      while (bufSize < 100000000 and continueParsing) {
        std::unique_ptr<bam1_t, BamRecordDeleter> b(bam_init1());
        const int res=sam_itr_next(threadInfo->htsfp.get(), regionIter.get(), b.get());
        bufSize+= b->l_data;
        totalSize+= b->l_data;

        if (res < 0) { // or totalReads < 15000) {
          continueParsing = false;
          cerr << "Ending parsing of " << region << " with " << totalSize << " data and " << chunkNumber << " iterations." << '\n';
          break;
        }
        /*
        cout << "read " << bam_get_qname(b) << '\n';
        if (strcmp("m64043_200714_124814/138021154/ccs", bam_get_qname(b)) == 0) {
          cout << "poblem" << '\n';
          }*/
        endpos=bam_endpos(b.get());
        reads.push_back(std::move(b));
        ++totalReads;
      }
      cerr << "Reading " << (*threadInfo->contigNames)[curSeq] << ", chunk " << chunkNumber << ".\t" << reads.size() << "/" << totalReads << " reads/total" << '\n';
      ++chunkNumber;
      pthread_mutex_unlock(threadInfo->semaphore);

      for (auto& b : reads) {
        IncrementCounts(b.get(), contigLength, nA, nC, nG, nT, nDel);
        endpos=bam_endpos(b.get());
        startpos=b->core.pos;
        (*(*threadInfo).totalReads)[curSeq]++;
        (*(*threadInfo).totalBases)[curSeq]+= endpos-startpos;

        const int nCigar=b->core.n_cigar;
        uint32_t*cigar=bam_get_cigar(b.get());
        int frontClip=-1, backClip=-1;
        int frontClipLen=0;
        int backClipLen=0;
        if (nCigar > 1 and bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
          frontClip=0;
          frontClipLen=bam_cigar_oplen(cigar[0]);
        }
        else if (nCigar > 2 and bam_cigar_op(cigar[1]) == BAM_CSOFT_CLIP) {
          frontClip=1;
          frontClipLen=bam_cigar_oplen(cigar[1]);
        }
        if (frontClip > 0 and bam_cigar_op(cigar[nCigar-1]) == BAM_CSOFT_CLIP) {
          backClip=nCigar-1;
          backClipLen=bam_cigar_oplen(cigar[nCigar-1]);
        }
        else if (frontClip > 0 and nCigar > 2 and bam_cigar_op(cigar[nCigar-2]) == BAM_CSOFT_CLIP) {
          backClip=nCigar-2;
          backClipLen=bam_cigar_oplen(cigar[nCigar-2]);
        }
        if (frontClipLen > MIN_CLIP_LENGTH) {
          const int bin=startpos/BIN_LENGTH;
          (*threadInfo->clipBins)[curSeq][bin] += 1;
        }
        if (backClipLen > MIN_CLIP_LENGTH) {
          const int bin=endpos/BIN_LENGTH;
          (*threadInfo->clipBins)[curSeq][bin] += 1;
        }
        b.reset(nullptr);
      }
    }
    // Never compute in the last bin
    const int nBins=contigLength/BIN_LENGTH;

    for (int bin=0; bin < nBins; bin++) {
      const int start=bin*BIN_LENGTH;
      const int end=min((bin+1)*BIN_LENGTH, contigLength);
      int totCov=0;
      for (int bp=start; bp < end; bp++) {
        totCov+=nA[bp] + nC[bp] + nG[bp] + nT[bp] + nDel[bp];
      }
      (*threadInfo->covBins)[curSeq][bin] =totCov/BIN_LENGTH;
    }

    //
    // Detect SNVS
    //
    StoreSNVs(chromSeq, chromLen, threadInfo->mean,
	      nA, nC, nG, nT, nDel,
	      (*threadInfo->snvs)[curSeq]);
    cerr << "Stored " << (*threadInfo->snvs)[curSeq].size()
         << " snvs for " << (*threadInfo->contigNames)[curSeq] << '\n';

  }
  if (threadInfo->exit) {
    pthread_exit(NULL);
  }
}

int GetRefAndAlt(char refNuc, const vector<int> &counts, int &ref, int &alt) {
  vector<int> sCounts=counts;
  sort(sCounts.begin(), sCounts.end());

  assert(sCounts.size() >= 4);

  const int second = sCounts[2];
  const int first  = sCounts[3];
  int firstIndex=-1, secondIndex=-1;
  for (int i=0; i < 4; i++) {
    if (firstIndex == -1 and counts[i] == first) {
      firstIndex=i;
    }
    else if (secondIndex == -1 and counts[i] == second) {
      secondIndex=i;
    }
  }
  if (firstIndex == -1 or secondIndex == -1) {
    ref=0;
    alt=0;
    return 4;
  }
  const int refNucIndex=NucMap[refNuc];
  if (first == 0 or second/first < 0.2) {
    ref=0;
    alt=0;
    return 4;
  }
  if (firstIndex == refNucIndex) {
    ref=first;
    alt=second;
    return secondIndex;
  }
  else {
    ref=second;
    alt=first;
    return firstIndex;
  }
}

const char* nucs = "ACGTN";

static int pileup_blank(void *data, bam1_t *b) {
  return 0;
}

int EstimateCoverage(const string &bamFileName,
                     const vector<vector<int>> &allCovBins,
                     const vector<string> &chroms,
                     const vector<int> &lengths,
                     string &useChrom,
                     double &mean,
                     double &var) {
  size_t useChromIndex=0;
  if (useChrom == "") {
    int maxLen=0;
    assert(lengths.size() <= allCovBins.size());
    for (size_t i=0; i < lengths.size(); i++) {
      long totCov=0;
      for (size_t j=0; j < static_cast<int>(allCovBins[i].size()); j++ ) {
        totCov+=allCovBins[i][j];
      }
      if (totCov > 0 and lengths[i] > maxLen) {
        useChrom = chroms[i];
        maxLen=lengths[i];
        useChromIndex=i;
      }
    }
  }
  cerr << "Estimating coverage from " << useChrom << '\n';
  int contigLength=0;
  assert(chroms.size() <= lengths.size());
  for (size_t i=0; i < chroms.size(); i++) {
    if (chroms[i] == useChrom) {
      contigLength=lengths[i];
      useChromIndex=i;
      break;
    }
  }

  if (contigLength == 0) {
    cerr << "ERROR Could not estimate coverage." << '\n';
    exit(1);
  }
  if (allCovBins.size() > 0) {
    assert(allCovBins[useChromIndex].size() == lengths[useChromIndex]/100);
    const size_t lastBin=allCovBins[useChromIndex].size();
    if (lastBin == 0) {
      cerr << "ERROR. Could not estimate coverage using precomputed bins." << '\n';
      exit(1);
    }
    long totCov=0;
    assert(useChromIndex <= allCovBins.size());
    assert(lastBin <= allCovBins[useChromIndex].size());
    for (size_t binIndex=0; binIndex<lastBin; binIndex++) {
      totCov+=allCovBins[useChromIndex][binIndex];
    }
    mean=totCov/lastBin;
    //
    // Recompute summary stats using limited data
    int nSamples=0;
    totCov=0;
    long totCovSq=0;
    for (size_t binIndex=0; binIndex<lastBin; binIndex++) {
      if (allCovBins[useChromIndex][binIndex] > 0.25*mean and
          allCovBins[useChromIndex][binIndex] < 1.75 * mean) {
        totCov+=allCovBins[useChromIndex][binIndex];
        totCovSq+=allCovBins[useChromIndex][binIndex]*allCovBins[useChromIndex][binIndex];
        nSamples++;
      }
    }
    if (nSamples > 0) {
      mean=totCov/nSamples;
      var=totCovSq/nSamples-mean*mean;
      cerr << "Estimating coverage on precomputed bins " << nSamples << " ending at " << mean << '\t' << var << '\n';
    }
    else {
      cerr << "Could not estimate coverage using precomputed coverage on chrom " << useChrom << '\n';
      exit(1);
    }
  }
  else {
    std::unique_ptr<htsFile, HtslibFileDeleter> htsfp(hts_open(bamFileName.c_str(),"r"));
    vector<int> covBins;

    std::unique_ptr<hts_idx_t, HtslibIndexDeleter> bamidx(sam_index_load(htsfp.get(), bamFileName.c_str()));
    if (bamidx == nullptr) {
      cerr << "ERROR reading index" << '\n';
      exit(0);
    }

    const htsFormat *fmt = hts_get_format(htsfp.get());
    if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
      cerr << "Cannot determine format of input reads." << '\n';
      exit(1);
    }

    std::unique_ptr<bam_hdr_t, BamHeaderDeleter> samHeader(sam_hdr_read(htsfp.get()));

    std::unique_ptr<hts_itr_t, HtslibIteratorDeleter> regionIter(
      sam_itr_querys(bamidx.get(), samHeader.get(), useChrom.c_str()));

    vector<int> nA(contigLength, 0), nC(contigLength, 0), nT(contigLength, 0), nG(contigLength,0), nDel(contigLength, 0);


    bool continueParsing=true;
    int nSamples=0;
    int curEndPos=0;
    int curCovBin=0;
    long totalSize;

    while (continueParsing) {
      int bufSize=0;
      int nReads=0;
      while (bufSize < 100000000 and continueParsing) {
        std::unique_ptr<bam1_t, BamRecordDeleter> b(bam_init1());
        const int res=sam_itr_next(htsfp.get(), regionIter.get(), b.get());
        bufSize+= b->l_data;
        totalSize+= b->l_data;

        if (res < 0) {
          continueParsing = false;
          break;
        }
        if (IncrementCounts(b.get(), contigLength, nA, nC, nG, nT, nDel)) {
          curEndPos=bam_endpos(b.get());
        }
        ++nReads;
      }

      //
      // Compute coverage for bins
      const int lastBin=(curEndPos-30000)/100;
      assert((lastBin+1)*100 <= nA.size());
      assert((lastBin+1)*100 <= nC.size());
      assert((lastBin+1)*100 <= nG.size());
      assert((lastBin+1)*100 <= nT.size());
      assert((lastBin+1)*100 <= nDel.size());
      for (int binIndex=curCovBin; binIndex < lastBin; binIndex++) {
        int binTot=0;
        for (int nuc=binIndex*100; nuc < (binIndex+1)*100; nuc++) {
          binTot += nA[nuc] + nC[nuc]+ nG[nuc] + nT[nuc] + nDel[nuc];
        }
        covBins.push_back(binTot/100);
      }
      curCovBin=lastBin;

      //
      // Get summary statistics
      //
      double totCov=0;
      double totCovSq=0;
      //
      // First pass gets close to CN=2
      //

      assert(lastBin <= covBins.size());
      for (int binIndex=0; binIndex<lastBin; binIndex++) {
        totCov+=covBins[binIndex];
      }

      mean=totCov/lastBin;
      totCov=0;
      nSamples=0;
      for (int binIndex=0; binIndex<lastBin; binIndex++) {
        if (covBins[binIndex] > 0.25*mean and covBins[binIndex] < 1.75 * mean) {
          totCov+=covBins[binIndex];
          totCovSq+=covBins[binIndex]*covBins[binIndex];
          nSamples++;
        }
      }
      if (nSamples > 0) {
        mean=totCov/nSamples;
        var=totCovSq/nSamples-mean*mean;
        cerr << "Estimating coverage " << nReads
             << " ending at " << curEndPos << '\t'
             << mean << '\t' << var << '\n';
      }
      if (nSamples > 80000) {
        return 1;
      }
    }
  }
  return 1;
}

void InitParams(vector<vector<double>> &covCovTransP,
                vector<vector<double>> &covSnvTransP,
                vector<vector<double>> &snvSnvTransP,
                int nCovStates, int nSNVStates,
                double diag, double offDiag,
                vector<vector<double>> &emisP,
                int model, int maxCov, double mean, double var,
                vector<vector<vector<double>>> &binoP) {
  covCovTransP.resize(nCovStates);
  for (int i=0;i<nCovStates;i++) {
    covCovTransP[i].resize(nCovStates);
    for (int j=0;j<nCovStates;j++) {
      if(i==j) {
        covCovTransP[i][j]  = diag; //log(1 - std::exp(beta) );
      }
      else {
        covCovTransP[i][j] = offDiag;
      }
    }
  }
  const double snvOffDiag=(1-diag)/2;
  covSnvTransP.resize(nCovStates);
  for (int i=0; i < nCovStates; i++) {
    covSnvTransP[i].resize(nSNVStates);
    for (int j=0; j < nSNVStates; j++) {
      if (i == j + 1) {
        covSnvTransP[i][j] = diag;
      }
      else {
        covSnvTransP[i][j] = snvOffDiag;
      }
    }
  }

  snvSnvTransP.resize(nSNVStates);
  for (int i=0; i < nSNVStates; i++) {
    snvSnvTransP[i].resize(nSNVStates);
    for (int j=0; j < nSNVStates; j++) {
      if (i == j + 1) {
        snvSnvTransP[i][j] = diag;
      }
      else {
        snvSnvTransP[i][j] = snvOffDiag;
      }
    }
  }

  emisP.resize(nCovStates);
  for (int i=0;i<nCovStates;i++) {
    emisP[i].resize(maxCov+1);
    double stateSum=0;
    for (int j=0;j<=maxCov;j++) {
      if (model == POIS) {
        emisP[i][j]=LgPrpoiss( (int) i , j , (int) mean/2 );
        stateSum+=exp(LgPrpoiss( (int) i , j , (int) mean/2 ));
      }
      else  {
        emisP[i][j]=LgNegBinom((int)i, (int) j, mean/2, var/2);
        stateSum+=exp(LgNegBinom((int)i, (int) j, mean/2, var/2));
      }
    }
  }

  binoP.resize(3);
  for (size_t i=0; i < binoP.size(); i++) {
    //
    // Create a matrix for each state. Store up to max-cov
    //
    binoP[i].resize(maxCov+1);

    // Initialize matrix.
    //
    for (int j=0; j <= maxCov; j++) {
      binoP[i][j].resize(j+1);
    }
  }

  //
  // Initialize copy number 1
  int i=0;
  const double bino_one=0.9;
  binoP[i][0][0] = bino_one;
  for (int j2=1; j2 < maxCov; j2++) {
    binoP[i][0][0] = log(1/((1-bino_one)/j2));
    for (int k=0; k <= j2; k++) {
      binoP[i][j2][k] = LgBinom(0.1, k, j2);
    }
  }

  // CN=2, diploid
  i=1;
  //
  binoP[i][0][0] = log(bino_one);
  for (int j2=1; j2 < maxCov; j2++) {
    for (int k=0; k <= j2; k++) {
      binoP[i][j2][k] = LgBinom(0.5, k, j2);
    }
  }
  // CN=3, one extra copy
  i=2;
  binoP[i][0][0] = log(bino_one);
  for (int j2=1; j2 < maxCov; j2++) {
    for (int k=0; k <= j2 ; k++) {
      binoP[i][j2][k] = LgBinom(0.66, k, j2);
    }
  }
}

int hmcnc(Parameters& params) {

  double scale=2;
  int maxState=10;
  int averageReadLength=0;

  const string faiFileName{params.referenceName + ".fai"};
  vector<string> contigNames, allContigNames;
  vector<int>    contigLengths, allContigLengths;
  //
  // Determine what chroms to operate on.
  //
  ReadFai(faiFileName, allContigNames, allContigLengths);
  if (params.hmmChrom == "") {
    contigNames = allContigNames;
    contigLengths = allContigLengths;
  }
  else {
    contigNames.push_back(params.hmmChrom);
    for (size_t i=0; i < allContigNames.size(); i++) {
      if (allContigNames[i] == params.hmmChrom) {
        contigLengths.push_back(allContigLengths[i]);
        break;
      }
    }
    if (contigLengths.size() == 0) {
      cerr << "ERROR. Could not find contig for hmm " << params.hmmChrom << '\n';
      exit(1);
    }
  }

  vector<vector<int>> covBins;
  vector<vector<int>> clipBins;
  double mean;
  double var;
  int nStates;
  int maxCov;
  vector<double> startP;
  vector<vector<double>> covCovTransP, covSnvTransP, snvSnvTransP;
  vector<vector<double>> updateTransP;
  vector<vector<SNV>> snvs;
  vector<vector<int>> copyNumber;
  vector< vector<double>> fCov, bCov, fSNV, bSNV;
  vector<vector<double>> emisP;
  vector<vector<double>> updateEmisP;
  vector<vector<vector< double>>> binoP;
  vector<vector<double>> expCovCovTransP, expCovSnvTransP, expSnvSnvTransP, expEmisP;
  vector<int> nReads;
  vector<long> totalBases;

  if (params.covBedInFileName != "") {
    ReadCoverage(params.covBedInFileName, contigNames, covBins);
  }

  if (params.clipInFileName != "") {
    ReadCoverage(params.clipInFileName, contigNames, clipBins);
  }

  if (params.snvInFileName != "") {
    ReadSNVs(params.snvInFileName, contigNames, snvs);
  }

  if (params.paramInFile != "") {
     ReadParameterFile(params.paramInFile, nStates,
		       mean, var, maxState, maxCov,
		       startP, covCovTransP, emisP);
  }

  //
  // Get the header of the bam file
  //
  std::unique_ptr<htsFile, HtslibFileDeleter> htsfp;
  std::shared_ptr<bam_hdr_t> samHeader;
  std::shared_ptr<hts_idx_t> bamidx;


  if (params.bamFileName != "") {
    htsfp.reset(hts_open(params.bamFileName.c_str(),"r"));
    samHeader.reset(sam_hdr_read(htsfp.get()), BamHeaderDeleter{});;
    bamidx.reset(sam_index_load(htsfp.get(), params.bamFileName.c_str()), HtslibIndexDeleter{});
    if (bamidx == nullptr) {
      cerr << "ERROR reading index" << '\n';
      exit(0);
    }
    const htsFormat *fmt = hts_get_format(htsfp.get());
    if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
      cout << "Cannot determine format of input reads." << '\n';
      exit(1);
    }
  }

  //
  // Read index for random io
  //

  params.nproc = min(params.nproc, static_cast<int>(contigNames.size()));

  std::vector<pthread_t> threads(params.nproc);
  std::vector<pthread_attr_t> threadAttr(params.nproc);
  std::vector<ThreadInfo> threadInfo(params.nproc);

  pthread_mutex_t semaphore;
  pthread_mutex_init(&semaphore, NULL);

  int curSeq=0;
  vector<vector<Interval>> copyIntervals;
  copyIntervals.resize(contigNames.size());
  nReads.resize(contigNames.size());
  totalBases.resize(contigNames.size());
  double pModel=0;

  for (int procIndex = 0; procIndex < params.nproc ; procIndex++) {
    if (params.bamFileName != "") {
      threadInfo[procIndex].htsfp.reset(hts_open(params.bamFileName.c_str(),"r"));
    }
    threadInfo[procIndex].bamidx = bamidx;
    threadInfo[procIndex].samHeader=samHeader;
    threadInfo[procIndex].fai.reset(fai_load_format(params.referenceName.c_str(), FAI_FASTA));
    if (params.nproc > 1) {
      threadInfo[procIndex].exit = true;
    }
    else {
      threadInfo[procIndex].exit = false;
    }
    threadInfo[procIndex].lastSeq = &curSeq;
    threadInfo[procIndex].semaphore = &semaphore;
    threadInfo[procIndex].contigNames = &contigNames;
    threadInfo[procIndex].contigLengths = &contigLengths;
    threadInfo[procIndex].covBins = &covBins;
    threadInfo[procIndex].clipBins = &clipBins;
    threadInfo[procIndex].copyNumber = &copyNumber;
    threadInfo[procIndex].snvs = &snvs;
    threadInfo[procIndex].lepsi = lepsi;
    threadInfo[procIndex].scale = scale;
    threadInfo[procIndex].maxState = maxState;
    threadInfo[procIndex].maxCov = maxCov;
    threadInfo[procIndex].mean = mean;
    threadInfo[procIndex].var = var;
    threadInfo[procIndex].transP = &covCovTransP;
    threadInfo[procIndex].expTransP = &expCovCovTransP;
    threadInfo[procIndex].expEmisP = &expEmisP;
    threadInfo[procIndex].emisP = &emisP;
    threadInfo[procIndex].startP = &startP;
    threadInfo[procIndex].copyIntervals = &copyIntervals;
    threadInfo[procIndex].pModel=&pModel;
    threadInfo[procIndex].totalReads = &nReads;
    threadInfo[procIndex].totalBases = &totalBases;
  }

  //
  // If already trained, read those parameters.
  //
  copyNumber.resize(contigLengths.size());

  if (params.paramInFile == "") {
    //max cov value observed or upper cov bound -> max nState---------------
    nStates= std::min( maxState , MAX_CN  ) + 1; //+1 zeroth state
    MAX_CN=nStates+1;
    startP.resize(nStates);
    for(int i=0;i<(nStates);i++) {
      startP[i]=log(1./(nStates));
    }
  }

  //
  // Allocate coverage bins if not reading
  //
  if (params.snvInFileName == "") {
    snvs.resize(contigLengths.size());
  }
  if (params.clipInFileName == "") {
    clipBins.resize(contigLengths.size());
    for (size_t c=0; c < contigLengths.size(); c++ ) {
      clipBins[c].resize(contigLengths[c]/BIN_LENGTH);
    }
  }
  if (params.covBedInFileName == "") {
    covBins.resize(contigLengths.size());

    for (size_t c=0; c < contigLengths.size(); c++ ) {
      covBins[c].resize(contigLengths[c]/BIN_LENGTH);
      clipBins[c].resize(contigLengths[c]/BIN_LENGTH);
      copyNumber[c].resize(contigLengths[c]/BIN_LENGTH);
    }

    //
    // Compute coverage from bam.
    //
    const int parseChromNProc=min(4,params.nproc);
    if (params.nproc > 1) {
      for (int procIndex = 0; procIndex < parseChromNProc; procIndex++) {
        pthread_attr_init(&threadAttr[procIndex]);
        pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*)) ParseChrom, &threadInfo[procIndex]);
      }
      for (int procIndex = 0; procIndex < parseChromNProc ; procIndex++) {
        pthread_join(threads[procIndex], NULL);
      }
      for (int procIndex = 0; procIndex < parseChromNProc ; procIndex++) {
        pthread_attr_destroy(&threadAttr[procIndex]);
      }
    }
    else {
      ParseChrom(&threadInfo[0]);
    }
    long totalBaseSum=0;
    int totalReadsSum=0;
    totalBaseSum=accumulate(totalBases.begin(), totalBases.end(), 0);
    totalReadsSum=accumulate(nReads.begin(), nReads.end(), 0);
    averageReadLength=totalBaseSum/totalReadsSum;
    cerr << "Length cutoff of average read length " << averageReadLength << '\n';
    if (params.covBedOutFileName != "" ) {
      WriteCovBed(params.covBedOutFileName, contigNames, covBins);
    }
    if (params.clipOutFileName != "") {
      WriteCovBed(params.clipOutFileName, contigNames, clipBins);
    }
    if (params.snvOutFileName != "") {
      WriteSNVs(params.snvOutFileName, contigNames, snvs);
    }
  }

  EstimateCoverage(params.bamFileName, covBins, allContigNames, allContigLengths, params.useChrom, mean, var);

  //
  // Cap coverage where hmm does not bother calculating.
  //
  maxCov=(int)mean/2*(maxState+1);

  for (size_t c=0; c < covBins.size(); c++) {
    for (size_t b=0; b < covBins[c].size(); b++) {
      covBins[c][b] = min(covBins[c][b], maxCov);
    }
  }

  //
  // Filter SNVs that are too close
  //
  assert(snvs.size() <= contigNames.size());
  for (size_t c=0 ;c < contigNames.size(); c++) {
    vector<bool> tooClose(snvs[c].size(), false);
    for (size_t i=1; i < snvs[c].size(); i++ ) {
      if (snvs[c][i].pos - snvs[c][i-1].pos < 50) {
        tooClose[i] = true;
        tooClose[i-1] = true;
      }
    }
    int p=0;
    for (size_t i=0; i < snvs[c].size(); i++) {
      if (tooClose[i] == false) {
        snvs[c][p] = snvs[c][i];
        p++;
      }
    }
    snvs[c].resize(p);
  }

  cerr << "Computing copy-number up to " << MAX_CN << '\n';
  //----------------------------------------------------------------------

  if (mean == 0) {
    std::cerr << "mean is zero, Exiting" << '\n';
    return EXIT_FAILURE;
  }

  // trans prob, scale 3->2 by overlap of pdf.

  const poisson distribution1(3*mean/2);
  const double result3=pdf(distribution1, 3*mean/2);

  const poisson distribution2(2*mean/2);
  const double result2=pdf(distribution2, 3*mean/2);

  const double epsi23 = result2/result3;

  //*300

  //no epsi
  const double small=-30;
  //  double beta =  small + ( scale * log(epsi23))  ;
  const double beta=small;
  //mean no. of bins for cn=3 call
  const int nSNVStates=3;
  const double unif=log(1.0/nStates);
  if (params.paramInFile == "") {
    InitParams(covCovTransP, covSnvTransP, snvSnvTransP,
	       nStates, nSNVStates, log(1-exp(small)), log(exp(small)/(nStates-1)),
	       emisP, params.model, maxCov, mean, var, binoP);
    PrintModel("TRANS", covCovTransP, std::cerr);
    PrintModel("EMIS", emisP, std::cerr);
    //    printEmissions(emisP);

    //
    // Baum-Welch training.
    //

    double prevPX=0;
    assert(!emisP.empty());
    vector<vector<double> > prevTransP, prevEmisP;
    for (int i=0; i < 4; i++) {
      prevTransP=covCovTransP;
      prevEmisP=emisP;

      vector<double> stateWeightedTotCov(nStates, 0),
      stateWeightedTotVar(nStates, 0),
      stateTotProb(nStates, 0);
      vector<long> stateTotCov(nStates, 0), stateNCov(nStates, 0);

      expCovCovTransP.resize(nStates);
      expCovSnvTransP.resize(nStates);
      expSnvSnvTransP.resize(nStates);
      expEmisP.resize(nStates);
      for (int r=0; r < nStates; r++ ) {
        expCovCovTransP[r].resize(nStates,0);
        expCovSnvTransP[r].resize(nStates,0);
        expSnvSnvTransP[r].resize(nStates,0);
        fill(expCovCovTransP[r].begin(), expCovCovTransP[r].end(), 0);
        expEmisP[r].resize(emisP[0].size());
      }
      double px=0;
      int totalObs=0;
      curSeq=0;
      for (int procIndex = 0; procIndex < params.nproc; procIndex++) {
        pthread_attr_init(&threadAttr[procIndex]);
        pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*)) ThreadedBWE, &threadInfo[procIndex]);
      }
      for (int procIndex = 0; procIndex < params.nproc ; procIndex++) {
        pthread_join(threads[procIndex], NULL);
      }
      for (int procIndex = 0; procIndex < params.nproc ; procIndex++) {
        pthread_attr_destroy(&threadAttr[procIndex]);
      }
      px = pModel;

      if (prevPX != 0 and px - prevPX < 100 and i > 1) {
        cerr << "Ending iteration after " << i << " steps" << '\n';
        covCovTransP = prevTransP;
        emisP  = prevEmisP;
        break;
      }
      prevPX=px;

      vector<vector<double> > priorCovCov;
      priorCovCov.resize(covCovTransP.size());
      int nSites=covBins[0].size();
      BaumWelchM(startP, covCovTransP, emisP, binoP,
        params.model,
        stateTotCov, stateNCov,
        expCovCovTransP, expEmisP,
        priorCovCov,
        updateTransP, updateEmisP);

      PrintModel("TRANS", updateTransP, std::cerr);
      PrintModel("EMIS", updateEmisP, std::cerr);
      covCovTransP=updateTransP;
    }

    //
    // Eventually this needs to update for some multi-chrom code.
    //

    if (params.paramOutFile != "") {
      WriteParameterFile(params.paramOutFile, nStates, mean, var, maxState, maxCov, startP, covCovTransP, emisP);
    }
  }

  //
  // Now filter cn=1 and cn=3 intervals, or short calls
  //
  assert(copyIntervals.size() <= contigNames.size());
  assert(snvs.size() <= contigNames.size());
  for (size_t c=0; c < contigNames.size(); c++) {
    int snvStart=0;
    for (size_t i=0; i < copyIntervals[c].size(); i++) {
      const int curCN =copyIntervals[c][i].copyNumber;
      if ( curCN == 1 or curCN == 3) {
        while (snvStart < snvs[c].size() and snvs[c][snvStart].pos < copyIntervals[c][i].start) {
          snvStart++;
        }
        int snvEnd=snvStart;
        while (snvEnd < snvs[c].size() and snvs[c][snvEnd].pos < copyIntervals[c][i].end) {
          snvEnd++;
        }

        double pCN=0, pCN2=0;
        assert(snvEnd <= snvs[c].size());
        for (int cni=snvStart; cni < snvEnd; cni++ ) {
          int ref, alt;
          ref=snvs[c][cni].ref;
          alt=snvs[c][cni].alt;
          if (ref+alt >= maxCov) {
            alt=(int)(((float)maxCov)/(ref+alt) * alt);
            ref=(int)(((float)maxCov)/(ref+alt) * ref);
          }
          int totCov=ref+alt;
          if (totCov >= maxCov) {
            totCov = maxCov-1;
          }
          if (alt >= maxCov) {
            alt = maxCov-1;
          }
          pCN += binoP[curCN-1][totCov][alt];
          pCN2 += binoP[1][totCov][alt];
        }
        if (pCN < pCN2) {
          copyIntervals[c][i].filter = "FAIL";
        }
        if (averageReadLength > 0 and copyIntervals[c][i].end-copyIntervals[c][i].start *2 < averageReadLength) {
          copyIntervals[c][i].filter = "FAIL";
        }
	      snvStart=snvEnd;
      }
    }
  }

  ostream *outPtr;
  ofstream outFile;
  if (params.outFileName != "") {
    outFile.open(params.outFileName.c_str());
    outPtr = &outFile;
  }
  else {
    outPtr = &cout;
  }
  WriteVCF(*outPtr, params.referenceName, params.sampleName, contigNames, contigLengths, copyIntervals);

  return 0;
}

int hmcnc(int argc, const char* argv[]) {

  //
  // Initialize parameters from command line
  //
  Parameters params;
  CLI11_PARSE(params.CLI, argc, argv);
  return hmcnc(params);
}

int hmcnc_test() { return 42; }
