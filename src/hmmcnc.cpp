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

using boost::math::binomial;

using boost::math::poisson;
using boost::math::pdf;
using boost::math::negative_binomial_distribution;
using std::vector;
using std::cout;
using std::string;
using std::log;
using namespace std;

const int MIN_ECLIP=0;
const int MIN_FCLIP=0;
const int MIN_SV_LENGTH=1000;
const int MIN_DEL_LENGTH=5000;
const int MIN_MAPQ=10;
const int BIN_LENGTH=100;
const int MIN_CLIP_LENGTH=500;
const int MIN_CONTIG_SIZE=1000000;
const float MISMAP_RATE=0.01;
const double minNeg=-1*numeric_limits<double>::epsilon();
int MAX_CN=6;
double lepsi=-800;
int averageReadLength=0;
bool writeFail=true;

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


// -----------
// SNV
// -----------

SNV::SNV() {
  assert(0);
}

SNV::SNV(int p, int r, int a, int rc, int ac)
  : pos(p)
  , refNuc(r)
  , altNuc(a)
  , ref(rc)
  , alt(ac)
{ }

SNV::SNV(int p)
  : pos(p)
{ }

bool SNV::operator<(const SNV &rhs) const {
  return pos < rhs.pos;
}

// -----------

void StorePerChromAverageCoverage(vector<vector<int>  > &covBins, vector<double> &averageCoverage) {
  for (int i=0; i < covBins.size(); i++) {
    long totCov=0;
    int nonZero=0;
    for (auto &c: covBins[i]) {
      if (c > 0) {
      	totCov+=c;
      	nonZero++;
      }
    }
    if (nonZero > 0) {
      averageCoverage[i] = totCov / ((double)nonZero);
    }
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

// -----------
// Interval
// -----------

Interval::Interval() {
  assert(0);
}

Interval::Interval(int s, int e, int cn, float avg, double p)
  : start(s)
  , end(e)
  , copyNumber(cn)
  , averageCoverage(avg)
  , pVal(p)
  , filter("PASS")
  , altInfo("")
  , altSample("")
  , distanceToFrontClip(-1)
  , distanceToEndClip(-1)
  , nFrontClip(0)
  , nEndClip(0)

{ }

// -----------

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
  vector<vector<double>> *cl,*n;

  vector<vector<SNV>> *snvs;
  vector<vector<int>> *copyNumber;
  vector<vector<double>> *transP, *emisP, *expTransP, *expClipTransP, *expEmisP, *clTransP;
  vector<double> *startP;
  vector<vector<Interval>> *copyIntervals;
  vector<double> *chromCopyNumber;
  vector<vector<Interval>> *mergedNaiveIntervals;
  vector<vector<Interval>> *UnmergedNaiveIntervals;
  vector<vector<Interval>> *delT;

  int maxCov, maxState;
  bool exit;
  double mean;
  double var;
  double lepsi;
  double scale;
  double *pModel;
  string hmmChrom;
  vector<int> *totalReads;
  vector<long> *totalBases;
  vector<double> *averageCoverage;
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


double LgNegBinom(int cn, int cov, float Hmean, float Hvar) {
  double result=0;
  float r, p;


  if ((Hmean/Hvar)>=0.90){
    return (LgPrpoiss(cn,cov, (int) Hmean ));
  }
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
  if (cn==0) {//del_states
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
    cerr << "ERROR computing obs index." << '\n';
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
                        const std::vector<std::vector<double>> &emisP,
                        const std::vector<int> &obs,
                        std::vector<std::vector<double>> &f,
                        std::vector<std::vector<double>> &b) {
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
    f[j][0] = startP[j];
  }

  const double lgthird=log(1/3.);

  if (nCovStates == 1) {
    fill(f[0].begin(), f[0].end(), -1);
    fill(b[0].begin(), b[0].end(), -1);
    return 0;
  }
  int prevCovIdx=0, curCovIdx=0;
  int prevSNVIdx=0, curSNVIdx=0;
  double clipSum=0;
  for (int k=0; k < totObs; k++) {
    for (int i=0; i < nCovStates; i++) {
      double colSum=0;
      for (int j=0; j < nCovStates; j++) {
        assert(j== 0 or colSum != 0);
        assert(j < f.size());
        assert(k < f[j].size());
        assert(j < covCovTransP.size());
        assert(i < covCovTransP[j].size());
        clipSum = covCovTransP[i][j];
        colSum = PairSumOfLogP(colSum, f[j][k] + clipSum);
      }
      assert(obs[k] < emisP[i].size());
      f[i][k+1] = colSum + emisP[i][obs[k]];
    }
  }

  double finalCol=0;
  for (int j=0; j < nCovStates; j++) {
    finalCol = PairSumOfLogP(finalCol, f[j][totObs]+ log(1./nCovStates));
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
        clipSum = covCovTransP[i][j];
        colSum = PairSumOfLogP(colSum, b[j][k+1] + clipSum + emisP[j][obs[k]]);
      }
      b[i][k] = colSum;
    }
  }

  double bfinalCol=0;
  for (int j=0; j < nCovStates; j++) {
    bfinalCol = PairSumOfLogP(bfinalCol, b[j][1] + emisP[j][obs[0]] + log(1./nCovStates));
  }
  assert(finalCol==bfinalCol);


  return finalCol;
}




double ForwardBackwards(const vector<double> &startP,
                        const vector<vector<double>> &covCovTransP,
                        const std::vector<std::vector<double>> &clipCovCovTransP,
                        const std::vector<std::vector<double>> &emisP,
                        const std::vector<int> &obs,
                        std::vector<std::vector<double>> &f,
                        std::vector<std::vector<double>> &b,
                        std::vector<double> &Pn, std::vector<double> &Pcl) {
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
  //  and Pcl[0...nObs-1] , Pn[0...nObs-1]
  //  nObs=totObs
  
  f.resize(nCovStates);
  b.resize(nCovStates);

  for (int i = 0; i < nCovStates; i++) {
    f[i].resize(totObs+1);
    fill(f[i].begin(), f[i].end(), 0);
    b[i].resize(totObs+1);
    fill(b[i].begin(), b[i].end(), 0);
  }

  for (int j=0; j < nCovStates; j++) {
    f[j][0] = startP[j];
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
  double clipSum=0;
  for (int k=0; k < totObs; k++) {
    for (int i=0; i < nCovStates; i++) {
      double colSum=0;
      for (int j=0; j < nCovStates; j++) {
        assert(j== 0 or colSum != 0);
        assert(j < f.size());
        assert(k < f[j].size());
        assert(j < covCovTransP.size());
        assert(i < covCovTransP[j].size());
        assert(j < clipCovCovTransP.size());
        assert(i < clipCovCovTransP[j].size());
        clipSum = PairSumOfLogP(covCovTransP[i][j] + Pn[k] , clipCovCovTransP[i][j] + Pcl[k] );
        colSum = PairSumOfLogP(colSum, f[j][k] + clipSum);
      }
      assert(obs[k] < emisP[i].size()); //capped coverage
      f[i][k+1] = colSum + emisP[i][obs[k]];
    }
  }

  double finalCol=0;
  for (int j=0; j < nCovStates; j++) {
    finalCol = PairSumOfLogP(finalCol, f[j][totObs]+ log(1./nCovStates));
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
        clipSum = PairSumOfLogP(covCovTransP[i][j] + Pn[k] , clipCovCovTransP[i][j] + Pcl[k] );
        colSum = PairSumOfLogP(colSum, b[j][k+1] + clipSum + emisP[j][obs[k]]);
      }
      b[i][k] = colSum;
    }
  }

  double bfinalCol=0;
  for (int j=0; j < nCovStates; j++) {
    bfinalCol = PairSumOfLogP(bfinalCol, b[j][1] + emisP[j][obs[0]] + log(1./nCovStates));
  }
  assert(finalCol==bfinalCol);

  return finalCol;
}

void ApplyPriorToClipTransP(vector<vector<int> > &f,
			int nStates,
			vector<vector<double> > &prior,
			vector<vector<double>>  &expCovCovTransP) {
  //
  // For now the prior is hard-wired for human genomes.
  //
  prior.resize(nStates);
  for (int i=0; i < nStates; i++) {
    prior[i].resize(nStates, 0);
  }
  // Assume the following:
  //
  if (nStates < 4) {
    cerr << "ERROR number of states should be at least 4" << endl;
    exit(1);
  }
  //
  // Prior is to keep in the same state.
  //
  for (int contig=0; contig< f.size(); contig++) {
    int nBins=f[contig].size();
    for (int i=0; i < nStates; i++) {
      for (int j=0; j < nStates; j++ ) {
    		if (i == j) {
    		  expCovCovTransP[i][j] += 4000;
    		}
        else if(j==2){
          expCovCovTransP[i][j] += 1000;          
        }
    		else if(j==1 or j==3){
    		  expCovCovTransP[i][j] += 10;
    		}
        else
          expCovCovTransP[i][j] += 100;          
      }
    }
  }
  // Expext roughly this many transitions in a mammalian genome.
  //  expCovCovTransP[2][1] += 1000;
  //  expCovCovTransP[2][3] += 100;
}

void ApplyPriorToTransP(vector<vector<int> > &f,
      int nStates,
      vector<vector<double> > &prior,
      vector<vector<double>>  &expCovCovTransP) {
  //
  // For now the prior is hard-wired for human genomes.
  //
  prior.resize(nStates);
  for (int i=0; i < nStates; i++) {
    prior[i].resize(nStates, 0);
  }
  // Assume the following:
  //
  if (nStates < 4) {
    cerr << "ERROR number of states should be at least 4" << endl;
    exit(1);
  }
  //
  // Prior is to keep in the same state.
  //
  for (int contig=0; contig< f.size(); contig++) {
    int nBins=f[contig].size();
    for (int i=0; i < nStates; i++) {
      for (int j=0; j < nStates; j++ ) {
        if (i == j) {
          expCovCovTransP[i][j] += max(500000,nBins*500);
        }
        else if(j==2){
          expCovCovTransP[i][j] += max(100000,nBins*100);          
        }
        else if(j==1 or j==3){
          expCovCovTransP[i][j] += 100;
        }
        else
          expCovCovTransP[i][j] += 1000;          
      }
    }
  }
  // Expext roughly this many transitions in a mammalian genome.
  //  expCovCovTransP[2][1] += 1000;
  //  expCovCovTransP[2][3] += 100;
}

double BaumWelchEOnChrom(const vector<double> &startP,
			 vector<vector<double>> &covCovTransP, vector<vector<double>> &clipCovCovTransP,
			 vector<vector<double>> &emisP,
			 vector<int> &obs,
			 vector<vector<double>> &f,
			 vector<vector<double>> &b,
			 vector<vector<double>> &expCovCovNoClipTransP, vector<vector<double>> &expCovCovClipTransP,
			 vector<vector<double>> &expEmisP,
       vector<double> &Pn, vector<double> &Pcl) {

  const int nStates = static_cast<int>(startP.size());
  const int nObs = obs.size();
  const int nclObs = Pcl.size();
  const int nNObs = Pn.size();

  assert(nNObs==nclObs);
  assert(nObs==nclObs);

  const double pxNoclip = ForwardBackwards(startP, covCovTransP, emisP, obs, f, b);

  const double pxClip = ForwardBackwards(startP, covCovTransP, clipCovCovTransP, emisP, obs, f, b, Pn , Pcl);

  std::cerr<<"\npxNoclip: "<<pxNoclip<<"\npxClip: "<<pxClip<<endl;

  double pNoClipEdge = 0;
  double pClipEdge = 0;


  for (int k=1; k< nObs-1; k++) {
    //
    // Calculate total probability of all transitions at this step.
    //
    //
    double logClipSum = 0;
    double logNoClipSum = 0;
    double transPsum = 0;
    
    
    for (size_t i=0; i < covCovTransP.size(); i++) {
      for (size_t j=0; j < covCovTransP[i].size(); j++) {
        pNoClipEdge = f[i][k] + covCovTransP[i][j] + Pn[k] + emisP[j][obs[k]] + b[j][k+1];
        pClipEdge   = f[i][k] + clipCovCovTransP[i][j] + Pcl[k] + emisP[j][obs[k]] + b[j][k+1];
    		expCovCovNoClipTransP[i][j] += exp(pNoClipEdge - pxClip);
    		expCovCovClipTransP[i][j]   += exp(pClipEdge - pxClip);
      }
    }

  }
  //
  // For now do not tune parameters for emission.
  //
  return pxClip;
}

void AssignNearestClip(vector<vector<int > > &clipBins,
		       double averageCoverage,
		       int maxSearch,
		       vector<string> &contigNames,
		       vector<vector<Interval> > &intervals) {
  int minClip=averageCoverage/4;
  for (int contig=0; contig < intervals.size(); contig++) {
    int lastClip=0;
    int curInterval=0;
    for (int i=0; i < intervals[contig].size(); i++) {
      int maxClip=0;
      int maxClipPos=-1;
      int intvStart = intervals[contig][i].start/100;
      int intvEnd   = intervals[contig][i].end/100;
      // Do a bit of logic to determine when to stop the search
      int searchEnd = min((int)clipBins[contig].size(),
			  min(max(intvStart, intvEnd -maxSearch),
			      intvStart+maxSearch));

      for (int clipIndex=max(0, intvStart - maxSearch); clipIndex < searchEnd; clipIndex++) {
        if (clipBins[contig][clipIndex] > maxClip and clipBins[contig][clipIndex] > minClip) {
          maxClip = clipBins[contig][clipIndex];
          maxClipPos=clipIndex;
        }
      }
      if (maxClipPos != -1) {
        intervals[contig][i].distanceToFrontClip=abs(maxClipPos - intvStart);
        intervals[contig][i].nFrontClip=maxClip;

/*        cerr << "Found start clip for " << contigNames[contig] << "\t"
             << i << "\t"
             << intervals[contig][i].distanceToFrontClip << "\t"
             << intervals[contig][i].nFrontClip << "\t"
             << contigNames[contig] << ":"
             << intervals[contig][i].start << "-"
             << intervals[contig][i].end  << "\t"
             << intervals[contig][i].copyNumber << endl;
 */ 
      }

      // Look for clipping at the end of the interval
      maxClip=0;
      maxClipPos=-1;
      for (int clipIndex = max(searchEnd, intvEnd - maxSearch);
	         clipIndex < min((int) clipBins[contig].size(), intvEnd + maxSearch); clipIndex++) {
        if (clipBins[contig][clipIndex] > maxClip and clipBins[contig][clipIndex] > minClip) {
          maxClip = clipBins[contig][clipIndex];
          maxClipPos=clipIndex;
        }
      }
      if (maxClipPos != -1) {
        intervals[contig][i].distanceToEndClip=abs(maxClipPos - intvEnd);
        intervals[contig][i].nEndClip=maxClip;

/*        cerr << "Found End clip for " << contigNames[contig] << "\t"
             << i << "\t"
             << intervals[contig][i].distanceToEndClip << "\t"
             << intervals[contig][i].nEndClip << "\t"
             << contigNames[contig] << ":"
             << intervals[contig][i].start << "-"
             << intervals[contig][i].end <<  "\t"
             << intervals[contig][i].copyNumber << endl;
  
*/
      }
    }
  }
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


bool compareInterval(Interval i1, Interval i2) {return (i1.start < i2.start);}

bool compareIntervalLength(Interval i1, Interval i2) {return (i1.end - i1.start) < (i2.end - i2.start);}

void mergeIntervals(vector<Interval> & intervals, vector<Interval> &mergedIntervals, string contig) {
  
  const int n = intervals.size();
  
  std::sort(intervals.begin(), intervals.end() , compareInterval);

  for (int i = 0; i < n - 1; i++) {

    if ( (intervals[i].start >= intervals[i + 1].start && intervals[i].start <= intervals[i + 1].end) || 
      (intervals[i].end>= intervals[i + 1].start && intervals[i].end<= intervals[i + 1].end)) {

      intervals[i + 1].start = std::min(intervals[i].start, intervals[i + 1].start);
      intervals[i + 1].end= std::max(intervals[i].end, intervals[i + 1].end);
      // Remove previous interval
      intervals[i].start = -1;
      intervals[i].end = -1;
      //count number of del cigars 
      intervals[i+1].averageCoverage = intervals[i].averageCoverage + 1.0;
    }
  }
  for (int i = 0; i < n; i++) {
    if (!(intervals[i].start == -1 && intervals[i].end== -1)) {
      mergedIntervals.push_back(Interval(intervals[i].start,intervals[i].end , 0, intervals[i].averageCoverage , 0.0  ));
    }
  }
}



void mergeNaiveIntervals(vector<Interval> &intervals, vector<Interval> &mergedIntervals, string contig) {

  //std::sort(intervals.begin(), intervals.end() , compareInterval);
  const int n = intervals.size();
  for (int i = 0; i < n - 1; i++) {

    if ( ((intervals[i].start >= intervals[i + 1].start && intervals[i].start <= intervals[i + 1].end) || (intervals[i].end>= intervals[i + 1].start && intervals[i].end<= intervals[i + 1].end)) && intervals[i].copyNumber == intervals[i+1].copyNumber ) {
      intervals[i + 1].start = std::min(intervals[i].start, intervals[i + 1].start);
      intervals[i + 1].end= std::max(intervals[i].end, intervals[i + 1].end);
      // Remove previous interval
      intervals[i].start = -1;
      intervals[i].end= -1;
    }
  }
  for (int i = 0; i < n; i++) {
    if (!(intervals[i].start == -1 && intervals[i].end== -1)) {
      mergedIntervals.push_back(Interval(intervals[i].start,intervals[i].end , intervals[i].copyNumber, 0.0, 0.0  ));
    }
  }
}



void NaiveCaller(vector<int> &covBins, vector<Interval> & UnmergedNaiveIntervals, double mean ){

  const int bins = covBins.size();
  const double Hmean = mean/2;

  for (int i=0;i<bins;i++){
    int cn=round(covBins[i]/Hmean);
    if (cn==0 && (covBins[i] >= Hmean/2) ){cn=1;}
    UnmergedNaiveIntervals.push_back( Interval(i*BIN_LENGTH, (i+1)*BIN_LENGTH, cn, (float) covBins[i], 0.0));
  }
}

void intersectDelCall( vector<Interval> &mergedIntervals, vector<Interval> & copyIntervals, double mean)
{

  int i = 0, j = 0;
  int n = mergedIntervals.size(), m = copyIntervals.size();
  //std::sort(copyIntervals.begin(), copyIntervals.end() , compareInterval);

  while (i < n && j < m)
  {
    if (copyIntervals[j].copyNumber > 1 )
    {
        j++;
    }
    else
    {

      // Left bound for intersecting segment
      int l = max(mergedIntervals[i].start, copyIntervals[j].start);
      // Right bound for intersecting segment
      int r = min(mergedIntervals[i].end, copyIntervals[j].end);
      int inter = r-l;
      if (inter > 0)
      {
        //float lr = (float)((inter)/(mergedIntervals[i].end - mergedIntervals[i].start));
        //float rr = (float)((inter)/(copyIntervals[i].end - copyIntervals[i].start));
        copyIntervals[j].start = mergedIntervals[i].start;
        copyIntervals[j].end = mergedIntervals[i].end;
        copyIntervals[j].altInfo += ":DF";
        copyIntervals[j].altSample += ":1" ;        
        double diff = ((double) copyIntervals[j].averageCoverage) - ((double) mergedIntervals[i].averageCoverage) ;
        if ( diff  > 2 ) //support for CN=1 can be much lower than (mean*0.25)   )
        {
          copyIntervals[j].copyNumber = 1;
          copyIntervals[j].averageCoverage = diff;
        }
        else{
          copyIntervals[j].copyNumber = 0;
          copyIntervals[j].averageCoverage = diff;
        }
        
        i++;j++;

      }
      else if (mergedIntervals[i].end < copyIntervals[j].end)
      {
          i++;
      }
      else if (mergedIntervals[i].end > copyIntervals[j].end)
      {
        copyIntervals[j].altInfo += ":DF";
        copyIntervals[j].altSample += ":0" ;
        j++;
      }
    }
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
		const vector<vector<double>> &expNoClipTransP, const vector<vector<double>> &expClipTransP,
		vector<vector<double>> &expEmisP,
		vector<vector<double>> &covCovPrior,
		vector<vector<double>> &updateNoClipTransP, vector<vector<double>> &updateClipTransP,
		vector<vector<double>> &updateEmisP) {

  //
  // M step.
  //
  const int nStates=static_cast<int>(startP.size());
  updateNoClipTransP.resize(nStates);
  updateClipTransP.resize(nStates);


  for (int j=0; j < nStates; j++) {
    double noClipColSum=0;
    double clipColSum=0;
    updateNoClipTransP[j].resize(nStates);
    updateClipTransP[j].resize(nStates);
    assert(nStates <= expNoClipTransP[j].size());
    assert(nStates <= expClipTransP[j].size());

    for (int k=0; k< nStates; k++) {
      noClipColSum += expNoClipTransP[j][k];
      clipColSum   += expClipTransP[j][k];
    }
    
    //cerr << j;
    for (int k=0; k < nStates; k++) {
      //      assert(expClipTransP[j][k] - clipColSum < 0); 
      //      assert(expNoClipTransP[j][k] - noClipColSum < 0);
      updateClipTransP[j][k] = log( expClipTransP[j][k] / clipColSum);
      updateNoClipTransP[j][k] = log( expNoClipTransP[j][k] / noClipColSum);
    }
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
  if (readLength < BIN_LENGTH or b->core.qual < MIN_MAPQ or b->core.flag & 0x800) {
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
int IncrementCounts(bam1_t *b, int contigLength,
                    vector<int> &nA,
                    vector<int> &nC,
                    vector<int> &nG,
                    vector<int> &nT,
                    vector<int> &nDel,
                    vector<Interval> &delT) {
  const int readLength = b->core.l_qseq;
  if (readLength < BIN_LENGTH or b->core.qual < MIN_MAPQ or b->core.flag & 0x800) {
    return 0;
  }

 // string ctg = b->core.tid;
//  string rdnm = b->core.l_qname;

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
      int stop=refPos+opLen;
      if (opLen >= MIN_DEL_LENGTH ){
        delT.push_back(Interval( refPos, stop,0,0.0,0.0)); //, ctg, rdnm ))
      }
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
    if (threadInfo->hmmChrom != "" and (*(*threadInfo).contigNames)[curSeq] != threadInfo->hmmChrom) {
      break;
    }
    if (curSeq >= threadInfo->contigNames->size()) {
      break;
    }
    if ((*(*threadInfo).averageCoverage)[curSeq] == 0) {
      break;
    }

    vector<vector<double>> f, b, expCovNoClipTransP, expCovClipTransP, expEmisP;
    expCovNoClipTransP.resize(threadInfo->transP->size());
    expCovClipTransP.resize(threadInfo->transP->size());
    
    for (int i=0; i < expCovNoClipTransP.size(); i++) {
      expCovNoClipTransP[i].resize(expCovNoClipTransP.size(), 0);
      expCovClipTransP[i].resize(expCovNoClipTransP.size(), 0);

    }
    expEmisP.resize(threadInfo->emisP->size());
    for (int i=0; i < expEmisP.size(); i++){
      expEmisP[i].resize((*(threadInfo->emisP))[i].size());
    }

    pChrom = BaumWelchEOnChrom(*threadInfo->startP,
                                *threadInfo->transP, *threadInfo->clTransP,
                                *threadInfo->emisP,
                                (*threadInfo->covBins)[curSeq],
                                f, b,
                                expCovNoClipTransP, expCovClipTransP,
                                expEmisP,
                                (*threadInfo->n)[curSeq], (*threadInfo->cl)[curSeq]);

    //
    // Update expected transitions
    //
    pthread_mutex_lock(threadInfo->semaphore);

    if ((*threadInfo->chromCopyNumber)[curSeq] >= 1.5 and (*threadInfo->chromCopyNumber)[curSeq] < 2.5)
    {
      for (size_t i=0; i < threadInfo->transP->size(); i++) {
	      for (size_t j=0; j < (*threadInfo->transP)[i].size(); j++) {
	        (*threadInfo->expTransP)[i][j] += expCovNoClipTransP[i][j];
          	(*threadInfo->expClipTransP)[i][j] += expCovClipTransP[i][j];          
        }
      }
      for (size_t i=0; i < threadInfo->emisP->size(); i++) {
	      for (size_t j=0; j < (*threadInfo->emisP)[i].size(); j++) {
	        (*threadInfo->expEmisP)[i][j] += expEmisP[i][j];
	      }
      }
      *threadInfo->pModel  += pChrom;
    }
    pthread_mutex_unlock(threadInfo->semaphore);
    StorePosteriorMaxIntervals((*threadInfo->covBins)[curSeq],
			       f, b,
			       (*threadInfo->copyIntervals)[curSeq]);
    
    /*
    debug
    for (size_t j=0; j<f[0].size(); j++){
      cout<<(*threadInfo->contigNames)[curSeq]<<"\t"<<j<<"\t";
      for (size_t i=0; i<f.size(); i++){
        cout<<f[i][j]<<"\t";
      }
      cout<<endl;
    }
    cout<<endl;
    cerr << "Stored " << (*threadInfo->copyIntervals)[curSeq].size()
         << " copy intervals for " << (*threadInfo->contigNames)[curSeq] << endl;
  //debug
  }
}

void ParseChrom(ThreadInfo *threadInfo) {

  while (*(threadInfo->lastSeq) < (*(*threadInfo).contigNames).size()) {
    //
    // Grab current chrom to process
    //
    pthread_mutex_lock(threadInfo->semaphore);

    const int curSeq = *((*threadInfo).lastSeq);
    if (curSeq >= threadInfo->contigNames->size()) {
      pthread_mutex_unlock(threadInfo->semaphore);
      break;
    }

    *(threadInfo->lastSeq) = *(threadInfo->lastSeq) + 1;
    (*(*threadInfo).totalReads)[curSeq] = 0;
    (*(*threadInfo).totalBases)[curSeq] = 0;

    //
    // Deal with race condition by double checking curSeq;
    //

    (*threadInfo).procChroms.push_back(curSeq);
    //    (*threadInfo).covBins->push_back(vector<int>());
    //    (*threadInfo).snvs->push_back(vector<SNV>());
    pthread_mutex_unlock(threadInfo->semaphore);

    const int contigLength=threadInfo->contigLengths->at(curSeq);

    vector<int> nA(contigLength, 0), nC(contigLength, 0), nT(contigLength, 0), nG(contigLength,0), nDel(contigLength, 0);

    stringstream regionStrm;
    string contigName=(*(*threadInfo).contigNames)[curSeq];
    regionStrm << (*(*threadInfo).contigNames)[curSeq];// << ":1-" << contigLength;

    const string region=regionStrm.str();

    std::unique_ptr<hts_itr_t, HtslibIteratorDeleter> regionIter(
      sam_itr_querys(threadInfo->bamidx.get(), threadInfo->samHeader.get(), region.c_str()));

    int tid = sam_hdr_name2tid(threadInfo->samHeader.get(), contigName.c_str());
    if (tid < 0) {
      cerr << "ERROR. Contig " << contigName << " is not a reference sequence in the input bam." << std::endl;
      exit(1);
    }
    uint64_t idxTotalBases, idxTotalReads;
    hts_idx_get_stat(threadInfo->bamidx.get(), tid, &idxTotalBases, &idxTotalReads);

    int chromLen;
    char *chromSeq = fai_fetch(threadInfo->fai.get(), region.c_str(), &chromLen);

    bool continueParsing=true;
    vector<std::unique_ptr<bam1_t, BamRecordDeleter>> reads; //(bam_init1());
    long totalSize=0;
    int chunkNumber=0;
    uint64_t totalBases=0;
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
        endpos=bam_endpos(b.get());
        if ((b->core.flag & BAM_FSUPPLEMENTARY) == 0 && (b->core.flag & BAM_FSECONDARY) == 0) {
	  totalBases+=b->core.l_qseq;
          reads.push_back(std::move(b));
        }
      }
      cerr << "Reading " << (*threadInfo->contigNames)[curSeq] << ", chunk " << chunkNumber << ".\t" << reads.size() << "/" << totalBases << "/" << idxTotalBases << " reads/net/total" << '\n';
      ++chunkNumber;
      pthread_mutex_unlock(threadInfo->semaphore);

      for (auto& b : reads) {
        IncrementCounts(b.get(), contigLength, nA, nC, nG, nT, nDel, (*(threadInfo)->delT)[curSeq]);
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

        if (nCigar > 1 and bam_cigar_op(cigar[nCigar-1]) == BAM_CSOFT_CLIP) {
          backClip=nCigar-1;
          backClipLen=bam_cigar_oplen(cigar[nCigar-1]);
        }
        else if (nCigar > 2 and bam_cigar_op(cigar[nCigar-2]) == BAM_CSOFT_CLIP) {
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
    if ((*threadInfo->totalReads)[curSeq] > 0) {
      (*threadInfo->averageCoverage)[curSeq] = (*threadInfo->totalBases)[curSeq]/((double)(*threadInfo->totalReads)[curSeq]);
    }
    else {
      (*threadInfo->averageCoverage)[curSeq] = 0;
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

void ChromCopyNumber(const vector<vector<int>> &allCovBins,
		     double mean,
		     vector<double> &chromCopy) {
  chromCopy.resize(allCovBins.size(), 0);
  for (auto c=0; c < allCovBins.size(); c++) {
    if (allCovBins[c].size() > 0) {
      double totChromCov= std::accumulate(allCovBins[c].begin(), allCovBins[c].end(), 0);
      totChromCov /= allCovBins[c].size();
      chromCopy[c] = 2*totChromCov/mean;
    }
    else {
      chromCopy[c] = 0;
    }
  }
}


int EstimateCoverage(const string &bamFileName,
                     const vector<vector<int>> &allCovBins,
                     const vector<string> &chroms,
                     const vector<int> &lengths,
                     string &useChrom,
                     double &mean,
                     double &var)
{
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
    int binCov=0;
    int nSamples=0;
    for (size_t binIndex=0; binIndex<lastBin; binIndex++) {
      binCov=allCovBins[useChromIndex][binIndex];
      totCov+=binCov;
      if (binCov > 0) { nSamples++;}
    }
    if (nSamples > 0) {
      mean=totCov/((float)nSamples);
    }
    else {
      cerr << "Could not estimate coverage using precomputed coverage on chrom " << useChrom << " because there were 0 covered bins\n";
      exit(1);
    }
    //
    // Recompute summary stats using limited data
    nSamples=0;
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
      cerr << "Estimating coverage on precomputed bins " << nSamples << "\nmean\tvar\n " << mean << '\t' << var << '\n';
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
                vector<vector<double>> &clipCovCovTransP,                
                vector<vector<double>> &covSnvTransP,
                vector<vector<double>> &snvSnvTransP,
                int nCovStates, int nSNVStates,
                double diag, double offDiag,
                double beta, double clipBeta, double epsi12, double epsi23,
                vector<vector<double>> &emisP,
                int model, int maxCov, double mean, double var,
                vector<vector<vector<double>>> &binoP)
{
  covCovTransP.resize(nCovStates);
  clipCovCovTransP.resize(nCovStates);

  const double Diag1 = log(  1 -  (  exp(beta)   * (nCovStates-3)) + exp(epsi23) + exp(epsi12)  );

  const double Diag2 = log(  (1 -  (  exp(beta)   * (nCovStates-2)) /2)) ;

  const double Diag0 = log(  1 -  (  exp(beta)   * (nCovStates-2)) +  exp(epsi12)  );

  const double clDiag1 = log(  1 -  (  exp(clipBeta)   * (nCovStates-3)) + exp(epsi23) + exp(epsi12)  );

  const double clDiag2 = log(  (1 -  (  exp(clipBeta)   * (nCovStates-2)) /2)) ;
  
  const double clDiag0 = log(  1 -  (  exp(clipBeta)   * (nCovStates-2)) +  exp(epsi12)  );

  const double neutralScaler =  0 - log(6) ;

  const double stayScaler = neutralScaler + log(5);

  const double unif = log(1.0/nCovStates);

  const double large = -30;
  
  const double offDiagE = log(exp(large)/(nCovStates-1));

  const double DiagE = log(1- exp(large));// + log (nCovStates-1)) );
  
  //offDiag = unif;

  //diag = unif;

  for (int i=0;i<nCovStates;i++) {
    covCovTransP[i].resize(nCovStates);
    clipCovCovTransP[i].resize(nCovStates);
    for (int j=0;j<nCovStates;j++) {
      if (i==0)
      {//leaving del state
        if (i==j){
          covCovTransP[i][j] = DiagE;// - log(2) ;// Diag0 - log(2);
          clipCovCovTransP[i][j] = diag ;//- log(2) ; //clDiag0 - log(2);
        }
        else if(j==1){
          covCovTransP[i][j] = offDiagE ;// epsi12;
          clipCovCovTransP[i][j] = offDiag ; //epsi12;
        }
        else if(j==2){
          covCovTransP[i][j] = offDiagE;// - log(2);////Diag0 - log(2);
          clipCovCovTransP[i][j] =  offDiag ;//- log(2) ;//clDiag0 - log(2);
        }
        else{
          covCovTransP[i][j] = offDiagE ;// //beta;
          clipCovCovTransP[i][j] =  offDiag ;//clipBeta;
        }
      }

      else if(i==2)
      {//leaving neutral state
        if(i==j){
          covCovTransP[i][j] = DiagE;// -9.9999999999999994e-10;// DiagE ;//Diag1;
          clipCovCovTransP[i][j] =  diag ;//clDiag1;
        }
        else if(j==1){
          covCovTransP[i][j] = offDiagE;////epsi12;
          clipCovCovTransP[i][j] = offDiag ;//epsi12;
        }
        else if (j==3){
          covCovTransP[i][j] = offDiagE;////epsi23;
          clipCovCovTransP[i][j] = offDiag ;//epsi23;
        }
        else{
          covCovTransP[i][j] = offDiagE ;////beta; 
          clipCovCovTransP[i][j] =  offDiag;//clipBeta; 
        }
      }
      else
      {
        if(i==j){ 
          covCovTransP[i][j] = DiagE;// - log(2); //Diag2;
          clipCovCovTransP[i][j] =  diag ;//- log(2);//clDiag2; 
        }
        else if(j==2){
          covCovTransP[i][j] = offDiagE ;//- log(2);////Diag2;
          clipCovCovTransP[i][j] =  offDiag;// - log(2) ;//clDiag2;
        }
        else{
          covCovTransP[i][j] = offDiagE ;////beta;
          clipCovCovTransP[i][j] =  offDiag ;//clipBeta;
        }
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

// ------------
// Parameters
// ------------

Parameters::Parameters()
  : sampleName{"sample"}
  , CLI{"Hidden Markov Copy Number Caller\n"}
{
  //
  // Positional args
  //
  CLI.add_option("reference", referenceName,
    "Read reference from this FASTA file.")->
    type_name("FILE")->
    required();

  //
  // Input options
  //
  const std::string inputGroupName{"Input"};
  CLI.add_option("-a", bamFileName,
    "Read alignments from this BAM file and calculate depth on the fly.")->
    group(inputGroupName)->
    type_name("FILE");

  CLI.add_option("-b", covBedInFileName,
    "Read depth bed from this file (skip calculation of depth).")->
    group(inputGroupName)->
    type_name("FILE");

  CLI.add_option("-s", snvInFileName,
    "Read SNVs from this file (when not estimating from a BAM).")->
    group(inputGroupName)->
    type_name("FILE");

  CLI.add_option("-p", paramInFile,
    "Read parameter file (do not train with Baum-Welch).")->
    group(inputGroupName)->
    type_name("FILE");

  CLI.add_option("-l", clipInFileName,
    "Read clipping signature file (when not estimating from a BAM). ")->
    group(inputGroupName)->
    type_name("FILE");

  //
  // Depth calculation options
  //
  const std::string depthGroupName{"Depth Calculation"};
  CLI.add_option("-e",lepsi,
    "Value of log-epsilon. [-800]")->
    group(depthGroupName);

  CLI.add_option("-m", modelString,
    "Coverage model to use: Poisson (pois), or negative binomial (nb). [nb]")->
    group(depthGroupName);

  CLI.add_option("-t", nproc,
    "Number of threads. [4]")->
    group(depthGroupName);

  CLI.add_option("-c", useChrom,
    "Use this contig to estimate coverage. By default, longest contig.")->
    group(depthGroupName);

  //
  // Output options
  //
  const std::string outputGroupName{"Output"};
  CLI.add_option("-o", outFileName,
    "Output vcf to this file. Write to stdout if not provided.")->
    group(outputGroupName)->
    type_name("FILE");

  CLI.add_option("--sample", sampleName,
    "Sample name in the vcf ['sample']")->
    group(outputGroupName)->
    ignore_case();

  CLI.add_option("-M", mergeBins,
    "Merge consecutive bins with the same copy number.")->
    group(outputGroupName)->
    type_name("");

  CLI.add_option("-C", hmmChrom,
    "Only run hmm on this chrom.")->
    group(outputGroupName);

  CLI.add_option("-B", covBedOutFileName,
    "Write coverage bed to this file.")->
    group(outputGroupName)->
    type_name("FILE");

  CLI.add_option("-P", paramOutFile,
    "Write trained parameter file.(4 iterations)")->
    group(outputGroupName)->
    type_name("FILE");

  CLI.add_option("--readLength", averageReadLength,
		 "Set the average read length, special treatment of calls under this length.")->
     group(outputGroupName);

  CLI.add_option("-L", clipOutFileName,
    "Stores the number of reads with clipping > 500 bases in each bin.")->
    group(outputGroupName)->
    type_name("FILE");

  CLI.add_option("-S", snvOutFileName,
    "Write SNVs to this file.")->
    group(outputGroupName)->
    type_name("FILE");

  CLI.add_option("--writeFail", writeFail,
    "Write calls flagged as FAIL.")->
    group(outputGroupName)->type_name("");

  CLI.add_option("--bed", outBedName,
    "Output calls in bed format to this file.")->
    group(outputGroupName)->
    type_name("FILE");



  //
  // Post-parsing sanity checks
  //
  CLI.callback([this]() {
    if (this->modelString == "pois") {
      this->model = POIS;
    }
    if (this->covBedInFileName != "" and this->covBedOutFileName != "") {
      cerr << "ERROR. Cannot specify -b and -B.\n";
      exit(1);
    }
    if (this->covBedInFileName == "" and this->bamFileName == "") {
      cerr << "ERROR. Must specify either a coverage file or a bam file\n";
      exit(1);
    }
  });
}

// ------------

int hmcnc(Parameters& params) {

  double scale=2;
  int maxState=10;

  const string faiFileName{params.referenceName + ".fai"};
  vector<string> contigNames, allContigNames;
  vector<int>    contigLengths, allContigLengths;
  //
  // Determine what chroms to operate on.
  //
  ReadFai(faiFileName, allContigNames, allContigLengths);

  //  if (params.hmmChrom == "") {
  contigNames = allContigNames;
  contigLengths = allContigLengths;
  /*
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
  */

  vector<vector<int>> covBins, origCovBins;
  vector<vector<int>> clipBins;
  double mean;
  double var;
  int nStates;
  int maxCov;
  vector<double> startP;
  vector<vector<double>> covCovTransP, covSnvTransP, snvSnvTransP, clipCovCovTransP;
  vector<vector<double>> updateTransP, updateClipTransP;
  vector<vector<SNV>> snvs;
  vector<vector<int>> copyNumber;
  vector< vector<double>> fCov, bCov, fSNV, bSNV;
  vector<vector<double>> emisP;
  vector<vector<double>> updateEmisP;
  vector<vector<vector< double>>> binoP;
  vector<vector<double>> expCovCovTransP, expCovCovClipTransP, expCovSnvTransP, expSnvSnvTransP, expEmisP;
  vector<int> nReads;
  vector<long> totalBases;
  vector<double> averageCoverage;

  if (params.covBedInFileName != "") {
    if (params.snvInFileName == ""){
      cerr << "ERROR: SNV file (-s) must be specified if precomputed bins (-b) are used." << '\n';
      exit(0);
    }
    else{
    ReadCoverage(params.covBedInFileName, contigNames, covBins);
    ReadSNVs(params.snvInFileName, contigNames, snvs);
    }
  }

  if (params.clipInFileName != "") {
    ReadCoverage(params.clipInFileName, contigNames, clipBins);
  }

  if (params.paramInFile != "") {
     ReadParameterFile(params.paramInFile, nStates,
		       mean, var, maxState, maxCov,
		       startP, covCovTransP, clipCovCovTransP, emisP);
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
  std::vector<float> contigAvgCoverage(contigNames.size());
  pthread_mutex_t semaphore;
  pthread_mutex_init(&semaphore, NULL);

  int curSeq=0;

  //////////////////////////////////////////////

  vector<vector<Interval>> delT;
  delT.resize(contigNames.size());

  
  vector<vector<double>> Pcl, Pn;

  Pn.resize(contigNames.size());
  Pcl.resize(contigNames.size());

  vector<vector<Interval>> MdelT;
  MdelT.resize(contigNames.size());


  vector<vector<Interval>> copyIntervals;
  copyIntervals.resize(contigNames.size());


  vector<vector<Interval>> UnmergedNaiveIntervals;
  UnmergedNaiveIntervals.resize(contigNames.size());


  vector<vector<Interval>> mergedNaiveIntervals;
  mergedNaiveIntervals.resize(contigNames.size());

  vector<double> chromCopyNumber;
  chromCopyNumber.resize(contigNames.size());

  //////////////////////////////////////////////////

  nReads.resize(contigNames.size());
  totalBases.resize(contigNames.size());
  averageCoverage.resize(contigNames.size(), 0);

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
    threadInfo[procIndex].hmmChrom = params.hmmChrom;
    threadInfo[procIndex].transP = &covCovTransP;
    threadInfo[procIndex].clTransP = &clipCovCovTransP;
    threadInfo[procIndex].cl = &Pcl;
    threadInfo[procIndex].n = &Pn;
    threadInfo[procIndex].expTransP = &expCovCovTransP;
    threadInfo[procIndex].expClipTransP = &expCovCovClipTransP;
    threadInfo[procIndex].expEmisP = &expEmisP;
    threadInfo[procIndex].emisP = &emisP;
    threadInfo[procIndex].startP = &startP;
    threadInfo[procIndex].copyIntervals = &copyIntervals;
    threadInfo[procIndex].UnmergedNaiveIntervals = &UnmergedNaiveIntervals;
    threadInfo[procIndex].mergedNaiveIntervals = &mergedNaiveIntervals;
    threadInfo[procIndex].pModel=&pModel;
    threadInfo[procIndex].totalReads = &nReads;
    threadInfo[procIndex].totalBases = &totalBases;
    threadInfo[procIndex].averageCoverage = &averageCoverage;
    threadInfo[procIndex].chromCopyNumber = &chromCopyNumber;
    threadInfo[procIndex].delT = &delT;
  }

  //
  // If already trained, read those parameters.
  //
  copyNumber.resize(contigLengths.size());

  if (params.paramInFile == "") {
    //max cov value observed or upper cov bound -> max nState---------------
    nStates= std::min( maxState , MAX_CN  ) + 1; //+1 zeroth state
    //MAX_CN=nStates+1;
    startP.resize(nStates);
    for(int i=0;i<(nStates);i++) {
      if (i==2)
        startP[i]=log(0.5);
      else
        startP[i]=log(0.5) - log(nStates-1);
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
      Pn[c].resize(contigLengths[c]/BIN_LENGTH);
      Pcl[c].resize(contigLengths[c]/BIN_LENGTH);

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

    totalBaseSum=accumulate(totalBases.begin(), totalBases.end(), totalBaseSum);
    totalReadsSum=accumulate(nReads.begin(), nReads.end(), totalReadsSum);
    if (totalReadsSum > 0) {
      averageReadLength=totalBaseSum/totalReadsSum;
    }
    else {
      averageReadLength=0;
    }
    cerr << "Length cutoff of average read length " << averageReadLength << '\n';
    if (params.covBedOutFileName != "" ) {
      WriteCovBed(params.covBedOutFileName, contigNames, covBins);
    }
    if (params.clipOutFileName != "") {
      WriteClipBed(params.clipOutFileName, contigNames, clipBins, Pn, Pcl);
    }
    if (params.snvOutFileName != "") {
      WriteSNVs(params.snvOutFileName, contigNames, snvs);
    }
  }
  if (params.covBedInFileName != "") {
    StorePerChromAverageCoverage(covBins, averageCoverage);
  }

  EstimateCoverage(params.bamFileName, covBins, allContigNames, allContigLengths, params.useChrom, mean, var);

  //debug
  mean = 36;
  var = 100;

  if ((mean/var)>=0.90 and (mean/var)<=1.10){
    params.model= POIS;
    MODEL_TYPE model=POIS;
    std::cerr<<"Mean is approximately equal to variance, Model switched to Poisson."<<std::endl;
  }


  ChromCopyNumber(covBins, mean,chromCopyNumber);


  double clippingSum = 0;
  vector<int> clipCounts;
  for (auto c=0 ;c < contigNames.size(); c++) {
    if (chromCopyNumber[c] > 1.5 and chromCopyNumber[c] < 2.5) {
      NaiveCaller(covBins[c], UnmergedNaiveIntervals[c], mean );
      mergeNaiveIntervals(UnmergedNaiveIntervals[c], mergedNaiveIntervals[c], contigNames[c] );



    }
    else { 
      cerr << "Not using naive depth on " << contigNames[c] << " copy number " << chromCopyNumber[c] << endl;
    }
    for (int i=0; i < clipBins[c].size(); i++){
      if (clipBins[c][i]>0){
        clippingSum+=clipBins[c][i];
        clipCounts.push_back(clipBins[c][i]);
      }
    }
  }



  double clipMean, clipVar;
  double clipHmean = mean/2;
  bool   useClip;
  if (clipCounts.size() > 1) {    
    clipMean = (clippingSum/clipCounts.size());
    for (size_t i=0 ; i< clipCounts.size(); i++){
      clipVar += (clipCounts[i] - clipMean) * (clipCounts[i] - clipMean);
    }
    clipVar=(clipVar/(clipCounts.size()-1));
    useClip  = true;
  }
  else {
    clipMean=1;
    clipVar=3;
    useClip=false;
  }
    


  cerr<<"Clip Mean: "<<clipMean<<"\nClip Var: "<<clipVar<<"\nuseClip: "<<useClip<<endl;
  cerr<<"Cov Mean: "<<mean<<"\nCov Var: "<<var<<endl;


  //assert(clipMean*4<clipHmean);
  //
  // priors for clipping
  //
  // 100 clipping events WG
  // 10E7 bins 
  //

  double PNeutral = log((10E7 - 100)/10E7);
  double PClipped = log(100/10E7);

  cerr<<"WG Pn: "<<PNeutral<<"\nWG Pcl: "<<PClipped<<"\nEst. ~100 clips in 10E7 bins"<<endl;


  // 
  // update to NegB
  //

  poisson distributionClip(clipMean);
  poisson distributionHClip(clipHmean);

  int clip_count;
  double Pneutral, Pclipped, denom;
  int clipMax = MAX_CN * mean; // cap excessive clipping
  for (auto c=0 ;c < contigNames.size(); c++) {
    Pn[c].resize(clipBins[c].size());
    Pcl[c].resize(clipBins[c].size());
    for (auto b=0 ;b < clipBins[c].size(); b++) {
        clip_count = min(clipMax , clipBins[c][b]); 
        Pneutral = log(pdf(distributionClip, clip_count)) + PNeutral;
        Pclipped = log(pdf(distributionHClip, clip_count)) + PClipped;
        denom = PairSumOfLogP( Pneutral , Pclipped );

      	Pn[c][b]  = Pneutral - denom;
      	Pcl[c][b] = Pclipped - denom;

    }
  }

/*
  for (auto c=0 ;c < contigNames.size(); c++) {
    for (int i=0;i<Pcl[c].size();i++){
    std:cout<<Pn[c][i]<<"\t"<<Pcl[c][i]<<endl;
    }
  }
*/










  //
  // Cap coverage where hmm does not bother calculating.
  //
  maxCov=(int)mean/2*(maxState+1);
  origCovBins=covBins;
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

  const poisson distribution2(2*mean/2);
  const double result32=pdf(distribution2, 3*mean/2);
  const double result12=pdf(distribution2, mean/2);

  const poisson distribution3(3*mean/2);
  const poisson distribution1(mean/2);

  const double result33=pdf(distribution3, 3*mean/2);
  const double result11=pdf(distribution1, mean/2);

  const double epsi23 = log(result32)-log(result33);
  const double epsi12 = log(result12)-log(result11);


  const double lepsi23_nb = LgNegBinom(2, (int) (3*mean/2) , (float) (mean/2), (float)(var/2)  ) - LgNegBinom(3, (int) (3*mean/2) , (float) (mean/2), (float)(var/2)  );
  const double lepsi21_nb = LgNegBinom(2, (int) (mean/2) , (float) (mean/2), (float)(var/2)  ) - LgNegBinom(1, (int) (mean/2) , (float) (mean/2), (float)(var/2)  );



  const double eps = log(1);

  if( contigNames.size() > 1 or contigLengths[0] > MIN_CONTIG_SIZE ){
    vector<Interval> stats;
    double q = 0.99;
    int cn3quant=1;
    int cn1quant=1;
    double epsi21_emp=0;
    double epsi23_emp=0;
    double epsi3_emp=0;
    double epsi1_emp=0;
    //quant(mergedNaiveIntervals, 0.99, contigNames, stats);

    vector<Interval> lens1;
    vector<Interval> lens3;

    for (size_t c=0 ;c < contigNames.size(); c++) {
      for (int i=0; i < mergedNaiveIntervals[c].size(); i++){
        if (mergedNaiveIntervals[c][i].copyNumber==3){
          lens3.push_back(Interval( mergedNaiveIntervals[c][i].start, mergedNaiveIntervals[c][i].end, 3, (float)c, 0.0 ));
        }
        else if(mergedNaiveIntervals[c][i].copyNumber==1){
          lens1.push_back(Interval( mergedNaiveIntervals[c][i].start, mergedNaiveIntervals[c][i].end, 1, (float)c, 0.0 ));
        }
      }
    }
    if (!lens1.empty()){
      std::sort(lens1.begin(),lens1.end(),compareIntervalLength);
      const int idx1= (int)(q * lens1.size());
      const int idx1_5 = (int)(0.5 * lens1.size()); 
      stats.push_back(lens1[idx1]);
      cerr<<"Median CN=1 length: "<<lens1[idx1_5].end - lens1[idx1_5].start<<"\n";
      cerr<<"99th percentile CN=1 length: "<<stats[0].end-stats[0].start<<"\n";
      cn1quant = (int) (stats[0].end-stats[0].start)/BIN_LENGTH;
      const size_t contigg1 = (size_t) stats[0].averageCoverage;
      const int end1 = (int) stats[0].end/BIN_LENGTH;
      std::cerr<<contigNames[contigg1]<<"\t"<<end1-cn1quant<<"\t"<<end1<<std::endl;
      assert(UnmergedNaiveIntervals[contigg1].size()==covBins[contigg1].size());
      for (int i=end1-1; i>=end1-cn1quant;i-- ){
        assert(UnmergedNaiveIntervals[contigg1][i].averageCoverage > 0);
        epsi21_emp +=  LgNegBinom( 2 , (int) covBins[contigg1][i], (float) (mean/2), (float)(var/2) ) ;
        epsi1_emp +=  LgNegBinom( 1 , (int) covBins[contigg1][i], (float) (mean/2), (float)(var/2) ) ;
      }
      const double lepsi21_emp = epsi21_emp - epsi1_emp;
      std::cerr<<"empirical lepsi21: "<<lepsi21_emp<<std::endl;
    }
    else{
      cn1quant = 1;
      cerr<<"No CN=1 found in Naive intervals."<<endl;
    }

    if (!lens3.empty()){
      std::sort(lens3.begin(),lens3.end(),compareIntervalLength);
      const int idx3= (int)(q * lens3.size());
      const int idx3_5 = (int)(0.5 * lens1.size());
      stats.push_back(lens3[idx3]);
      cerr<<"Median CN=3 length: "<<lens3[idx3_5].end - lens3[idx3_5].start<<"\n";
      cerr<<"99th percentile CN=1 length: "<<stats[1].end-stats[1].start<<"\n";
      const int cn3quant = (int) (stats[1].end-stats[1].start)/BIN_LENGTH;
      const size_t contigg3 = (size_t) stats[1].averageCoverage;
      const int end3 = (int) stats[1].end/BIN_LENGTH;  
      std::cerr<<contigNames[contigg3]<<"\t"<<end3-cn3quant<<"\t"<<end3<<std::endl;
      assert(UnmergedNaiveIntervals[contigg3].size()==covBins[contigg3].size());

      for (int i=end3-1; i>=end3-cn3quant;i-- ){
        assert(UnmergedNaiveIntervals[contigg3][i].averageCoverage > (float)(mean) );
        epsi23_emp += LgNegBinom( 2 , (int) covBins[contigg3][i], (float) (mean/2), (float)(var/2) );
        epsi3_emp += LgNegBinom( 3 , (int) covBins[contigg3][i], (float) (mean/2), (float)(var/2) );
      }
      const double lepsi23_emp = epsi23_emp - epsi3_emp;
      std::cerr<<"empirical lepsi23: "<<lepsi23_emp<<endl;
    }
    else{
      cn3quant = 1;
      cerr<<"No CN=3 found in Naive intervals."<<endl;
    }






    const double scaler3 = (double) cn3quant;
    const double scaler1 = (double) cn1quant;
    std::cerr<<"scaler3: "<<scaler3<<"\tscaler1: "<<scaler1<<std::endl;
    const double beta = eps + (scaler3 * epsi23);  //log(nStates-1)
    const double beta_nb = eps + (scaler3 * lepsi23_nb);  //log(nStates-1)

    const double beta1 = eps + (scaler1 * epsi12);  //log(nStates-1)
    const double beta_nb1 = eps + (scaler1 * lepsi21_nb);  //log(nStates-1)


    std::cerr<<"beta_p3: "<<beta<<" beta_nb3: "<<beta_nb<<std::endl;
    std::cerr<<"beta_p1: "<<beta1<<" beta_nb1: "<<beta_nb1<<std::endl;

    /*
    //genome wide ratio

      for (size_t c=0 ;c < contigNames.size(); c++) {
        std::cerr<<contigNames[c]<<" first "<<LgNegBinom( 2 , (int) UnmergedNaiveIntervals[c][0].averageCoverage, (float) (mean/2), (float)(var/2) )<<"\t"<< UnmergedNaiveIntervals[c][0].averageCoverage<<std::endl;
        for (size_t i =0; i < UnmergedNaiveIntervals.size(); i++){
          if(UnmergedNaiveIntervals[c][i].copyNumber==3){
            assert(UnmergedNaiveIntervals[c][i].averageCoverage > (float)(mean) );
            epsi23_emp +=  LgNegBinom( 2 , (int) UnmergedNaiveIntervals[c][i].averageCoverage, (float) (mean/2), (float)(var/2)  );
            epsi3_emp +=  LgNegBinom( 3 , (int) UnmergedNaiveIntervals[c][i].averageCoverage, (float) (mean/2), (float)(var/2)  );

          }
          else if(UnmergedNaiveIntervals[c][i].copyNumber==1){
            assert(UnmergedNaiveIntervals[c][i].averageCoverage > 0);
            epsi21_emp +=  LgNegBinom( 2 , (int) UnmergedNaiveIntervals[c][i].averageCoverage, (float) (mean/2), (float)(var/2)  );
            epsi1_emp +=  LgNegBinom( 1 , (int) UnmergedNaiveIntervals[c][i].averageCoverage, (float) (mean/2), (float)(var/2)  );

          }
        }

      }
      */

      



  }





  const double beta_new = 100 * lepsi23_nb;

  const double clipBeta = 100 * lepsi23_nb;


  std::cerr<<"negBin lepsi23: "<<lepsi23_nb<<"\nnegBin lepsi21: "<<lepsi21_nb<<std::endl;
  std::cerr<<"poisson lepsi23: "<<epsi23<<"\npoisson lepsi21: "<<epsi12<<endl;

  std::cerr<<"Using neutral beta: "<<beta_new<<std::endl;
  std::cerr<<"Using clipped beta: "<<clipBeta<<std::endl;


  const int nSNVStates=3;
  const double unif=log(1.0/nStates);


  const double small=-5;

  if (params.paramInFile == "") {
    InitParams(covCovTransP, clipCovCovTransP, covSnvTransP, snvSnvTransP,
	       nStates, nSNVStates, 
         log(1-exp(small)), log(exp(small)/(nStates-1)),
         beta_new, clipBeta, 100*epsi12 , 100*epsi23,
	       emisP, params.model, maxCov, mean, var, binoP);
    
    cerr<<"\nNeutral"<<endl;
    printModel(covCovTransP, &cerr);
    cerr<<"\nClipped"<<endl;    
    printModel(clipCovCovTransP, &cerr);
    
    //    printEmissions(emisP);


    // Baum-Welch training.
    //

    double prevPX=0;
    assert(!emisP.empty());
    vector<vector<double> > prevTransP, prevClipTransP, prevEmisP;
    for (int i=0; i < 4; i++)
    {
      prevTransP = covCovTransP;
      prevEmisP = emisP;
      prevClipTransP = clipCovCovTransP;

      vector<double> stateWeightedTotCov(nStates, 0),
      stateWeightedTotVar(nStates, 0),
      stateTotProb(nStates, 0);
      vector<long> stateTotCov(nStates, 0), stateNCov(nStates, 0);

      expCovCovTransP.resize(nStates);
      expCovCovClipTransP.resize(nStates);
      expCovSnvTransP.resize(nStates);
      expSnvSnvTransP.resize(nStates);
      expEmisP.resize(nStates);
      for (int r=0; r < nStates; r++ ) {
        expCovCovTransP[r].resize(nStates,0);
        expCovCovClipTransP[r].resize(nStates,0);
        expCovSnvTransP[r].resize(nStates,0);
        expSnvSnvTransP[r].resize(nStates,0);
        fill(expCovCovTransP[r].begin(), expCovCovTransP[r].end(), 0);
        fill(expCovCovClipTransP[r].begin(), expCovCovClipTransP[r].end(), 0);
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
        clipCovCovTransP = prevClipTransP;
        emisP  = prevEmisP;
        break;
      }
      prevPX=px;

      vector<vector<double> > priorCovCov;
      priorCovCov.resize(covCovTransP.size());
      int nSites=covBins[0].size();

      //
      // Not used for now.
      //
      vector<vector<double> > prior;
      ApplyPriorToTransP(covBins,
			 covCovTransP.size(),
			 prior,
			 expCovCovTransP);
      ApplyPriorToClipTransP(covBins,
			 covCovTransP.size(),
			 prior,
			 expCovCovClipTransP);

      BaumWelchM(startP, covCovTransP, emisP, binoP,
        params.model,
        stateTotCov, stateNCov,
        expCovCovTransP, expCovCovClipTransP, 
        expEmisP,
        priorCovCov,
        updateTransP, updateClipTransP, updateEmisP);

      printModel(updateTransP, &cerr);
      printModel(updateClipTransP, &cerr);
      covCovTransP=updateTransP;
      clipCovCovTransP=updateClipTransP;


      if (params.paramOutFile != "") {
        string s = std::to_string(i);
        string outName = params.paramOutFile + "." + s;
        WriteParameterFile( outName , nStates, mean, var, maxState, maxCov, startP, covCovTransP, clipCovCovTransP, emisP);
      }


    }



    //
    // Eventually this needs to update for some multi-chrom code.
    //


  }
  else {
    //
    // Use parameters from file to compute copyIntervals
    //
    vector<vector<double>> f;
    vector<vector<double>> b;
    for (size_t i = 0; i < covBins.size(); ++i) {
      ForwardBackwards(startP, covCovTransP, clipCovCovTransP, emisP, covBins[i], f, b, Pn[i] , Pcl[i]);
      StorePosteriorMaxIntervals(covBins[i], f, b, copyIntervals[i]);
    }
  }

  //
  // Now filter cn=1 and cn=3 intervals, or short calls
  //
  // skip filtering of composite calls
  //
  //


  assert(copyIntervals.size() <= contigNames.size());
  assert(snvs.size() <= contigNames.size());
  AssignNearestClip(clipBins,
		    mean,
		    5,
		    contigNames,
		    copyIntervals);

  for (size_t c=0; c < contigNames.size(); c++) {
    int snvStart=0;
    for (size_t i=0; i < copyIntervals[c].size(); i++) {
      const int curCN = copyIntervals[c][i].copyNumber;
      
      if (i==0){
        if ( (copyIntervals[c][i+1].start - copyIntervals[c][i].end) < 2  ) //first call ith overlap with i+1 call
          continue;
      }
      if (i == copyIntervals[c].size()-1){
        if ( (copyIntervals[c][i].start - copyIntervals[c][i-1].end) < 2  ) //last call ith overlap with i-1 call
          continue;
      }
      if ( (copyIntervals[c][i].start - copyIntervals[c][i-1].end) < 2  and  (copyIntervals[c][i+1].start - copyIntervals[c][i].end) < 2 ) {
        continue;
      }

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
        copyIntervals[c][i].altInfo += ":BN";
        stringstream strm;
        strm << pCN-pCN2;
        copyIntervals[c][i].altSample += ":" + strm.str();
        if (pCN < pCN2) {
          copyIntervals[c][i].filter = "FAIL";
        }

        if (averageReadLength > 0 and copyIntervals[c][i].end-copyIntervals[c][i].start *2 < averageReadLength) {
          copyIntervals[c][i].filter = "FAIL";
        }

      	// Give it a chance to recover with clipping
      	if (copyIntervals[c][i].nFrontClip > MIN_FCLIP or copyIntervals[c][i].nEndClip > MIN_ECLIP) {
      	  copyIntervals[c][i].filter = "PASS";
      	}
        /*
      	if (copyIntervals[c][i].end - copyIntervals[c][i].start < MIN_SV_LENGTH) {
      	  copyIntervals[c][i].filter = "FAIL";
      	}
        */
	      snvStart=snvEnd;
      }
    }

 // Add gap-based deletion call from delT structure (sort, aggregate, consensus).
  // remove overlapping depth-based deletion call.

    mergeIntervals(delT[c], MdelT[c], contigNames[c] );

    intersectDelCall(MdelT[c], copyIntervals[c], mean );


  }
  //
  // Restore original copy number for high depth calls.
  //
  for (auto c=0; c < copyIntervals.size(); c++) {
    for (auto i=0; i < copyIntervals[c].size(); i++) {

      int intvStart = copyIntervals[c][i].start/100;
      int intvEnd   = copyIntervals[c][i].end/100;
      if (copyIntervals[c][i].copyNumber >= MAX_CN-2) {
		double intvCoverage = std::accumulate(&origCovBins[c][intvStart], &origCovBins[c][intvEnd], 0);
		double cnReal = 2*intvCoverage / ((intvEnd-intvStart)*mean);
		int cn=round(cnReal);
		//cerr << "Prev CN " << copyIntervals[c][i].copyNumber << " CUR: " << cn << endl;
		copyIntervals[c][i].copyNumber = cn;
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

  WriteVCF(*outPtr, params.referenceName, params.sampleName, contigNames, contigLengths, copyIntervals, writeFail);

  std::cerr<<"hmcnc done."<<endl;

  if (params.outBedName != "") {
      const string del_out = "del_cigar." + params.outBedName;    
      const string naive_out = "naive." + params.outBedName;
      WriteBed( copyIntervals, params.outBedName, contigNames);
      WriteBed( mergedNaiveIntervals, naive_out, contigNames);
      WriteBed( delT , del_out, contigNames);

    }


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
