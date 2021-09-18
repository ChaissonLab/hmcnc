#include "htslib/hts.h"
#include "htslib/sam.h"
#include <string>
#include <vector>
#include <iostream>
#include <ostream>
#include <stdlib.h>
#include <fstream>
#include <ctype.h>
#include <stdio.h>
#include <assert.h>
#include <algorithm>
#include "htslib/faidx.h"
#include <thread>
#include <map>
#include <sstream>
#include <iomanip>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/negative_binomial.hpp>
#include <boost/math/distributions/binomial.hpp>


#include <istream>
#include <limits>
#include <numeric>

using boost::math::binomial;

using boost::math::poisson;
using boost::math::pdf;
using boost::math::negative_binomial_distribution;
using std::vector;
using std::cout;
using std::endl;
using std::string;
using std::log;

int BIN_LENGTH=100;
int MIN_CLIP_LENGTH=500;
int MAX_CN=6;
float MISMAP_RATE=0.01;
using namespace std;
double lepsi=-800;
double minNeg=-1*numeric_limits<double>::epsilon();

void Moments(vector<double> &v, double &ex, double &var) {
  ex=0;
  for (int i=0; i < v.size(); i++) {
    ex+=v[i]*i;
  }
  var=0;
  for (int i=0; i < v.size(); i++) {
    var+= v[i]*(i-ex)*(i-ex);
  }
}


typedef enum { POIS, NEG_BINOM  } MODEL_TYPE;
class SNV {
public:
  char refNuc, altNuc;
  int ref;
  int alt;
  int pos;
  SNV() { assert(0);}
  SNV(int p, int r, int a, int rc, int ac) : pos(p), refNuc(r), altNuc(a), ref(rc), alt(ac) {}
  SNV(int p) : pos(p) { ref=0;alt=0;refNuc=0;altNuc=0;}
  int operator<(const SNV &rhs) const {
    return pos < rhs.pos;
  }
};

void Reset(vector<vector<double> > &v) {
  for (int i=0; i < v.size(); i++) {
    fill(v[i].begin(), v[i].end(), 0);
  }
}

double oned = 1.0;
double ONE = log(oned - 1e-10);
double epsilon = 0.00000000001;
double PairSumOfLogP(double a, double b) {
  double res = b;
  if (a != 0) {
    double m = max(a, b);

    assert(a <= 0);
    double diff = min(a,b) - m;
    double e = exp(diff);
    double lg=log(1+e);
    res      = m + lg;
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

double SumOfLogP(vector<double> &vals) {
  if (vals.size() == 0) {
    // Return 0 for error.
    return 0;
  }
  double maxVal = vals[0];
  for (size_t i=1; i < vals.size(); i++) { maxVal=max(maxVal,vals[i]); }
  double expSum=0;
  for (size_t i=0; i < vals.size(); i++) {
    expSum+= exp(vals[i]-maxVal);
  }
  return maxVal+log(expSum);
}

class Interval {
public:
  Interval(int s, int e, int cn, float avg, double p) : start(s), end(e), copyNumber(cn), averageCoverage(avg), pVal(p) { filter="PASS";}
  Interval() { assert(0);}
  int start, end, copyNumber;
  float averageCoverage;
  double pVal;
  string filter;
};

using namespace std;
class ThreadInfo {
public:
  htsFile *htsfp;
  hts_idx_t *bamidx;
  bam_hdr_t *samHeader;			
  faidx_t *fai;
  int *lastSeq;
  string refFileName;
  pthread_mutex_t *semaphore;
  vector<string> *contigNames;
  vector<int>    *contigLengths;
  vector<int> procChroms;
  vector<vector<int> > *covBins;
  vector<vector<int> > *clipBins;  
  vector<vector<SNV> > *snvs;
  vector<vector<int> > *copyNumber;
  vector<vector<double> > *transP, *emisP, *expTransP, *expEmisP;
  vector<double> *startP;
  vector<vector<Interval > > *copyIntervals;
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

//
// Not super 
void ReadCoverage(string &covFileName,
		  vector<string> &contigNames,
		  vector<vector<int> > &covBins) {
  ifstream covFile(covFileName.c_str());
  string chrom="", curChrom;
  int start, end;
  int cov;
  int last=0;
  int length;
  covFile.seekg(0, std::ios::end);    // go to the end
  length = covFile.tellg();           // report location (this is the length)
  covFile.seekg(0, std::ios::beg);    // go back to the beginning
  char *buffer = new char[length];    // allocate memory for a buffer of appropriate dimension
  covFile.read(buffer, length);       // read the whole file into the buffer
  cerr << "read cov buffer of len " << length << endl;
  covFile.close();  
  int i=0;
  string contigName("");
  int curContig=0;
  if (length > 0) {
    covBins.push_back(vector<int>() );
  }
  while (i < length) {
    while (i < length and isspace(buffer[i])) { i++; }
    int c=i;    
    while (i < length and isspace(buffer[i]) == false)  { i++; }
    if (i < length) {
      if (i-c > contigName.size()) { contigName.resize(i-c);}
      contigName = string(&buffer[c], i-c);
      ++i;
      start=atoi(&buffer[i]);
      while(i < length and buffer[i] != '\t') { i++;}
      i++;
      end=atoi(&buffer[i]);
      while(i < length and buffer[i] != '\t') { i++;}
      cov=atoi(&buffer[i]);
      i++;      
      //      sscanf(&buffer[i], "%d	%d	%d", &start, &end, &cov);
      if (contigName != contigNames[curContig]) {
	covBins.push_back(vector<int>());
	curContig++;
	cerr << "i " << i << "\t" << curContig << endl;	
      }
      covBins[curContig].push_back(cov);
    }
    while (i < length and buffer[i] != '\n') { i++;};
  }
}


void WriteSNVs(string &snvFileName,
	       vector<string> &contigNames,	       
	       vector<vector<SNV > > &snvs) {

  ofstream snvOut(snvFileName.c_str());
  for (int c=0; c < contigNames.size(); c++) {
    for (int i=0; i < snvs[c].size(); i++) {
      snvOut << contigNames[c] << "\t" << snvs[c][i].pos << "\t" << snvs[c][i].refNuc << "\t" << snvs[c][i].altNuc << "\t" << snvs[c][i].ref << "\t" << snvs[c][i].alt << endl;
    }
  }
}

void ReadSNVs(string &snvFileName,
	      vector<string> &contigNames,
	      vector<vector<SNV> > &snvs) {
  snvs.resize(contigNames.size());
  int curContig=0;
  ifstream snvIn(snvFileName);
  string line;
  string chrom;
  int pos, ref, alt;
  char refNuc, altNuc, t;
  
  while (curContig < contigNames.size()) {
    snvIn >> chrom >> pos >> refNuc >> altNuc >> ref >> alt;
    if (chrom == "" or snvIn.eof()) { break;}
    while (curContig < contigNames.size() and chrom != contigNames[curContig]) { curContig++;}
    if (curContig < contigNames.size()) {
      snvs[curContig].push_back(SNV(pos, refNuc, altNuc, ref, alt));
    }
  }
}

void WriteCovBed(string &covFileName,
		 vector<string> &contigNames,
		 vector<vector<int> > &covBins) {
  ofstream covFile(covFileName.c_str());
  for (int c=0; c < contigNames.size(); c++) {
    for (int i=0; i < covBins[c].size(); i++) {
      covFile << contigNames[c] << "\t" << i*100 << "\t" << (i+1)*100 << "\t" << covBins[c][i] << endl;
    }
  }
}

static void printModel(vector<vector<double> > &transP, ostream *strm)
{
  *strm << "\nTRANS: \n";
  for (int r=0;r<transP.size();r++)
    {
      *strm << r <<":";
      for (int c=0;c<transP[r].size();c++)
        {
	  *strm << "\t"<< transP[r][c];
        }
      *strm << endl;
    }
  *strm << endl;
}

static void printEmissions(vector<vector<double> > &emisP, ostream *strm)

{
  *strm << "\nEMIS: \n";
  for (int i=0; i < emisP[0].size(); i++) {
    *strm << std::setw(7) << i;
  }
  *strm << endl;
  for (size_t r=0;r<emisP.size(); r++) 
    {
      for (size_t c=0;c<emisP[r].size(); c++ )
	{
	  *strm << std::setw(7) << std::setprecision(2) << emisP[r][c];
	  if (c+1 < emisP[r].size()) 
	    {
	      *strm << "\t";
	    }	  
	}
      *strm << endl;
    }
}


double max_over_row(vector<vector<double> > &v , size_t col ,size_t nStates ){

  double maxi=-1 * (std::numeric_limits<double>::max()) ;
  for(size_t i=0;i< nStates;i++){
    maxi=std::max(v[i][col], maxi);
  }
  return maxi;
}


double max_over_rows(vector<vector<double> > &v , size_t col ,vector<vector<double> > &v2 , size_t nextState,size_t nStates ){
  double maxi2=-1 * (std::numeric_limits<double>::max()) ;

  for(size_t i=0;i< nStates;i++){
    maxi2=std::max( v[i][col] + v2[i][nextState] , maxi2);

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
      
  if (Hmean==0) {//no alignment in contig
    if (cn!= 0)
      result=lepsi;
    else
      result=0;    
  }
  else if (cov >= MAX_CN*Hmean) {//max_obs filtered previously
    if(cn!= MAX_CN)
      result=lepsi;
    else
      result=0;
  }
  else if(cn==0){//del_states
    //
    // Use poisson for 0-state assuming it's a random mismap process.
    poisson distribution(MISMAP_RATE*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0) 
      {
	result=lepsi;
      }
    else 
      {
	result=max(lepsi, log(prob));
      }
  }
  else{
    Hmean*=cn;
    Hvar*=cn;
    p=Hmean/Hvar;
    r=Hmean*(Hmean/Hvar)/(1-Hmean/Hvar);
    
    negative_binomial_distribution<double> distribution(r,p);    
    double prob=pdf(distribution, cov);
    if (prob == 0) 
      {
	result=lepsi;
      }
    else
      {
	result=max(lepsi, log(prob));

      }
  }
  return result;  
}

double LgBinom(double p, int s, int n) {
  binomial snv(n, p);
  double pVal;
  pVal=pdf(snv,s);
  if (pVal == 0) {
    return lepsi;
  }
  else {
    double lP=log(pVal);
    double retVal=max(lP, lepsi);
    assert(isnan(lP) == 0);
    return retVal;
  }

}

double LgPrpoiss(int cn,  int cov, int Hmean) {
  double result=0;

  // Boundary case
  if (Hmean==0){//no alignment in contig
    if(cn!= 0)
      // Cannot emit 0 reads from non-zero states.
      result=lepsi;
    else
      // Can only emit 0
      result=0;
  }
  else if (cov >= MAX_CN*Hmean){//max_obs filtered previously
    //
    // Have state that catches ultra-high probability states.
    if(cn != MAX_CN)
      result=lepsi;
    else
      result=0;
  }
  else if(cn==0){//del_states
    poisson distribution(MISMAP_RATE*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0) {
      result = lepsi;      
    }
    else {
      result = max(lepsi, log(prob));
    }
  }
  else {
    poisson distribution(cn*Hmean);
    double prob=pdf(distribution, cov);
    if (prob == 0)
      {
	result=lepsi;
      }
    else 
      {
	result=max(lepsi, log(prob));	
      }
  }
  return result;
}




static void correctModel(vector<vector<double> > &transP,
                         int nStates)
{
  double sum;
  for (int i=0;i<nStates;i++)
    {
      sum = 0;
      for (int j=0;j<nStates;j++)
	sum+= std::exp(transP[i][j]);
      for (int j=0;j<nStates;j++)
	transP[i][j]= log(std::exp(transP[i][j])/sum);
    }
}//correctModel

void CombineEmissions( vector<int> &obs,
		       vector<SNV> &snvs,
		       vector<uint8_t> &isCov,
		       vector<int> &obsIndex) {
  int totObs=obs.size() + snvs.size();
  isCov.resize(totObs);
  fill(isCov.begin(), isCov.end(), false);
  obsIndex.resize(totObs);
  fill(obsIndex.begin(), obsIndex.end(), -1);
  int c=0,s=0,p=0, obsi=0;
  while (c < obs.size() or s < snvs.size()) {
    int curPos=c*100;
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
    cerr << "ERROR computing obs incex." << endl;
    exit(1);
  }
}


double CSEmisP(	int state,
		int pos,
		vector<int> &obs,
		vector<SNV> &snvs,
		vector<uint8_t> &isCov,
		vector<int> &obsIndex,
		vector<vector<double> > &emisP, 
		vector<vector<vector<double> > > &binoP ) {
  if (isCov[pos]) {
    int covIndex=obsIndex[pos];
    return emisP[state][obs[covIndex]];
  }
  else {
    int snvIndex = obsIndex[pos];
    int maxCov = emisP[0].size();
    int ref=snvs[snvIndex].ref;
    int alt=snvs[snvIndex].alt;
    int nucCov = ref+alt;
    if (nucCov > maxCov) {
      double f=((double)maxCov)/nucCov;
      ref=ref*f;
      alt=alt*f;
      nucCov=ref+alt;
    }
    assert(nucCov < binoP[state].size());
    int m=min(ref,alt);
    assert(m < binoP[state][nucCov].size());
    return binoP[state][nucCov][m];
  }
}  

		      

double ForwardBackwards( vector<double> &startP,
			 vector<vector<double> > &covCovTransP,			 
			 vector<vector<double> > &emisP,
			 vector<int> &obs,
			 vector<vector<double> > &f,
			 vector<vector<double> > &b) {

  
  int totObs    = obs.size();
  int nCovObs   = obs.size();
  int nCovStates = startP.size();

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
  
  double lgthird=log(1/3.);
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

double BaumWelchEOnChrom(vector<double> &startP,
			 vector<vector<double> > &covCovTransP,			 
			 vector<vector<double> > &emisP,
			 vector<int> &obs,
			 vector<vector<double> > &f,
			 vector<vector<double> > &b,
			 vector<vector<double> > &expCovCovTransP,
			 vector<vector<double> > &expEmisP) {
  int nStates = startP.size();
  int nObs    = obs.size();
  double px;
  px = ForwardBackwards( startP, covCovTransP, emisP,
			 obs,
			 f, b);
  /*
  ofstream fb("fb.tsv");
  for (int k=0; k < f[0].size()-1; k++) {
    fb << k << "\t";
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
      if (j +1 < nStates) { fb << "\t";}
      if (maxfb < p) {
	maxfb=p;
	maxi=j;
      }
    }
    fb << "\t" << maxi << "\t" << maxfb/fbSum << "\t" << obs[k] << endl;

  }
  fb.close();
*/
    
  for (int k=1; k< nObs-1; k++) {
    double logSum=0;
    //
    // Calculate total probability of all transitions at this step.
    //
    for (int i=0; i < covCovTransP.size(); i++) {
      covCovTransP[i].resize(covCovTransP[i].size());
      for (int j=0; j < covCovTransP[i].size(); j++) {
	logSum = PairSumOfLogP(logSum, f[i][k] + covCovTransP[i][j] + emisP[j][obs[k+1]] + b[j][k+1]);
      }
    }
    for (int i=0; i < covCovTransP.size(); i++) {
      for (int j=0; j < covCovTransP[i].size(); j++) {
	double pEdge=f[i][k] + covCovTransP[i][j] + emisP[j][obs[k+1]] + b[j][k+1];
	assert(isnan(exp(pEdge-logSum)) == false);
	expCovCovTransP[i][j] += exp(pEdge - logSum);	  
      }
    }
  }
  for (int ei=0; ei < expEmisP.size(); ei++) {
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
  for (int ei=0; ei < expEmisP.size(); ei++) {
    double rowSum=0;
    rowSum=std::accumulate(expEmisP[ei].begin(),expEmisP[ei].end(), rowSum);
    for (int pi=0; pi < expEmisP[ei].size(); pi++) {
      expEmisP[ei][pi] /= rowSum;
    }
  }
  
  return px;
  
}

void PrintIntervals(string chrom,
		    vector<Interval> &intv,
		    ostream &out) {
  for (int i=0; i < intv.size(); i++) {
    out << chrom << "\t" << intv[i].start << "\t" << intv[i].end << "\t" << intv[i].copyNumber << "\t" << intv[i].averageCoverage << "\t" << intv[i].pVal << endl;
  }
}
void StorePosteriorMaxIntervals(vector<int> &cov,
				vector<vector<double> > &f,
				vector<vector<double> > &b,
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
  for (i=0; i < f[0].size()-1; i++) {
    double colMax=f[0][i] + b[0][i];
    int colMaxCN=0;
    assert(i < cov.size());
    totCov+=cov[i];
    colSum=f[0][i] + b[0][i];
    for (int j = 1; j < f.size(); j++) {
      double v=f[j][i] + b[j][i];
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
      //      out << chrom << "\t" << prevStart*100 << "\t" << i*100 << "\t" << prevCN << "\t" << i-prevStart << "\t" << totCov/(i-prevStart) << "\t" << nSNV << endl;
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
string version="0.8";
string reference;
void WriteVCF(ostream &out,
	      string &refName,
	      string &sampleName,
	      vector<string> &contigNames,
	      vector<int> &contigLengths,
	      vector<vector<Interval> > &intervals) {
  out << "##fileformat=VCFv4.1" << endl
      << "##source=hmmcnc_v" << version << endl
      << "##reference=" << reference << endl;
  for (int i = 0; i < contigNames.size(); i++) {
    out << "##contig=<ID=" << contigNames[i] << ",length=" << contigLengths[i]
        << ">" << endl;
  }

  out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of "
    "structural variant\">"
      << endl
      << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of "
    "the structural variant described in this record\">"
      << endl
      << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in "
    "length between REF and ALT alleles\">"
      << endl
      << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise "
    "structural variation\">"
      << endl;
  out << "##FORMAT=<ID=CN,Number=1,Type=String,Description=\"CopyNumber\">" 
      << endl
      << "##FORMAT=<ID=PP,Number=R,Type=Float,Description=\"Relative posterior "
    "probability (phred)\">"
      << endl
      << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at "
    "this position for this sample\">"
      << endl
      << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName
      << endl;
  for (int c = 0; c < contigNames.size(); c++) {
    for (int i = 0; i < intervals[c].size(); i++) {
      if (intervals[c][i].copyNumber != 2) {

        std::string cntype = (intervals[c][i].copyNumber > 2) ? "DUP" : "DEL";

        out << contigNames[c] << "\t" << intervals[c][i].start
            << "\t.\t<CNV>\t<CNV>\t30\t" << intervals[c][i].filter << "\t"
            << "SVTYPE=" << cntype << ";"
            << "END=" << intervals[c][i].end
            << ";SVLEN=" << intervals[c][i].end - intervals[c][i].start
            << ";IMPRECISE\t"
            << "CN:PP:DP\t" << intervals[c][i].copyNumber << ":"
            << intervals[c][i].pVal << ":" << intervals[c][i].averageCoverage
            << endl;
      }
    }
  }  
}
		       

void UpdateEmisP(vector<vector<double> > &emisP,
		 vector<vector<double> > &expEmisP,
		 int model=NEG_BINOM) {
  int nCovStates=emisP.size();
  for (int i=1;i<nCovStates;i++) {
    double mean, var;
    Moments(expEmisP[i], mean, var);
    int maxCov=expEmisP[i].size()-1;
    double stateSum=0;
    emisP[i].resize(expEmisP[i].size());
    for (int j=0;j<=maxCov;j++) {
      if (model == POIS or mean > var) {	
	emisP[i][j]=LgPrpoiss( (int) i , j , (int) mean );
	stateSum+=exp(LgPrpoiss( (int) i , j , (int) mean ));
      }
      else  {
	emisP[i][j]=LgNegBinom((int)i, (int) j, mean, var);
	stateSum+=exp(LgNegBinom((int)i, (int) j, mean, var));	
      }      
    }
  }
}

void BaumWelchM(vector<double> &startP,
		vector<vector<double> > &transP,
		vector<vector<double> > &emisP,
		vector<vector<vector< double> > > &binoP,
		int model,
		vector<long> &stateTotCov,
		vector<long> &stateNCov,
		vector<vector<double> > &expTransP,
		vector<vector<double> > &expEmisP,
		vector<vector<double> > &covCovPrior,
		vector<vector<double> > &updateTransP,
		vector<vector<double> > &updateEmisP) {

  //
  // M step.
  //
  int nStates=startP.size();
  updateTransP.resize(nStates);
  vector<double> colSums;

  cerr << "Update trans: " << endl;
  cerr << "p\t";
  for (int j=0; j < nStates; j++) {
    cerr << std::setw(8) << j << "\t";
  }
  cerr << endl;
  for (int j=0; j < nStates; j++) {
    double colSum=0;
    updateTransP[j].resize(nStates);
    for (int k=0; k< nStates; k++) {
      colSum +=expTransP[j][k]; //PairSumOfLogP(colSum, expTransP[j][k]);
    }
    colSums.push_back(colSum);
    cerr << j;
    for (int k=0; k < nStates; k++) {
      updateTransP[j][k] = log(expTransP[j][k]/colSum); //min(ONE, expTransP[j][k] - colSum);
      cerr << "\t" << std::setw(8) << updateTransP[j][k];      
    }
    cerr << endl;
  }
  //
  // M step for emissions -- use summary statistics to recompute emission values
  // Eventually to really be BW, this should fit a negative binomial dist to the data.
  //
  updateEmisP.resize(nStates);  

  vector<double> stateMean(nStates);
  for (int i=0; i < nStates; i++) {
    if (stateNCov[i] > 0) {
      stateMean[i]=stateTotCov[i]/stateNCov[i];
    }
    else {
      stateMean[i] = 0;
    }
  }
  int maxCov=emisP[0].size();
  UpdateEmisP(updateEmisP, expEmisP, model);
}
		
  

void viterbi( vector<double> &startP,
	      vector<vector<double> > &covCovTransP,
	      vector<vector<double> > &covSnvTransP,
	      vector<vector<double> > &snvSnvTransP,
	      vector<vector<double> > &emisP,
	      vector<vector<vector<double> > > &binoP,
	      vector<int> &cov,
	      vector<SNV> &snvs,
	      vector<uint8_t> &isCov,
	      vector<int> &obsIndex,
	      vector<int> & viterbiPath) {    
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
  int operator<(const CountNuc &rhs) {
    return count < rhs.count;
  }
};


int StoreSNVs(char *contigSeq,
	      int contigLength,
	      float mean,
	      vector<int> &nA, vector<int> &nC, vector<int> &nG, vector<int> &nT, vector<int> &nDel,
	      vector<SNV> &snvs) {

  vector<int* > fPtr(5);
  //
  // Easier for typing
  //

  const char *nucs="ACGTd";
  fPtr[0] = &nA[0];
  fPtr[1] = &nC[0];
  fPtr[2] = &nG[0];
  fPtr[3] = &nT[0];
  fPtr[4] = &nDel[0];
  vector<CountNuc > counts;
  counts.resize(5);
  long total=0;
  if (mean==0) {
    for (int i=0; i < contigLength; i++) {
      for (int j=0; j < fPtr.size(); j++) {
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
    char refNuc=toupper(contigSeq[i]);
    if (counts[4].index != 4 and
	counts[3].index != 4 and
	refNuc != 'N' and
	counts[3].count > 0.25*mean and
	counts[4].count > 0.25*mean and
	counts[3].count < 2*mean )  {
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
	cerr << "Stored " << snvs.size() << " at " << i << endl;
      }      
    }
  }
  return 1;
}



int IncrementCounts(bam1_t *b,
		    int contigLength,
		     vector<int> &nA, vector<int> &nC, vector<int> &nG, vector<int> &nT, vector<int> &nDel) {
  int readLength = b->core.l_qseq;			
  if (readLength < BIN_LENGTH or b->core.qual < 10 or b->core.flag & 0x800) {
    return 0;
  }
	  
  vector<char> seq(readLength);
  uint8_t *q = bam_get_seq(b);
  for (int i=0; i < readLength; i++) {seq[i]=seq_nt16_str[bam_seqi(q,i)];	}
  uint32_t* cigar = bam_get_cigar(b);
  int refLen = bam_cigar2rlen(b->core.n_cigar, cigar);
  int qPos=0;
  int refPos = b->core.pos;
  int ci;
  int regionOffset=refPos;
  bool first=true;
  for (ci=0; ci < b->core.n_cigar; ci++) {
    int opLen=bam_cigar_oplen(cigar[ci]);
    int op=bam_cigar_op(cigar[ci]);
	  
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
	  if (refPos >= contigLength){ 
	    break;
	  }
	  if (refPos >= 1) {
	    first=false;
	    char nuc=toupper(seq[qPos]);
	    assert(regionOffset < nA.size());
	    /*
	    if (regionOffset == 20504515 or refPos == 20504515 ) {
	      cout << regionOffset << "\t" << nuc << "\t" << nC[regionOffset] << "\t" << nT[regionOffset] << "\t" << bam_get_qname(b) << endl;
	      }*/
	    if (nuc == 'A') { nA[regionOffset]++;}
	    if (nuc == 'C') { nC[regionOffset]++;}
	    if (nuc == 'G') { nG[regionOffset]++;}
	    if (nuc == 'T') { nT[regionOffset]++;}
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

    int curSeq = *((*threadInfo).lastSeq);
    *(threadInfo->lastSeq) = *(threadInfo->lastSeq) + 1;
    pthread_mutex_unlock(threadInfo->semaphore);      
    
    if (curSeq >= threadInfo->contigNames->size()) {
      break;
    }
    vector<vector<double > > f,b, expCovCovTransP, expEmisP;
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
    for (int i=0; i < threadInfo->transP->size(); i++) {
      for (int j=0; j < (*threadInfo->transP)[i].size(); j++) {
	(*threadInfo->expTransP)[i][j] += expCovCovTransP[i][j];
      }
    }
    for (int i=0; i < threadInfo->emisP->size(); i++) {
      for (int j=0; j < (*threadInfo->emisP)[i].size(); j++) {
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

    int curSeq = *((*threadInfo).lastSeq);
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
    
    int contigLength=threadInfo->contigLengths->at(curSeq);


    vector<int> nA(contigLength, 0), nC(contigLength, 0), nT(contigLength, 0), nG(contigLength,0), nDel(contigLength, 0);
    
    stringstream regionStrm;
    regionStrm << (*(*threadInfo).contigNames)[curSeq];// << ":1-" << contigLength;

    string region=regionStrm.str();


    hts_itr_t *regionIter = sam_itr_querys(threadInfo->bamidx, threadInfo->samHeader, region.c_str());
    int chromLen;
    char *chromSeq = fai_fetch(threadInfo->fai, region.c_str(), &chromLen);
    
    bool continueParsing=true;
    vector<bam1_t*> reads; //(bam_init1());
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
	bam1_t *b = bam_init1();	
	int res=sam_itr_next(threadInfo->htsfp, regionIter, b);
	bufSize+= b->l_data;
	totalSize+= b->l_data;
	
	if (res < 0) { // or totalReads < 15000) {
	  continueParsing = false;
	  cerr << "Ending parsing of " << region << " with " << totalSize << " data and " << chunkNumber << " iterations." << endl;
	  break;
	}
	/*
	cout << "read " << bam_get_qname(b) << endl;
	if (strcmp("m64043_200714_124814/138021154/ccs", bam_get_qname(b)) == 0) {
	  cout << "poblem" << endl;
	  }*/
	endpos=bam_endpos(b);
	reads.push_back(b);
	++totalReads;
      }
      cerr << "Reading " << (*threadInfo->contigNames)[curSeq] << ", chunk " << chunkNumber << ".\t" << reads.size() << "/" << totalReads << " reads/total" << endl;
      ++chunkNumber;
      pthread_mutex_unlock(threadInfo->semaphore);

      for (int readIndex=0; readIndex < reads.size(); readIndex++) {
	bam1_t *b=reads[readIndex];
	IncrementCounts(b, contigLength, nA, nC, nG, nT, nDel);
	endpos=bam_endpos(b);
	startpos=b->core.pos;
	(*(*threadInfo).totalReads)[curSeq]++;
	(*(*threadInfo).totalBases)[curSeq]+= endpos-startpos;
	
	int nCigar=b->core.n_cigar;
	uint32_t*cigar=bam_get_cigar(b);
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
	  int bin=startpos/BIN_LENGTH;
	  (*threadInfo->clipBins)[curSeq][bin] += 1;
	}
	if (backClipLen > MIN_CLIP_LENGTH) {
	  int bin=endpos/BIN_LENGTH;
	  (*threadInfo->clipBins)[curSeq][bin] += 1;	  	  
	}
	bam_destroy1(b);
      }
    }
    // Never compute in the last bin
    int nBins=contigLength/BIN_LENGTH;

    for (int bin=0; bin < nBins; bin++) {
      int start=bin*BIN_LENGTH;
      int end=min((bin+1)*BIN_LENGTH, contigLength);
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
    cerr << "Stored " << (*threadInfo->snvs)[curSeq].size() << " snvs for " << (*threadInfo->contigNames)[curSeq] << endl;
    
  }
  if (threadInfo->exit) {
    pthread_exit(NULL);
  }
}
    
vector<int> NucMap;
int GetRefAndAlt(char refNuc, vector<int> &counts, int &ref, int &alt) {
  vector<int> sCounts=counts;
  sort(sCounts.begin(), sCounts.end());
  
  int second = sCounts[2];
  int first  = sCounts[3];
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
    ref=0; alt=0;
    return 4;
  }
  int refNucIndex=NucMap[refNuc];
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
int EstimateCoverage(string &bamFileName, vector<vector<int> > &allCovBins, vector<string> &chroms, vector<int> &lengths, string &useChrom, double &mean, double &var) {
  int useChromIndex=0;
  if (useChrom == "") {
    int maxLen=0;
    for (int i=0; i < lengths.size(); i++) {
      long totCov=0;
      for (int j=0; j < allCovBins[i].size(); j++ ) {
	totCov+=allCovBins[i][j];
      }
      if (totCov > 0 and lengths[i] > maxLen) {
	useChrom = chroms[i];
	maxLen=lengths[i];
	useChromIndex=i;
      }
    }
  }
  cerr << "Estimating coverage from " << useChrom << endl;
  int contigLength=0;
  for (int i=0; i < chroms.size(); i++) {
    if (chroms[i] == useChrom) {
      contigLength=lengths[i];
      useChromIndex=i;
      break;
    }
  }
  
  if (contigLength == 0) {
    cerr << "ERROR Could not estimate coverage." <<endl;
    exit(1);
  }
  if (allCovBins.size() > 0) {
    assert(allCovBins[useChromIndex].size() == lengths[useChromIndex]/100);
    int lastBin=allCovBins[useChromIndex].size();
    if (lastBin == 0) {
      cerr << "ERROR. Could not estimate coverage using precomputed bins." << endl;
      exit(1);
    }
    long totCov=0;
    for (int binIndex=0; binIndex<lastBin; binIndex++) {
      totCov+=allCovBins[useChromIndex][binIndex];
    }
    mean=totCov/lastBin;
    //
    // Recompute summary stats using limited data
    int nSamples=0;
    totCov=0;
    long totCovSq=0;
    for (int binIndex=0; binIndex<lastBin; binIndex++) {
      if (allCovBins[useChromIndex][binIndex] > 0.25*mean and allCovBins[useChromIndex][binIndex] < 1.75 * mean) {
	totCov+=allCovBins[useChromIndex][binIndex];
	totCovSq+=allCovBins[useChromIndex][binIndex]*allCovBins[useChromIndex][binIndex];
	nSamples++;
      }
    }
    if (nSamples > 0) {
      mean=totCov/nSamples;
      var=totCovSq/nSamples-mean*mean;
      cerr << "Estimating coverage on precomputed bins " << nSamples << " ending at " << mean << "\t" << var << endl;
    }
    else {
      cerr << "Could not estimate coverage using precomputed coverage on chrom " << useChrom << endl;
      exit(1);
    }
  }
  else {
    htsFile *htsfp;
    vector<int> covBins;
    htsfp = hts_open(bamFileName.c_str(),"r");
  
    hts_idx_t *bamidx;
    if ((bamidx = sam_index_load(htsfp, bamFileName.c_str())) == 0) {
      cerr << "ERROR reading index" << endl;
      exit(0);
    }
    
    const htsFormat *fmt = hts_get_format(htsfp);
    if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
      cerr << "Cannot determine format of input reads." << endl;
      exit(1);
    }

    bam_hdr_t *samHeader;			
    samHeader = sam_hdr_read(htsfp);
    
    hts_itr_t *regionIter = sam_itr_querys(bamidx, samHeader, useChrom.c_str());

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
	bam1_t *b = bam_init1();	
	int res=sam_itr_next(htsfp, regionIter, b);
	bufSize+= b->l_data;
	totalSize+= b->l_data;
	
	if (res < 0) {
	  continueParsing = false;
	  break;
	}
	if (IncrementCounts(b, contigLength, nA, nC, nG, nT, nDel)) {
	  curEndPos=bam_endpos(b);
	}
	bam_destroy1(b);
	++nReads;
      }

      //
      // Compute coverage for bins
      int lastBin=(curEndPos-30000)/100;
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
	cerr << "Estimating coverage " << nReads << " ending at " << curEndPos << "\t" << mean << "\t" << var << endl;
      }
      if (nSamples > 80000) {
	return 1;
      }
    }
  }
  return 1;
}



void WriteParameterFile(string &fileName,
			int nStates,
			double covMean,
			double covVar,
			int maxState,
			int maxCov,
			vector<double > &startP,
			vector<vector<double > >&transP,
			vector<vector<double > >&emisP) {
  ofstream outFile(fileName.c_str());
  outFile << "nStates\t" << nStates << endl
	  << "covMean\t" << covMean << endl
	  << "covVar\t" << covVar  << endl
	  << "maxState\t" << maxState << endl
	  << "maxCov\t" << maxCov << endl
	  << "startP" << endl;
  for (int i=0; i < startP.size(); i++) {
    outFile << startP[i];
    if (i+1 < startP.size()) { outFile << "\t";}
    outFile << endl;
  }
  outFile << "transP\t" << transP.size() << "\t" << transP[0].size() << endl;
  for (int i=0; i < transP.size(); i++) {
    for (int j=0; j < transP[j].size(); j++) {
      outFile << transP[i][j];
      if (i+1 < transP.size()) { outFile << "\t";}
    }
    outFile << endl;
  }
  outFile << "emisP\t" << emisP.size() << "\t" << emisP[0].size() << endl;
  for (int i=0; i < emisP.size(); i++) {
    for (int j=0; j < emisP[j].size(); j++) {
      outFile << emisP[i][j];
      if (i+1 < emisP.size()) { outFile << "\t";}
    }
    outFile << endl;
  }
}  
  

void ReadParameterFile(string &fileName,
		       int &nStates,
		       double &covMean,
		       double &covVar,
		       int &maxState,
		       int &maxCov,
		       vector<double > &startP,
		       vector<vector<double > >&transP,
		       vector<vector<double > >&emisP) {
  ifstream inFile(fileName.c_str());
  string spacer;
  inFile >> spacer >> nStates;
  inFile >> spacer >>covMean;
  inFile >> spacer >> covVar;
  inFile >> spacer >> maxState;
  inFile >> spacer >> maxCov;
  inFile >> spacer;
  double val;
  for (int i=0; i < nStates; i++) {
    inFile >> val;
    startP.push_back(val);
  }
  int nr, nc;
  inFile >> spacer >> nr >> nc;
  transP.resize(nr);
  for (int i=0; i < nr; i++) {
    for (int j=0; j < nc; j++) {
      inFile>> val;
      transP[i].push_back(val);
    }
  }
  inFile >> spacer >> nr >> nc;
  emisP.resize(nr);
  for (int i=0; i < nr; i++) {    
    for (int j=0; j < nc; j++) {
      inFile >> val;
      emisP[i].push_back(val);
    }
  }
}  

    
  
  

  
void ReadFai(string faiFileName, vector<string> &contigNames, vector<int> &contigLengths)   {
  ifstream faiIn(faiFileName.c_str());
  if (faiIn.good() == false) {
    cerr << "ERROR. Reference is not indexed, or could not open .fai file" << endl;
    exit(1);
  }
    
  while (faiIn) {
    string line;
    getline(faiIn, line);
    stringstream strm(line);
    if (line != "") {
      string contig;
      int length;
	
      strm >> contig;
      strm >> length;      
      contigNames.push_back(contig);
      contigLengths.push_back(length);
    }
  }
}


void InitParams(vector<vector<double> > &covCovTransP,
		vector<vector<double> > &covSnvTransP,		
		vector<vector<double> > &snvSnvTransP,
		int nCovStates, int nSNVStates, double diag, double offDiag,
		vector<vector<double> > &emisP, int model, int maxCov, double mean, double var,
		vector<vector<vector<double> > > &binoP) {
  covCovTransP.resize(nCovStates);
  for (size_t i=0;i<nCovStates;i++) {
    covCovTransP[i].resize(nCovStates);
    for (size_t j=0;j<nCovStates;j++) {
      if(i==j) {
	covCovTransP[i][j]  = diag; //log(1 - std::exp(beta) );
      }
      else {
	covCovTransP[i][j] = offDiag;
      }
    }
  }
  double snvOffDiag=(1-diag)/2;
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
  for (size_t i=0;i<nCovStates;i++) {
    emisP[i].resize(maxCov+1);
    double stateSum=0;
    for (size_t j=0;j<=maxCov;j++) {
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
    for (size_t j=0; j <= maxCov; j++) {
      binoP[i][j].resize(j+1);
    }
  }
    
  //
  // Initialize copy number 1
  int i=0;
  double bino_one=0.9;
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

void PrintHelp() {
    cerr << "usage: hmmcnc reference.fa" << endl
         << "   -a alignments    Read alignments from this file and calculate depth on the fly." << endl
         << "   -b bed           Read depth bed from this file. Skip calculation of depth." << endl
         << "   -s snv-file      Read SNVs from this file (when not estimating from a BAM)" << endl
         << "   -p parameter     Read parameter file (do not train with Baum-Welch)" << endl
	 << " Options controlling depth calculation " << endl      
	 << "   -e value (float)   Value of log-epsilon (-500)." << endl
	 << "   -m value [pois|nb] Coverage model to use, Poisson (pois), or negative binomial (nb). Default nb." << endl
	 << "   -x value Max state to allow (10)" << endl
	 << " -t value (int)     Number of threads (4) " << endl            
	 << " -c contig          Use this contig to estimate coverage. By default, longest contig." << endl      
	 << " Options controlling output:" << endl      
         << " -o file            Output vcf to this file (stdout)." << endl
         << " --Sample           use this sample name in the vcf (sample)" << endl
	 << " -C contig          Only run hmm on this chrom." << endl
         << " -B bed             Write coverage bed to this file." << endl
         << " -P model           Write trained parameter file." << endl
	 << " -M (flag)          Merge consecutive bins with the same copy number." << endl
	 << " -S snvs            Write SNVs to this file." << endl
	 << " -h help            Print this help message." << endl;

}
int main(int argc, const char* argv[]) {
  int nproc=4;
  double scale=2;
  
  if (argc < 3) {
    PrintHelp();
    exit(1);
  }
  NucMap.resize(256,4);
  NucMap[(int)'A']=0;
  NucMap[(int)'C']=1;
  NucMap[(int)'G']=2;
  NucMap[(int)'T']=3;	
  int maxState=10;	
  string bamFileName="";
  if (strcmp(argv[1], "-h") == 0) {
    PrintHelp();
    exit(1);
  }
  string referenceName=argv[1];

  MODEL_TYPE model=NEG_BINOM;
  string useChrom="";
  string hmmChrom="";
  string covBedInFileName="", covBedOutFileName="";
  bool   mergeBins=false;
  string outBedName="";
  string outFileName="";
  string snvFile="";
  string paramInFile="";
  string paramOutFile="";
  string hmmContig="";
  string snvInFileName="", snvOutFileName="";
  string clipInFileName="", clipOutFileName="";
  string sampleName="sample";
  int averageReadLength=0;
  if (argc > 2) {
    int argi=2;
    while (argi < argc) {
      if (strcmp(argv[argi], "-a") == 0) {
	++argi;
	bamFileName=argv[argi];
      }
      else if (strcmp(argv[argi], "-s") == 0) {
	++argi;
	snvInFileName=argv[argi];
      }
      else if (strcmp(argv[argi], "-S") == 0) {
	++argi;
	snvOutFileName=argv[argi];
      }      
      else if (strcmp(argv[argi], "-t") == 0) {
	++argi;
	nproc=atoi(argv[argi]);
      }
      else if (strcmp(argv[argi], "-p") == 0) {
	++argi;
	paramInFile=argv[argi];
      }
      else if (strcmp(argv[argi], "-P") == 0) {
	++argi;
	paramOutFile=argv[argi];
      }      
      else if (strcmp(argv[argi], "-b") == 0) {
	++argi;
	covBedInFileName = argv[argi];
      }
      else if (strcmp(argv[argi], "-B") == 0) {
	++argi;
	covBedOutFileName = argv[argi];
      }      
      else if (strcmp(argv[argi], "-e") == 0) {
	++argi;
	lepsi=atof(argv[argi]);
      }
      else if (strcmp(argv[argi], "--scale") == 0) {
	++argi;
	scale=atof(argv[argi]);
      }
      else if (strcmp(argv[argi], "-m") == 0) {
	++argi;
	if (strcmp(argv[argi], "pois") == 0) {
	  model=POIS;
	}
      }
      else if (strcmp(argv[argi], "-l") == 0) {
	++argi;
	clipInFileName=argv[argi];
      }
      else if (strcmp(argv[argi], "-L") == 0) {
	++argi;
	clipOutFileName=argv[argi];
      }      
      else if (strcmp(argv[argi], "-o") == 0) {
	++argi;
	outFileName = argv[argi];
      }
      /*      else if (strcmp(argv[argi], "-h") == 0) {
	++argi;
	hmmContig = argv[argi];
	}*/
      else if (strcmp(argv[argi], "-c") == 0) {
	++argi;
	useChrom = argv[argi];
      }
      else if (strcmp(argv[argi], "-C") == 0) {
	++argi;
	hmmChrom = argv[argi];
      }
      else  if (strcmp(argv[argi], "-M") == 0) {
	mergeBins=true;
      }
      else if (strcmp(argv[argi], "--earlyExit") == 0) {
	++argi;
      }
      else if (strcmp(argv[argi], "--sample") == 0) {
	++argi;
	sampleName=argv[argi];
      }
      else if (strcmp(argv[argi], "-h") == 0) {
	PrintHelp();
	exit(1);
      }
      
      else {	
	PrintHelp();
	cerr << "Invalid argument " << argv[argi] << endl;
	exit(1);
      }
      ++argi;
    }
  }
  //
  // Check some command line sanity.
  //

  if (covBedInFileName != "" and covBedOutFileName != "") {
    cerr << "ERROR. Cannot specify -b and -B." << endl;
    exit(1);
  }
  if (covBedInFileName == "" and bamFileName == "") {
    cerr << "ERROR. Must specify either a coverage file or a bam file" << endl;
    exit(1);
  }
  
  
  string faiFileName=referenceName + ".fai";
  vector<string> contigNames, allContigNames;
  vector<int>    contigLengths, allContigLengths;
  //
  // Determine what chroms to operate on.
  //
  ReadFai(faiFileName, allContigNames, allContigLengths);
  if (hmmChrom == "") {
    contigNames = allContigNames;
    contigLengths = allContigLengths;
  }
  else {
    contigNames.push_back(hmmChrom);
    for (int i=0; i < allContigNames.size(); i++) {
      if (allContigNames[i] == hmmChrom) {
	contigLengths.push_back(allContigLengths[i]);
	break;
      }
    }
    if (contigLengths.size() ==  0) {
      cerr << "ERROR. Could not find contig for hmm " << hmmContig << endl;
      exit(1);
    }
  }
  
  vector<vector<int> > covBins;
  vector<vector<int> > clipBins;
  double mean;
  double var;
  int nStates;
  int maxCov;
  vector<double> startP;
  vector<vector<double> > covCovTransP, covSnvTransP, snvSnvTransP;
  vector<vector<double> > updateTransP;
  vector<vector<SNV> > snvs;
  vector<vector<int > > copyNumber;
  vector< vector<double> > fCov, bCov, fSNV, bSNV;
  vector<vector<double> > emisP;
  vector<vector<double> > updateEmisP;
  vector<vector<vector< double> > > binoP;
  vector<vector<double> > expCovCovTransP, expCovSnvTransP, expSnvSnvTransP, expEmisP;
  vector<int> nReads;
  vector<long> totalBases;
  if (covBedInFileName != "") {
    ReadCoverage(covBedInFileName, contigNames, covBins);
  }
  if (clipInFileName != "") {
    ReadCoverage(clipInFileName, contigNames, clipBins);
  }

  if (snvInFileName != "") {
      ReadSNVs(snvInFileName, contigNames, snvs);
  }
  
  if (paramInFile != "") {
     ReadParameterFile(paramInFile, nStates,
		       mean, var, maxState, maxCov,
		       startP, covCovTransP, emisP);
  }
  
  //
  // Get the header of the bam file
  //
  htsFile *htsfp=NULL;
  bam_hdr_t *samHeader=NULL;
  hts_idx_t *bamidx=NULL;

  
  if (bamFileName != "") {
    htsfp = hts_open(bamFileName.c_str(),"r");
    samHeader = sam_hdr_read(htsfp);
    if ((bamidx = sam_index_load(htsfp, bamFileName.c_str())) == 0) {
      cerr << "ERROR reading index" << endl;
      exit(0);
    }
    const htsFormat *fmt = hts_get_format(htsfp);
    if (fmt == NULL or (fmt->format != sam and fmt->format != bam)) {
      cout << "Cannot determine format of input reads." << endl;
      exit(1);
    }    
  }


  
  // 
  // Read index for random io
  //

  nproc = min(nproc, (int) contigNames.size());

  pthread_t *threads = new pthread_t[nproc];
  vector<ThreadInfo> threadInfo(nproc);		
  pthread_mutex_t semaphore;		
  pthread_mutex_init(&semaphore, NULL);
  pthread_attr_t *threadAttr = new pthread_attr_t[nproc];
  int curSeq=0;
  vector<vector<Interval> > copyIntervals;
  copyIntervals.resize(contigNames.size());
  nReads.resize(contigNames.size());
  totalBases.resize(contigNames.size());
  double pModel=0;

  for (int procIndex = 0; procIndex < nproc ; procIndex++) {
    pthread_attr_init(&threadAttr[procIndex]);
    if (bamFileName != "") {
      threadInfo[procIndex].htsfp = hts_open(bamFileName.c_str(),"r");
      
    }
    else {
      threadInfo[procIndex].htsfp = NULL;
    }
    threadInfo[procIndex].bamidx = bamidx;
    threadInfo[procIndex].samHeader=samHeader;
    threadInfo[procIndex].fai = fai_load_format(referenceName.c_str(), FAI_FASTA);
    if (nproc > 1) 
      threadInfo[procIndex].exit = true;    
    else
      threadInfo[procIndex].exit = false;       
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

  if (paramInFile == "") {
    //max cov value observed or upper cov bound -> max nState---------------
    nStates= std::min( maxState , MAX_CN  ) + 1; //+1 zeroth state
    MAX_CN=nStates+1;
    startP.resize(nStates);
    for(size_t i=0;i<(nStates);i++){
      startP[i]=log(1./(nStates));
    }   
  }

  //
  // Allocate coverage bins if not reading
  //
  if (snvInFileName == "") {
    snvs.resize(contigLengths.size());
  }
  if (clipInFileName == "") {
    clipBins.resize(contigLengths.size());
    for (int c=0; c < contigLengths.size(); c++ ) {
      clipBins[c].resize(contigLengths[c]/BIN_LENGTH);
    }    
  }
  if (covBedInFileName == "") {
    covBins.resize(contigLengths.size());

    for (int c=0; c < contigLengths.size(); c++ ) {
      covBins[c].resize(contigLengths[c]/BIN_LENGTH);
      clipBins[c].resize(contigLengths[c]/BIN_LENGTH);
      copyNumber[c].resize(contigLengths[c]/BIN_LENGTH);       
    }

    //
    // Compute coverage from bam.
    //
    int parseChromNProc=min(4,nproc);
    if (nproc > 1) {
      for (int procIndex = 0; procIndex < parseChromNProc; procIndex++) {
	pthread_attr_init(&threadAttr[procIndex]);	
	pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*)) ParseChrom, &threadInfo[procIndex]);
      }
      
      for (int procIndex = 0; procIndex < parseChromNProc ; procIndex++) {
	pthread_join(threads[procIndex], NULL);
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
    cerr << "Length cutoff of average read length " << averageReadLength << endl;
    if (covBedOutFileName != "" ) {
      WriteCovBed(covBedOutFileName, contigNames, covBins);
    }
    if (clipOutFileName != "") {
      WriteCovBed(clipOutFileName, contigNames, clipBins);
    }
    if (snvOutFileName != "") {
      WriteSNVs(snvOutFileName, contigNames, snvs);
    }
  }

  EstimateCoverage(bamFileName, covBins, allContigNames, allContigLengths, useChrom, mean, var);

  
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
  for (int c=0 ;c < contigNames.size(); c++) {
    vector<bool> tooClose(snvs[c].size(), false);
    for ( int i=1; i < snvs[c].size(); i++ ){
      if (snvs[c][i].pos - snvs[c][i-1].pos < 50) {
	tooClose[i] = true;
	tooClose[i-1] = true;
      }
    }
    int p=0;
    for ( int i=0; i < snvs[c].size(); i++ ) {
      if (tooClose[i] == false) {
	snvs[c][p] = snvs[c][i];
	p++;
      }
    }
    snvs[c].resize(p);
  }

  cerr << "Computing copy-number up to " << MAX_CN << endl;
  //----------------------------------------------------------------------

  if (mean == 0) {
    std::cerr << "mean is zero, Exiting" << endl;
    return EXIT_FAILURE;
  }
    
  // trans prob, scale 3->2 by overlap of pdf.

  poisson distribution1(3*mean/2);
  double result3=pdf(distribution1, 3*mean/2);

  poisson distribution2(2*mean/2);
  double result2=pdf(distribution2, 3*mean/2);

  double epsi23 = result2/result3;

  //*300

  //no epsi
  double small=-30;
  //  double beta =  small + ( scale * log(epsi23))  ;
  double beta=small;
  //mean no. of bins for cn=3 call
  int nSNVStates=3;
  double unif=log(1.0/nStates);
  if (paramInFile == "") {
    InitParams(covCovTransP, covSnvTransP, snvSnvTransP,
	       nStates, nSNVStates, log(1-exp(small)), log(exp(small)/(nStates-1)),
	       emisP, model, maxCov, mean, var, binoP);
    printModel(covCovTransP, &cerr);
    //    printEmissions(emisP);    

    //
    // Baum-Welch training.
    //
    
    double prevPX=0;
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
      for (int procIndex = 0; procIndex < nproc; procIndex++) {  
	pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*)) ThreadedBWE, &threadInfo[procIndex]);
      }
      for (int procIndex = 0; procIndex < nproc ; procIndex++) {
	pthread_join(threads[procIndex], NULL);
      }
      px = pModel;

      if (prevPX != 0 and px - prevPX < 100 and i > 1) {
	cerr << "Ending iteration after " << i << " steps" << endl;       
	covCovTransP = prevTransP;
	emisP  = prevEmisP;
	break;
      }
      prevPX=px;

      vector<vector<double> > priorCovCov;
      priorCovCov.resize(covCovTransP.size());
      int nSites=covBins[0].size();
      BaumWelchM(startP, covCovTransP, emisP, binoP,
		 model,
		 stateTotCov, stateNCov,
		 expCovCovTransP, expEmisP,
		 priorCovCov,
		 updateTransP, updateEmisP);
      
      printModel(updateTransP, &cerr);
      covCovTransP=updateTransP;      
    }

    //
    // Eventually this needs to update for some multi-chrom code.
    //
   
    
    if (paramOutFile != "") {
      WriteParameterFile(paramOutFile, nStates, mean, var, maxState, maxCov, startP, covCovTransP, emisP);
    }
  }

  //
  // Now filter cn=1 and cn=3 intervals, or short calls
  //
  for (int c=0; c < contigNames.size(); c++) {
    int snvStart=0;
    for (int i=0; i < copyIntervals[c].size(); i++) {
      int curCN =copyIntervals[c][i].copyNumber;
      if ( curCN == 1 or curCN == 3) {
	while (snvStart < snvs[c].size() and snvs[c][snvStart].pos < copyIntervals[c][i].start) {
	  snvStart++;
	}
	int snvEnd=snvStart;
	while (snvEnd < snvs[c].size() and snvs[c][snvEnd].pos < copyIntervals[c][i].end) {
	  snvEnd++;
	}

	double pCN=0, pCN2=0;
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
  if (outFileName != "") {
    outFile.open(outFileName.c_str());
    outPtr = &outFile;
  }
  else {
    outPtr = &cout;
  }
  WriteVCF(*outPtr, referenceName, sampleName, contigNames, contigLengths, copyIntervals);
}
