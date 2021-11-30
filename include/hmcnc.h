#ifndef HMCNC_H
#define HMCNC_H

#include <string>
#include <vector>

#include "CLI11.hpp"

// ----------------------
// data structs/enums
// ----------------------

typedef enum { POIS, NEG_BINOM  } MODEL_TYPE;

struct Interval {
  int start;
  int end;
  int copyNumber;
  float averageCoverage;
  double pVal;
  std::string filter;
  std::string altInfo;
  std::string altSample;
  int distanceToFrontClip;
  int distanceToEndClip;
  int nFrontClip;
  int nEndClip;
  Interval();
  Interval(int s, int e, int cn, float avg, double p);
};

struct SNV {
  int pos;
  char refNuc = '\0';
  char altNuc = '\0';
  int ref = 0;
  int alt = 0;

  SNV();
  SNV(int p, int r, int a, int rc, int ac);
  SNV(int p);

  bool operator<(const SNV &rhs) const;
};

struct Parameters {

  // positional arg
  std::string referenceName;

  // options
  std::string bamFileName;
  std::string snvInFileName;
  std::string snvOutFileName;
  std::string paramInFile;
  std::string paramOutFile;
  std::string covBedInFileName;
  std::string covBedOutFileName;
  std::string outBedName;
  std::string clipInFileName;
  std::string clipOutFileName;
  std::string outFileName;
  std::string useChrom;
  std::string hmmChrom;

  int nproc = 4;
  MODEL_TYPE model = NEG_BINOM;
  bool mergeBins=false;
  std::string sampleName;

  CLI::App CLI;
  std::string modelString;

  Parameters();
};

// ----------------------
// Algorithm methods
// ----------------------

double BaumWelchEOnChrom(const std::vector<double> &startP,
                         std::vector<std::vector<double>> &covCovTransP,
                         std::vector<std::vector<double>> &emisP,
                         std::vector<int> &obs,
                         std::vector<std::vector<double>> &f,
                         std::vector<std::vector<double>> &b,
                         std::vector<std::vector<double>> &expCovCovTransP,
                         std::vector<std::vector<double>> &expEmisP);

void BaumWelchM(const std::vector<double> &startP,
                const std::vector<std::vector<double>> &transP,
                const std::vector<std::vector<double>> &emisP,
                const std::vector<std::vector<std::vector<double>>> &binoP,
                int model,
                const std::vector<long> &stateTotCov,
                const std::vector<long> &stateNCov,
                const std::vector<std::vector<double>> &expTransP,
                std::vector<std::vector<double>> &expEmisP,
                std::vector<std::vector<double>> &covCovPrior,
                std::vector<std::vector<double>> &updateTransP,
                std::vector<std::vector<double>> &updateEmisP);

void CombineEmissions(const std::vector<int> &obs,
                      const std::vector<SNV> &snvs,
                      std::vector<uint8_t> &isCov,
                      std::vector<int> &obsIndex);

double CSEmisP(int state, int pos,
               const std::vector<int> &obs,
               const std::vector<SNV> &snvs,
               const std::vector<uint8_t> &isCov,
               const std::vector<int> &obsIndex,
               const std::vector<std::vector<double>> &emisP,
               const std::vector<std::vector<std::vector<double>>> &binoP);

double ForwardBackwards(const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &covCovTransP,
                        const std::vector<std::vector<double>> &emisP,
                        const std::vector<int> &obs,
                        std::vector<std::vector<double>> &f,
                        std::vector<std::vector<double>> &b);

int GetRefAndAlt(char refNuc, const std::vector<int> &counts,
                 int &ref, int &alt);

double LgBinom(double p, int s, int n);

double LgNegBinom(int cn, int cov, float Hmean, float Hvar);

double LgPrpoiss(int cn,  int cov, int Hmean);

void Moments(const std::vector<double> &v, double &ex, double &var);

double PairSumOfLogP(double a, double b);

bool compareInterval(Interval i1, Interval i2);

void mergeIntervals(std::vector<Interval> &intervals, std::vector<Interval> &mergedIntervals, std::string contig);

void intersectDelCall( std::vector<Interval> &mergedIntervals, std::vector<Interval> & copyIntervals, double mean);

void mergeNaiveIntervals(std::vector<Interval> &intervals, std::vector<Interval> &mergedIntervals, std::string contig);

int quant( std::vector<std::vector<Interval>> &mergedIntervals ,double q, std::vector<std::string> &contigNames );

void NaiveCaller(std::vector<int> &covBins, std::vector<Interval> &NaiveIntervals, double mean );

void StorePosteriorMaxIntervals(const std::vector<int> &cov,
				                        const std::vector<std::vector<double>> &f,
				                        const std::vector<std::vector<double>> &b,
				                        std::vector<Interval> &intervals);

int StoreSNVs(char *contigSeq, int contigLength, float mean,
              std::vector<int> &nA, std::vector<int> &nC,
              std::vector<int> &nG, std::vector<int> &nT,
              std::vector<int> &nDel, std::vector<SNV> &snvs);

double SumOfLogP(const std::vector<double> &vals);

void UpdateEmisP(std::vector<std::vector<double>> &emisP,
                 std::vector<std::vector<double>> &expEmisP,
                 int model=NEG_BINOM);



// --------------------------
// I/O methods
// --------------------------

void ReadCoverage(std::istream &covFile,
                  const std::vector<std::string> &contigNames,
                  std::vector<std::vector<int>> &covBins);
void ReadCoverage(const std::string &covFileName,
                  const std::vector<std::string> &contigNames,
                  std::vector<std::vector<int>> &covBins);

void ReadFai(std::istream &faiFile,
             std::vector<std::string> &contigNames,
             std::vector<int> &contigLengths);
void ReadFai(const std::string faiFileName,
             std::vector<std::string> &contigNames,
             std::vector<int> &contigLengths);

void ReadParameterFile(std::istream &file, int &nStates, double &covMean,
                       double &covVar, int &maxState, int &maxCov,
                       std::vector<double> &startP,
                       std::vector<std::vector<double>> &transP,
                       std::vector<std::vector<double>> &emisP);
void ReadParameterFile(const std::string &fileName, int &nStates, double &covMean,
                       double &covVar, int &maxState, int &maxCov,
                       std::vector<double> &startP,
                       std::vector<std::vector<double>> &transP,
                       std::vector<std::vector<double>> &emisP);

void ReadSNVs(std::istream &snvFile,
              const std::vector<std::string> &contigNames,
              std::vector<std::vector<SNV>> &snvs);
void ReadSNVs(const std::string &snvFileName,
              const std::vector<std::string> &contigNames,
              std::vector<std::vector<SNV>> &snvs);

void WriteCovBed(std::ostream &covFile,
		             const std::vector<std::string> &contigNames,
		             const std::vector<std::vector<int>> &covBins);
void WriteCovBed(const std::string &covFileName,
		             const std::vector<std::string> &contigNames,
		             const std::vector<std::vector<int>> &covBins);

void WriteParameterFile(std::ostream &file, int nStates, double covMean,
                        double covVar, int maxState, int maxCov,
                        const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &transP,
                        const std::vector<std::vector<double>> &emisP);
void WriteParameterFile(const std::string &fileName, int nStates, double covMean,
                        double covVar, int maxState, int maxCov,
                        const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &transP,
                        const std::vector<std::vector<double>> &emisP);

void WriteSNVs(std::ostream &snvFile,
               const std::vector<std::string> &contigNames,
               const std::vector<std::vector<SNV>> &snvs);

void WriteSNVs(const std::string &snvFileName,
               const std::vector<std::string> &contigNames,
               const std::vector<std::vector<SNV>> &snvs);

void WriteBed(  const std::vector <std::vector<Interval>> &intv, 
                std::ostream &out, 
                const std::vector<std::string> &contigNames);

void WriteBed(const std::vector<std::vector<Interval>> &intv,
            const std::string &bedFileName,   
            const std::vector<std::string> &contigNames);



void WriteVCF(std::ostream &out, const std::string &refName,
              const std::string &sampleName,
              const std::vector<std::string> &contigNames,
              const std::vector<int> &contigLengths,
              const std::vector<std::vector<Interval>> &intervals,
	      bool writeFail);

// --------------------------
// main application runners
// --------------------------

// inject parameters for testing
int hmcnc(Parameters& params);

// initialize parameters from command line
int hmcnc(int argc, const char* argv[]);

int hmcnc_test();

#endif // HMCNC_H
