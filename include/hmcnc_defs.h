#ifndef HMCNC_DEFS_H
#define HMCNC_DEFS_H

#include <limits>
#include <string>
#include <vector>

#include "CLI11.hpp"

// ----------------------
// globals
// ----------------------

const int BIN_LENGTH=100;
const int MIN_CLIP_LENGTH=500;
const float MISMAP_RATE=0.01;
const double minNeg=-1*std::numeric_limits<double>::epsilon();
extern int MAX_CN;
extern double lepsi;

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

#endif // HMCNC_DEFS_H
