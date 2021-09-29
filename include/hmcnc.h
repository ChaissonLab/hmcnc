#ifndef HMCNC_H
#define HMCNC_H

#include <iostream>

#include "hmcnc_defs.h"

// ------------------------
// algorithm methods
// ------------------------

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

double max_over_row(const std::vector<std::vector<double>> &v,
                    size_t col, size_t nStates);
double max_over_rows(const std::vector<std::vector<double>> &v, size_t col,
                     const std::vector<std::vector<double>> &v2, size_t nextState,
                     size_t nStates);

// --------------------------
// main application runners
// --------------------------

// inject parameters for testing
int hmcnc(Parameters& params);

// initialize parameters from command line
int hmcnc(int argc, const char* argv[]);

int hmcnc_test();

#endif // HMCNC_H
