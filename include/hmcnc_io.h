#ifndef HMCNC_IO_H
#define HMCNC_IO_H

#include <iosfwd>
#include <string>
#include <vector>

#include "hmcnc_defs.h"

// --------------------------
// serialization
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

void ReadParameterFile(std::istream &file,
                       int &nStates,
                       double &covMean,
                       double &covVar,
                       int &maxState,
                       int &maxCov,
                       std::vector<double> &startP,
                       std::vector<std::vector<double>> &transP,
                       std::vector<std::vector<double>> &emisP);
void ReadParameterFile(const std::string &fileName,
                       int &nStates,
                       double &covMean,
                       double &covVar,
                       int &maxState,
                       int &maxCov,
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

void WriteParameterFile(std::ostream &file,
                        int nStates,
                        double covMean,
                        double covVar,
                        int maxState,
                        int maxCov,
                        const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &transP,
                        const std::vector<std::vector<double>> &emisP);
void WriteParameterFile(const std::string &fileName,
                        int nStates,
                        double covMean,
                        double covVar,
                        int maxState,
                        int maxCov,
                        const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &transP,
                        const std::vector<std::vector<double>> &emisP);

void WriteSNVs(std::ostream &snvFile,
               const std::vector<std::string> &contigNames,
               const std::vector<std::vector<SNV>> &snvs);
void WriteSNVs(const std::string &snvFileName,
               const std::vector<std::string> &contigNames,
               const std::vector<std::vector<SNV>> &snvs);

void WriteVCF(std::ostream &out,
	      const std::string &refName,
	      const std::string &sampleName,
	      const std::vector<std::string> &contigNames,
	      const std::vector<int> &contigLengths,
	      const std::vector<std::vector<Interval>> &intervals);

#endif // HMCNC_IO_H
