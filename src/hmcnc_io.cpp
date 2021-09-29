#include "../include/hmcnc_io.h"

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include <boost/algorithm/string.hpp>

void ReadCoverage(std::istream &covFile,
                  const std::vector<std::string> &contigNames,
                  std::vector<std::vector<int>> &covBins) {
  std::map<std::string, std::vector<int>> tempCovBins;

  int length = 0;
  std::string line;
  std::vector<std::string> fields;
  while (std::getline(covFile, line)) {
    fields.clear();
    boost::split(fields, line, boost::is_any_of("\t "));
    if (fields.size() < 4) {
      std::cerr << "ERROR. Invalid BED input: '" << line << "'\n";
      exit(EXIT_FAILURE);
    }

    const std::string &contig = fields.at(0);
    const int cov = std::stoi(fields.at(3));
    tempCovBins[contig].push_back(cov);
    length += line.size();
  }

  std::cerr << "read cov buffer of len " << length << '\n';

  covBins.clear();
  for (const auto &contig : contigNames) {
    const auto found = tempCovBins.find(contig);
    if (found == tempCovBins.cend()) {
      covBins.push_back(std::vector<int>{});
    } else {
      covBins.push_back(std::move(found->second));
    }
  }
}

void ReadCoverage(const std::string &covFileName,
                  const std::vector<std::string> &contigNames,
                  std::vector<std::vector<int>> &covBins) {
  std::ifstream covFile{covFileName.c_str()};
  ReadCoverage(covFile, contigNames, covBins);
}


void ReadFai(std::istream &faiIn,
             std::vector<std::string> &contigNames,
             std::vector<int> &contigLengths)
{
  if (faiIn.good() == false) {
    std::cerr << "ERROR. Reference is not indexed, or could not open .fai file" << '\n';
    exit(EXIT_FAILURE);
  }

  while (faiIn) {
    std::string line;
    getline(faiIn, line);
    std::stringstream strm(line);
    if (line != "") {
      std::string contig;
      int length;

      strm >> contig;
      strm >> length;
      contigNames.push_back(contig);
      contigLengths.push_back(length);
    }
  }
}

void ReadFai(const std::string faiFileName,
             std::vector<std::string> &contigNames,
             std::vector<int> &contigLengths)   {
  std::ifstream faiIn{faiFileName.c_str()};
  ReadFai(faiIn, contigNames, contigLengths);
}

void ReadParameterFile(std::istream &inFile,
                       int &nStates,
                       double &covMean,
                       double &covVar,
                       int &maxState,
                       int &maxCov,
                       std::vector<double> &startP,
                       std::vector<std::vector<double>> &transP,
                       std::vector<std::vector<double>> &emisP) {
  std::string spacer;
  inFile >> spacer >> nStates;
  inFile >> spacer >> covMean;
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

void ReadParameterFile(const std::string &fileName,
                       int &nStates,
                       double &covMean,
                       double &covVar,
                       int &maxState,
                       int &maxCov,
                       std::vector<double> &startP,
                       std::vector<std::vector<double>> &transP,
                       std::vector<std::vector<double>> &emisP) {
  std::ifstream inFile{fileName.c_str()};
  ReadParameterFile(inFile, nStates, covMean, covVar, maxState, maxCov,
                    startP, transP, emisP);
}

void ReadSNVs(std::istream &snvIn,
              const std::vector<std::string> &contigNames,
              std::vector<std::vector<SNV>> &snvs) {
  snvs.resize(contigNames.size());
  size_t curContig=0;

  std::string line;
  std::string chrom;
  int pos, ref, alt;
  char refNuc, altNuc, t;

  while (curContig < contigNames.size()) {
    snvIn >> chrom >> pos >> refNuc >> altNuc >> ref >> alt;
    if (chrom == "" or snvIn.eof()) {
      break;
    }
    while (curContig < contigNames.size() and chrom != contigNames[curContig]) {
      curContig++;
    }
    if (curContig < contigNames.size()) {
      snvs[curContig].push_back(SNV{pos, refNuc, altNuc, ref, alt});
    }
  }
}

void ReadSNVs(const std::string &snvFileName,
              const std::vector<std::string> &contigNames,
              std::vector<std::vector<SNV>> &snvs) {
  std::ifstream snvIn{snvFileName};
  ReadSNVs(snvIn, contigNames, snvs);
}

void WriteCovBed(std::ostream &covFile,
		             const std::vector<std::string> &contigNames,
		             const std::vector<std::vector<int>> &covBins) {
  for (size_t c=0; c < contigNames.size(); c++) {
    assert(c < covBins.size());
    const auto &contigName = contigNames[c];
    const auto &contigBins = covBins[c];
    for (size_t i=0; i < contigBins.size(); i++) {
      covFile << contigName << '\t'
              << i*100 << '\t'
              << (i+1)*100 << '\t'
              << contigBins[i] << '\n';
    }
  }
}

void WriteCovBed(const std::string &covFileName,
		             const std::vector<std::string> &contigNames,
		             const std::vector<std::vector<int>> &covBins) {
  std::ofstream covFile{covFileName.c_str()};
  WriteCovBed(covFile, contigNames, covBins);
}

void WriteParameterFile(std::ostream &outFile,
                        int nStates,
                        double covMean,
                        double covVar,
                        int maxState,
                        int maxCov,
                        const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &transP,
                        const std::vector<std::vector<double>> &emisP) {
  outFile << "nStates\t" << nStates << '\n'
	        << "covMean\t" << covMean << '\n'
	        << "covVar\t" << covVar  << '\n'
	        << "maxState\t" << maxState << '\n'
	        << "maxCov\t" << maxCov << '\n'
	        << "startP" << '\n';
  for (size_t i=0; i < startP.size(); i++) {
    outFile << startP[i];
    if (i+1 < startP.size()) {
      outFile << '\t';
    }
    outFile << '\n';
  }

  assert(!transP.empty());
  outFile << "transP\t" << transP.size() << '\t' << transP[0].size() << '\n';
  for (size_t i=0; i < transP.size(); i++) {
    for (size_t j=0; j < transP[j].size(); j++) {
      outFile << transP[i][j];
      if (i+1 < transP.size()) {
        outFile << '\t';
      }
    }
    outFile << '\n';
  }

  assert(!emisP.empty());
  outFile << "emisP\t" << emisP.size() << '\t' << emisP[0].size() << '\n';
  for (size_t i=0; i < emisP.size(); i++) {
    for (size_t j=0; j < emisP[j].size(); j++) {
      outFile << emisP[i][j];
      if (i+1 < emisP.size()) {
        outFile << '\t';
      }
    }
    outFile << '\n';
  }
}

void WriteParameterFile(const std::string &fileName,
                        int nStates,
                        double covMean,
                        double covVar,
                        int maxState,
                        int maxCov,
                        const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &transP,
                        const std::vector<std::vector<double>> &emisP) {
  std::ofstream outFile{fileName.c_str()};
  WriteParameterFile(outFile, nStates, covMean, covVar, maxState, maxCov,
                     startP, transP, emisP);
}

void WriteSNVs(std::ostream &snvOut,
               const std::vector<std::string> &contigNames,
               const std::vector<std::vector<SNV>> &snvs) {
  for (size_t c=0; c < contigNames.size(); c++) {
    assert(c < snvs.size());
    for (size_t i=0; i < snvs[c].size(); i++) {
      snvOut << contigNames[c] << '\t'
             << snvs[c][i].pos << '\t'
             << snvs[c][i].refNuc << '\t'
             << snvs[c][i].altNuc << '\t'
             << snvs[c][i].ref << '\t'
             << snvs[c][i].alt << '\n';
    }
  }
}

void WriteSNVs(const std::string &snvFileName,
               const std::vector<std::string> &contigNames,
               const std::vector<std::vector<SNV>> &snvs) {
  std::ofstream snvOut{snvFileName.c_str()};
  WriteSNVs(snvOut, contigNames, snvs);
}




void WriteVCF(std::ostream &out,
	      const std::string &refName,
	      const std::string &sampleName,
	      const std::vector<std::string> &contigNames,
	      const std::vector<int> &contigLengths,
	      const std::vector<std::vector<Interval>> &intervals) {

  const std::string version{"0.8"};
  const std::string reference;

  out << "##fileformat=VCFv4.1" << '\n'
      << "##source=hmmcnc_v" << version << '\n'
      << "##reference=" << reference << '\n';
  for (size_t i = 0; i < contigNames.size(); i++) {
    out << "##contig=<ID=" << contigNames[i] << ",length=" << contigLengths[i]
        << ">" << '\n';
  }

  out << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of "
    "structural variant\">"
      << '\n'
      << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of "
    "the structural variant described in this record\">"
      << '\n'
      << "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in "
    "length between REF and ALT alleles\">"
      << '\n'
      << "##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise "
    "structural variation\">"
      << '\n';
  out << "##FORMAT=<ID=CN,Number=1,Type=String,Description=\"CopyNumber\">"
      << '\n'
      << "##FORMAT=<ID=PP,Number=R,Type=Float,Description=\"Relative posterior "
    "probability (phred)\">"
      << '\n'
      << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth at "
    "this position for this sample\">"
      << '\n'
      << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName
      << '\n';
  for (size_t c = 0; c < contigNames.size(); c++) {
    assert(c < intervals.size());
    for (size_t i = 0; i < intervals[c].size(); i++) {
      if (intervals[c][i].copyNumber != 2) {

        const std::string cntype = (intervals[c][i].copyNumber > 2) ? "DUP" : "DEL";

        out << contigNames[c] << '\t' << intervals[c][i].start
            << "\t.\t<CNV>\t<CNV>\t30\t" << intervals[c][i].filter << '\t'
            << "SVTYPE=" << cntype << ";"
            << "END=" << intervals[c][i].end
            << ";SVLEN=" << intervals[c][i].end - intervals[c][i].start
            << ";IMPRECISE\t"
            << "CN:PP:DP\t" << intervals[c][i].copyNumber << ":"
            << intervals[c][i].pVal << ":" << intervals[c][i].averageCoverage
            << '\n';
      }
    }
  }
}
