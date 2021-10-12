#include "../include/hmcnc.h"

#include <cassert>
#include <cstdlib>

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <boost/algorithm/string.hpp>

//
// Not super
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
             std::vector<int> &contigLengths) {
  contigNames.clear();
  contigLengths.clear();

  std::string line;
  std::vector<std::string> fields;
  while (std::getline(faiIn, line)) {
    fields.clear();
    boost::split(fields, line, boost::is_any_of("\t "));
    if (fields.size() < 2) {
      std::cerr << "ERROR. Invalid FAI input: '" << line << "'\n";
      exit(EXIT_FAILURE);
    }

    contigNames.push_back(fields[0]);
    contigLengths.push_back(std::stoi(fields[1]));
  }
}

void ReadFai(const std::string faiFileName,
             std::vector<std::string> &contigNames,
             std::vector<int> &contigLengths) {
  std::ifstream faiIn{faiFileName.c_str()};
  if (faiIn.good() == false) {
    std::cerr << "ERROR. Reference is not indexed, or could not open .fai file" << '\n';
    exit(1);
  }
  ReadFai(faiIn, contigNames, contigLengths);
}

void ReadParameterFile(std::istream& inFile, int &nStates, double &covMean,
                       double &covVar, int &maxState, int &maxCov,
                       std::vector<double> &startP,
                       std::vector<std::vector<double>> &transP,
                       std::vector<std::vector<double>> &emisP) {
  std::string spacer;
  std::string section;
  int nr, nc;

  //
  // header
  //
  inFile >> spacer >> nStates;
  inFile >> spacer >> covMean;
  inFile >> spacer >> covVar;
  inFile >> spacer >> maxState;
  inFile >> spacer >> maxCov;

  //
  // startP
  //
  inFile >> section;
  if (section != "startP") {
    std::cerr << "ERROR. Parameter file: expected startP section, found '"
              << section << "' instead.\n";
    exit(EXIT_FAILURE);
  }

  double val;
  for (int i=0; i < nStates; i++) {
    inFile >> val;
    startP.push_back(val);
  }
  if (startP.size() != nStates) {
    std::cerr << "ERROR. Parameter file: unexpected number of startP values.";
    exit(EXIT_FAILURE);
  }

  //
  // transP
  //
  inFile >> section >> nr >> nc;
  if (section != "transP") {
    std::cerr << "ERROR. Parameter file: expected transP section, found '"
              << section << "' instead.\n";
    exit(EXIT_FAILURE);
  }
  transP.resize(nr);
  for (int i=0; i < nr; i++) {
    for (int j=0; j < nc; j++) {
      inFile >> val;
      transP[i].push_back(val);
    }
  }

  //
  // emisP
  //
  inFile >> section >> nr >> nc;
  if (section != "emisP") {
    std::cerr << "ERROR. Parameter file: expected emisP section, found '"
              << section << "' instead.\n";
    exit(EXIT_FAILURE);
  }
  emisP.resize(nr);
  for (int i=0; i < nr; i++) {
    for (int j=0; j < nc; j++) {
      inFile >> val;
      emisP[i].push_back(val);
    }
  }
}

void ReadParameterFile(const std::string &fileName, int &nStates, double &covMean,
                       double &covVar, int &maxState, int &maxCov,
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

void WriteParameterFile(std::ostream &outFile, int nStates, double covMean,
                        double covVar, int maxState, int maxCov,
                        const std::vector<double> &startP,
                        const std::vector<std::vector<double>> &transP,
                        const std::vector<std::vector<double>> &emisP) {
  //
  // header
  //
  outFile << "nStates\t" << nStates << '\n'
	        << "covMean\t" << covMean << '\n'
	        << "covVar\t" << covVar  << '\n'
	        << "maxState\t" << maxState << '\n'
	        << "maxCov\t" << maxCov << '\n';

  //
  // startP
  //
  outFile << "startP" << '\n';
  for (const auto s : startP) {
    outFile << s << '\n';
  }

  bool firstColumn = true;

  //
  // transP
  //
  assert(!transP.empty());
  outFile << "transP\t" << transP.size() << '\t' << transP[0].size() << '\n';
  for (const auto& row : transP) {
    firstColumn = true;
    for (const auto tp : row) {
      if (!firstColumn) {
        outFile << '\t';
      }
      outFile << tp;
      firstColumn = false;
    }
    outFile << '\n';
  }

  //
  // emisP
  //
  assert(!emisP.empty());
  outFile << "emisP\t" << emisP.size() << '\t' << emisP[0].size() << '\n';
  for (const auto& row : emisP) {
    firstColumn = true;
    for (const auto ep : row) {
      if (!firstColumn) {
        outFile << '\t';
      }
      outFile << ep;
      firstColumn = false;
    }
    outFile << '\n';
  }
}

void WriteParameterFile(const std::string &fileName, int nStates, double covMean,
                        double covVar, int maxState, int maxCov,
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
  std:: ofstream snvOut{snvFileName.c_str()};
  WriteSNVs(snvOut, contigNames, snvs);
}

std::string version="0.8";
std::string reference;


void WriteVCF(std::ostream &out,
	      const std::string &refName,
	      const std::string &sampleName,
	      const std::vector<std::string> &contigNames,
	      const std::vector<int> &contigLengths,
	      const std::vector<std::vector<Interval> > &intervals,
	      bool writeFail=false) {
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
      << "##INFO=<ID=REGION,Number=1,Type=String,Description=\"Region of interval "
    "for easy copy\">"
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
      << "##FORMAT=<ID=BN,Number=1,Type=Float,Description=\"Likelihood ratio of CN=2 vs "
    "CN=1 or CN=3 for heterozygous snvs\">"   
      << '\n'
      << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sampleName
      << '\n';
  for (size_t c = 0; c < contigNames.size(); c++) {
    assert(c < intervals.size());
    for (size_t i = 0; i < intervals[c].size(); i++) {
      if (intervals[c][i].copyNumber != 2) {

        const std::string cntype = (intervals[c][i].copyNumber > 2) ? "DUP" : "DEL";
	if (intervals[c][i].filter == "FAIL" and writeFail == false) {
	  continue;
	}
        out << contigNames[c] << '\t' << intervals[c][i].start
            << "\t.\t<CNV>\t<CNV>\t30\t" << intervals[c][i].filter << '\t'
            << "SVTYPE=" << cntype << ";"
            << "END=" << intervals[c][i].end
            << ";SVLEN=" << intervals[c][i].end - intervals[c][i].start
	    << ";REGION="<<contigNames[c] << ":" << intervals[c][i].start << "-" << intervals[c][i].end
            << ";IMPRECISE\t"
            << "CN:PP:DP"
	    << intervals[c][i].altInfo << "\t"
	    << intervals[c][i].copyNumber << ":"
            << intervals[c][i].pVal << ":" << intervals[c][i].averageCoverage
	    << intervals[c][i].altSample
            << '\n';
      }
    }
  }
}
