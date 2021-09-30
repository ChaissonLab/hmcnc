#include <gtest/gtest.h>

#include "../include/hmcnc_io.h"

#include <algorithm>
#include <sstream>
#include <string>

std::string Normalized(std::string s) {
  std::replace(s.begin(), s.end(), ' ', '\t');
  return s;
}

std::string CoverageBedText()
{
    const std::string bedText{R"(
chr1 0 100 23
chr1 100 200 23
chr1 200 300 25
chr1 300 400 26
chr1 400 500 28
chr1 500 600 29
chr1 600 700 30
chr1 700 800 31
chr1 800 900 31
chr1 900 1000 31
chr1 1000 1100 31
chr1 1100 1200 31
chr1 1200 1300 31
chr1 1300 1400 31
chr1 1400 1500 31
chr1 1500 1600 32
chr2 0 100 17
chr2 100 200 18
chr2 200 300 18
chr2 300 400 20
)"};
    return bedText.substr(1);
}

TEST(hmcnc_io, normal_coverage_bed_can_do_roundtrip) {

  const std::string inputText = Normalized(CoverageBedText());
  std::istringstream bedIn{inputText};

  const std::vector<std::string> contigNames{"chr1", "chr2"};
  std::vector<std::vector<int>> covBins;
  ReadCoverage(bedIn, contigNames, covBins);

  std::ostringstream bedOut;
  WriteCovBed(bedOut, contigNames, covBins);

  EXPECT_EQ(Normalized(bedOut.str()), inputText);
}

std::string FaiText() {
      const std::string faiText{R"(
chr1 1600 0 70 71
chr2 400 1601 70 71
)"};
    return faiText.substr(1);
}

TEST(hmcnc_io, can_read_normal_fai) {

  std::istringstream faiIn{Normalized(FaiText())};
  std::vector<std::string> contigNames;
  std::vector<int> contigLengths;

  ReadFai(faiIn, contigNames, contigLengths);

  ASSERT_EQ(contigNames.size(), 2);
  EXPECT_EQ(contigNames[0], "chr1");
  EXPECT_EQ(contigNames[1], "chr2");

  ASSERT_EQ(contigLengths.size(), 2);
  EXPECT_EQ(contigLengths[0], 1600);
  EXPECT_EQ(contigLengths[1], 400);
}

std::string ParameterFileText() {
  return {
    "nStates\t7\n"
    "covMean\t38\n"
    "covVar\t114\n"
    "maxState\t10\n"
    "maxCov\t209\n"
    "startP\n"
    "-1.94591\n"
    "-1.94591\n"
    "-1.94591\n"
    "-1.94591\n"
    "-1.94591\n"
    "-1.94591\n"
    "-1.94591\n"
    "transP\t7\t7\n"
    "-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\n"
    "emisP\t7\t7\n"
    "-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\t-31.7918\n"
    "-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-31.7918\t-9.35918e-14\n"
  };
}

TEST(hmcnc_io, parameter_file_can_do_roundtrip) {

  const std::string inputText = Normalized(ParameterFileText());
  std::istringstream paramIn{inputText};

  // data
  int nStates = 42;
  double covMean = 42.;
  double covVar = 42.;
  int maxState = 42;
  int maxCov = 42;
  std::vector<double> startP;
  std::vector<std::vector<double>> transP;
  std::vector<std::vector<double>> emisP;
  ReadParameterFile(paramIn, nStates, covMean, covVar, maxState, maxCov,
                    startP, transP, emisP);

  std::ostringstream paramOut;
  WriteParameterFile(paramOut, nStates, covMean, covVar, maxState, maxCov,
                    startP, transP, emisP);

  EXPECT_EQ(Normalized(paramOut.str()), inputText);
}

