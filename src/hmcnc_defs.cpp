#include "../include/hmcnc_defs.h"

#include <cassert>
#include <cstdlib>

#include <iostream>

int MAX_CN=6;
double lepsi=-800;

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
{ }

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
    "    ** Need description for clipInFileName ** ")->
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
    "Write trained parameter file.")->
    group(outputGroupName)->
    type_name("FILE");

  CLI.add_option("-L", clipOutFileName,
    "    ** Need description for clipOutFileName **  ")->
    group(outputGroupName)->
    type_name("FILE");

  CLI.add_option("-S", snvOutFileName,
    "Write SNVs to this file.")->
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
      std::cerr << "ERROR. Cannot specify -b and -B.\n";
      exit(EXIT_FAILURE);
    }
    if (this->covBedInFileName == "" and this->bamFileName == "") {
      std::cerr << "ERROR. Must specify either a coverage file or a bam file\n";
      exit(EXIT_FAILURE);
    }
  });
}
