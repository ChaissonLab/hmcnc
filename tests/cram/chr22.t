  $ export DATA="${TESTDIR}"/../data

# --------------
# Setup
# --------------

Extract fields from input VCF
  $ grep -v "^#" "${DATA}"/chr22.vcf | cut -f 2 > "${CRAMTMP}"/expected_positions.txt
  $ grep -v "^#" "${DATA}"/chr22.vcf | cut -f 8 > "${CRAMTMP}"/expected_info.txt
  $ grep -v "^#" "${DATA}"/chr22.vcf | cut -f 10 | cut -f 1,3 -d ':' > "${CRAMTMP}"/expected_sample_info.txt

# --------------
# Normal run
# --------------

Run hmcnc on chr22 - generate VCF (and parameter file)
  $ "${__HMCNC_EXE}" "${DATA}"/chr22.fa -b "${DATA}"/chr22.bed -t 4 -P "${CRAMTMP}"/out.params -o "${CRAMTMP}"/out.vcf 2> /dev/null

Extract fields from result VCF
  $ grep -v "^#" "${CRAMTMP}"/out.vcf | cut -f 2 > "${CRAMTMP}"/observed_positions.txt
  $ grep -v "^#" "${CRAMTMP}"/out.vcf | cut -f 8 > "${CRAMTMP}"/observed_info.txt
  $ grep -v "^#" "${CRAMTMP}"/out.vcf | cut -f 10 | cut -f 1,3 -d ':' > "${CRAMTMP}"/observed_sample_info.txt

Compare
  $ diff "${CRAMTMP}"/expected_positions.txt "${CRAMTMP}"/observed_positions.txt
  $ diff "${CRAMTMP}"/expected_info.txt "${CRAMTMP}"/observed_info.txt
  $ diff "${CRAMTMP}"/expected_sample_info.txt "${CRAMTMP}"/observed_sample_info.txt

# ----------------------------
# Use input parameter file
# ----------------------------

Run hmcnc on chr22 - using existing parameter file
  $ "${__HMCNC_EXE}" "${DATA}"/chr22.fa -b "${DATA}"/chr22.bed -t 4 -p "${CRAMTMP}"/out.params -o "${CRAMTMP}"/from_params.vcf 2> /dev/null

Extract fields from result VCF
  $ grep -v "^#" "${CRAMTMP}"/from_params.vcf | cut -f 2 > "${CRAMTMP}"/from_params_positions.txt
  $ grep -v "^#" "${CRAMTMP}"/from_params.vcf | cut -f 8 > "${CRAMTMP}"/from_params_info.txt
  $ grep -v "^#" "${CRAMTMP}"/from_params.vcf | cut -f 10 | cut -f 1,3 -d ':' > "${CRAMTMP}"/from_params_sample_info.txt

Compare
  $ diff "${CRAMTMP}"/expected_positions.txt "${CRAMTMP}"/from_params_positions.txt
  $ diff "${CRAMTMP}"/expected_info.txt "${CRAMTMP}"/from_params_info.txt
  $ diff "${CRAMTMP}"/expected_sample_info.txt "${CRAMTMP}"/from_params_sample_info.txt
