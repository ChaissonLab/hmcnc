  $ export DATA="${TESTDIR}"/../data

Extract fields from input VCF
  $ grep -v "^#" "${DATA}"/chr22.vcf | cut -f 2 > expected_positions.txt
  $ grep -v "^#" "${DATA}"/chr22.vcf | cut -f 8 > expected_info.txt
  $ grep -v "^#" "${DATA}"/chr22.vcf | cut -f 10 | cut -f 1-3 -d ':' > expected_sample_info.txt

Run hmcnc
  $ "${__HMCNC_EXE}" "${DATA}"/chr22.fa -b "${DATA}"/chr22.bed -t 4 -o out.vcf 2> /dev/null

Extract fields from result VCF
  $ grep -v "^#" out.vcf | cut -f 2 > observed_positions.txt
  $ grep -v "^#" out.vcf | cut -f 8 > observed_info.txt
  $ grep -v "^#" out.vcf | cut -f 10 | cut -f 1-3 -d ':' > observed_sample_info.txt

Compare
  $ diff expected_positions.txt observed_positions.txt
  $ diff expected_info.txt observed_info.txt
  $ diff expected_sample_info.txt observed_sample_info.txt
