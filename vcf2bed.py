#!usr/bin/python3

import sys
import os
import subprocess
import datetime
import pprint as pp


def vcf_file(vcf_file, cp_ref_text):

	cp_variants = []
	vcf_variants = []

	with open (cp_ref_text) as cp_ref:

		for cp_line in cp_ref:
			if cp_line.startswith("C"):
				cp_fields = cp_line.strip().split("\t")
				cp_fields = [cp_fields[1], cp_fields[2], cp_fields[3]]
				print (cp_fields)
				cp_variants.append(cp_fields)
			else:
				continue

	with open (vcf_file) as vcf_lines: 

		for vcf_line in vcf_lines:
			if vcf_line.startswith("17"):
				vcf_fields = vcf_line.strip().split("\t")
				vcf_fields = [vcf_fields[1], vcf_fields[3], vcf_fields[4]]
				vcf_variants.append(vcf_fields)
			else:
				continue
				
	return cp_variants, vcf_variants

def main(vcf_file, cp_ref_text):

	check_BRCA1_SGE_variants(vcf_file, cp_ref_text)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])