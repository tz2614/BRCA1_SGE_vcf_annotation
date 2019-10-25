#!usr/bin/python3

import BRCA1_SGE_ref
import csv
import pytest
import os
import shutil

def test_vcf2csv(create_test_csv_dir):

	""" make the BRCA1_ref_vcf using the create_ref_vcf function from BRCA1_SGE_vcf_annotator, 
	and check the variants match between vcf and csv files"""

	test_csv_dir, new_ref_csv = create_test_csv_dir
	new_ref_vcf = BRCA1_SGE_ref.create_ref_vcf(new_ref_csv)
	variant_dict = {}

	with open (new_ref_vcf, "r") as vcf_file:
	
		for index, line in enumerate(vcf_file):
			if line.startswith("##"):
				continue

			elif line.startswith("#CHROM"):
				# for each variant in vcf add the variant into the variant_dict as {chrom:{pos{ref: [alt, score]}}}
				line_num = index

				if index > line_num:
					fields = line.strip().split("\t")
					chrom = fields[0]
					pos = fields[1]
					ref = fields[3]
					alt = fields[4]
					score = fields[7]
					variant_dict[chrom][pos][ref] = [alt, score]

	with open (new_ref_csv, "r") as csv_file:

		lines = csv.reader(csv_file)

		for index, line in enumerate(lines):
			#for each variant in csv, check that the same variant is present in vcf

			if line[0] != "BRCA1" and index > 2:		
				
				assert variant_dict[line[1]][line[2]][line[3]] == [line[16], line[17]], "variant not found in vcf"
				continue

