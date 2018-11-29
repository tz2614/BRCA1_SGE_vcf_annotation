#!usr/bin/python3

import BRCA1_SGE_ref
import csv
import pytest
import glob
import os
import shutil
import subprocess

# directories containing test datasets for the BRCA1_SGE_vcf_annotator program
working_files_dir = "/mnt/storage/home/zhengt/Competencies/CBI-6_Advanced_Clinical_Bioinformatics/BRCA1_SGE_vcf_annotator"

# other file types used to annotate vcf files with BRCA1_SGE annotations

ref_csv = "41586_2018_461_MOESM3_ESM.csv"
ref_vcf_hdr = "BRCA1_SGE_ref.vcf.hdr"

# declare current directory
current_directory = os.getcwd()

# new directories to be created for test modules
dir0 = os.path.abspath(os.path.join(current_directory, "pytest_dir0")) # test to see if you have write permission to current directory
test_csv_dir = os.path.abspath(os.path.join(current_directory, "pytest_BRCA1_csv")) # temp directory to store BRCA1_SGES_ref.csv/hdr files

def test_directory_does_not_exist(): 

	"""this checks whether the user has write permission to the current directory where the test directory will be stored"""
	assert os.makedirs(dir0), "The user do not have permission to create folders in current_directory"
	if os.path.exists(dir0):
		os.rmdir(dir0)

@pytest.mark.usefixtures("create_test_csv_dir")
def compare_vcf2csv(create_test_csv_dir):

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
				line_num = index

			elif index > line_num:
				fields = line.strip().split("\t")
				chrom = fields[0]
				pos = fields[1]
				ref = fields[3]
				alt = fields[4]
				score = fields[7]
				variant_dict[chrom][pos][ref] = [alt, score]

	with open (new_ref_csv, "r") as csv_file:

		lines = csv.reader(csv_file)

		for line in lines:

			if line[0] != "BRCA1":
				continue

			else:
				assert variant_dict[line[1]][line[2]][line[3]] = [line[16], line[17]], "variant not found in vcf"










