#!usr/bin/python3

import BRCA1_SGE_vcf_annotator
import pytest
import glob
import os
import shutil
import subprocess

# directories containing test datasets for the BRCA1_SGE_vcf_annotator program
BRCA1_rf_path = "/mnt/storage/home/zhengt/Competencies/CBI-6_Advanced_Clinical_Bioinformatics/BRCA1_test/"
working_files_dir = "/mnt/storage/home/zhengt/Competencies/CBI-6_Advanced_Clinical_Bioinformatics/BRCA1_SGE_vcf_annotator"

# vcf filenames used within the test directories
pos_vcf_filename = "pos_variants.vcf"
neg_vcf_filename = "neg_variants.vcf"

# annotated vcf filenames used within the test directories
pos_annotated_vcf_filename = "pos_variants.annotated.vcf"
neg_annotated_vcf_filename = "neg_variants.annotated.vcf"

# other file types used to annotate vcf files with BRCA1_SGE annotations

ref_csv = "41586_2018_461_MOESM3_ESM.csv"
ref_vcf = "BRCA1_SGE_ref.vcf"
ref_vcf_gz = "BRCA1_SGE_ref.vcf.gz"
ref_vcf_gz_tbi = "BRCA1_SGE_ref.vcf.gz.tbi"
ref_vcf_hdr = "BRCA1_SGE_ref.vcf.hdr"

# declare current directory
current_directory = os.getcwd()

# new directories to be created for test modules
directory0 = os.path.abspath(os.path.join(current_directory, "pytest_dir0")) # test to see if you have write permission to current directory
test_vcf_dir = os.path.abspath(os.path.join(current_directory, "pytest_BRCA1_vcf")) # temp directory to store all test vcf and BRCA1_annotated.vcf
test_work_dir = os.path.abspath(os.path.join(current_directory, "pytest_work")) # temp directory to store BRCA1_SGES_ref.csv/hdr files

@pytest.fixture(scope="function")
def create_test_vcfs_dir():

	"""
    Set up test directory
    """

	# if the directory already exist then error
	assert not os.path.exists(test_vcf_dir), "The test directory {} already exists".format(test_vcf_dir)

	# otherwise make it in preparation for testing
	os.makedirs(test_vcf_dir)
	
	org_pos_vcf = os.path.join(BRCA1_rf_path, pos_vcf_filename)
	org_neg_vcf = os.path.join(BRCA1_rf_path, neg_vcf_filename)

	new_pos_vcf = os.path.join(test_vcf_dir, pos_vcf_filename)
	new_neg_vcf = os.path.join(test_vcf_dir, neg_vcf_filename)

	shutil.copyfile(pos_vcf, new_pos_vcf)
	shutil.copyfile(neg_vcf, new_neg_vcf)
	
	yield test_vcf_dir, new_pos_vcf, new_neg_vcf  # This gets passed to the test function

	# This runs after the test function is complete to clean up
	# It deletes the test directory we created and anything it contains
	shutil.rmtree(test_vcf_dir)


@pytest.fixture(scope="function")
def create_test_csv_dir():

	"""
    Set up test directory
    """

	# if the directory already exist then error
	assert not os.path.exists(test_work_dir), "The test directory {} already exists".format(test_work_dir)

	# otherwise make it in preparation for testing
	os.makedirs(test_work_dir)
	
	org_ref_csv = os.path.join(test_work_dir, ref_csv)

	shutil.copyfile(org_ref_csv, new_ref_csv)
	
	yield test_csv_dir, new_ref_csv # This gets passed to the test function

	# This runs after the test function is complete to clean up
	# It deletes the test directory we created and anything it contains
	shutil.rmtree(test_work_dir)


@pytest.fixture(scope="function")
def create_work_dir():

	"""
	Set up working directory to store 
	BRCA1_SGE_ref.csv, BRCA1_SGE_ref.vcf, BRCA1_SGE_ref.vcf.gz, BRCA1_SGE_ref.gz.tbi and BRCA1_SGE_ref.hdr and reference BRCA1_SGE.csv

	"""

	# if the directory already exist then error
	assert not os.path.exists(test_work_dir), "The test directory {} already exists".format(test_work_dir)

	# otherwise make it in preparation for testing
	os.makedirs(test_work_dir)

	org_ref_hdr = os.path.join(working_files_dir, ref_vcf_hdr)
	org_ref_vcf = os.path.join(working_files_dir, ref_csv)
	org_ref_vcf_gz = os.path.join(working_files_dir, ref_vcf_gz)
	org_ref_vcf_gz_tbi = os.path.join(working_files_dir, ref_vcf_gz_tbi)
	
	new_ref_hdr = os.path.join(test_work_dir, ref_vcf_hdr)
	new_ref_vcf = os.path.join(test_work_dir, ref_csv)
	new_ref_vcf_gz = os.path.join(test_work_dir, ref_vcf_gz)
	new_ref_vcf_gz_tbi = os.path.join(test_work_dir, ref_vcf_gz_tbi)
	
	shutil.copyfile(org_ref_hdr, new_ref_hdr)
	shutil.copyfile(org_ref_vcf, new_ref_vcf)
	shutil.copyfile(org_ref_vcf_gz, new_ref_gz)
	shutil.copyfile(org_ref_vcf_gz_tbi, new_ref_vcf_gz_tbi)

	yield test_work_dir, new_ref_hdr, new_ref_vcf, new_ref_vcf_gz, new_ref_vcf_gz_tbi # This gets passed to the test function

	# This runs after the test function is complete to clean up
	# It deletes the test directory we created and anything it contains
	shutil.rmtree(test_work_dir)

