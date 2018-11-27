#!usr/bin/python3

import BRCA1_SGE_vcf_annotator
import pytest
import glob
import os
import shutil
import subprocess

# directories containing test datasets for the BRCA1_SGE_vcf_annotator program
BRCA1_rf_path = "/mnt/storage/home/zhengt/Competencies/CBI-6_Advanced_Clinical_Bioinformatics/BRCA1_test/"
working_files_directory = "/mnt/storage/home/zhengt/Competencies/CBI-6_Advanced_Clinical_Bioinformatics/BRCA1_SGE_vcf_annotator"

# vcf filenames used within the test directories
vcf_name1 = "sample1.vcf"
vcf_name2 = "sample2.vcf"
vcf_name3 = "sample3.vcf"
vcf_name4 = "sample4.vcf"

# annotated vcf filenames used within the test directories
annotated_vcf_name1 = "sample1.BRCA1_annotated.vcf"
annotated_vcf_name2 = "sample2.BRCA1_annotated.vcf"
annotated_vcf_name3 = "sample3.BRCA1_annotated.vcf"
annotated_vcf_name4 = "sample4.BRCA1_annotated.vcf"

# declare current directory
current_directory = os.getcwd()

# new directories to be created for each scenario
directory0 = os.path.abspath(os.path.join(current_directory, "pytest_sample0")) # test to see if you have write permission to current directory
test_directory = os.path.abspath(os.path.join(current_directory, "pytest_BRCA1")) # temp vcf and BRCA1_annotated.vcf for all test vcfs
working_files_directory
@pytest.fixture(scope="function")
def test_files_fixture():

	"""
    Set up test directory
    """

	# if the directory already exist then error
	assert not os.path.exists(test_directory), "The test directory {} already exists".format(directory1)

	# Otherwise make it in preparation for testing
	os.makedirs(test_directory)
	
	org_vcf_1 = os.path.join(BRCA1_rf_path, vcf_name1)
	org_vcf_2 = os.path.join(BRCA1_rf_path, vcf_name2)
	org_vcf_3 = os.path.join(BRCA1_rf_path, vcf_name3)
	org_vcf_4 = os.path.join(BRCA1_rf_path, vcf_name4)

	new_vcf_1 = os.path.join(test_directory, vcf_name1)
	new_vcf_2 = os.path.join(test_directory, vcf_name2)
	new_vcf_3 = os.path.join(test_directory, vcf_name3)
	new_vcf_4 = os.path.join(test_directory, vcf_name4)

	shutil.copyfile(org_vcf_1, new_vcf_1)
	shutil.copyfile(org_vcf_2, new_vcf_2)
	shutil.copyfile(org_vcf_3, new_vcf_3)
	shutil.copyfile(org_vcf_4, new_vcf_4)

	yield test_directory, new_vcf_1, new_vcf_2, new_vcf_3, new_vcf_4  # This gets passed to the test function

	# This runs after the test function is complete to clean up
	# It deletes the test directory we created and anything it contains
	shutil.rmtree(test_directory)

@pytest.fixture(scope="function")
def working_files_fixture(scope="function"):

	"""
	Set up working directory to store BRCA1_SGE_ref.csv, BRCA1_SGE_ref.txt, BRCA1_SGE_ref.vcf, BRCA1_SGE_ref.vcf.gz, BRCA1_SGE_ref.gz.tbi and BRCA1_SGE_ref.hdr files
	"""



def test_directory_does_not_exist(): 

	"""this checks whether the user has write permission to the current directory where the test directory will be stored"""
	assert os.makedirs(directory0), "The user do not have permission to create folders in current_directory"
	if os.path.exists(directory0):
		os.rmdir(directory0)

@pytest.mark.usefixtures("test_fixture")
def test_create_ref_vcf(test_fixture):

	""" make the BRCA1_ref_vcf using the create_ref_vcf function from BRCA1_SGE_vcf_annotator, and check the vcf have been generated correctly"""

	new_org_rf_path_1, new_bkup_rf_path_1, new_org_bam_1, new_org_md5_1, new_bkup_bam_1, new_bkup_md5_1 = scenario_1_fixture
		
	new_org_bam_list_1, new_org_md5_list_1 = md5sumscript.create_two_lists(new_org_rf_path_1)
	new_bkup_bam_list_1, new_bkup_md5_list_1 = md5sumscript.create_two_lists(new_bkup_rf_path_1)

	for bam_path in new_org_bam_list_1:
		assert bam_path.split(".")[-1] == "bam", "{} in original is NOT a bam file".format(bam_path)
	for bam_path in new_bkup_bam_list_1:
		assert bam_path.split(".")[-1] == "bam", "{} in backup is NOT a bam file".format(bam_path)
	for md5_path in new_org_md5_list_1:
		assert md5_path.split(".")[-1] == "md5" and md5_path.split(".")[-2] == "bam", "{} is NOT a bam associated md5 file".format(md5_path)
	for md5_path in new_bkup_md5_list_1:
		assert md5_path.split(".")[-1] == "md5" and md5_path.split(".")[-2] == "bam", "{} is NOT a bam associated md5 file".format(md5_path)

@pytest.mark.usefixtures("scenario_1_fixture", "scenario_9_fixture")
def test_check_md5_exist(scenario_1_fixture, scenario_9_fixture):

	"""using data from scenario 1 and 9, copy the test md5s and check the md5s exist in each test runfolder, 
	scenario 1 should return True as the md5 will be generated, whereas for scenario 9, it will be False, as the runfolder do not contain any md5"""

	new_org_rf_path_1, new_bkup_rf_path_1, new_org_bam_1, new_org_md5_1, new_bkup_bam_1, new_bkup_md5_1 = scenario_1_fixture
	new_org_rf_path_9, new_bkup_rf_path_9, new_org_bam_9, new_bkup_bam_9 = scenario_9_fixture

	org_checkfilepath_1 = md5sumscript.create_logfile(new_org_rf_path_1)
	bkup_checkfilepath_1 = md5sumscript.create_logfile(new_bkup_rf_path_1)
	org_checkfilepath_9 = md5sumscript.create_logfile(new_org_rf_path_9)
	
	assert md5sumscript.check_md5_exist(new_org_md5_1, org_checkfilepath_1), "ERROR, md5 DO NOT exist" # positive control
	assert md5sumscript.check_md5_exist(new_bkup_md5_1, bkup_checkfilepath_1), "ERROR, md5 DO NOT exist" # positive control
	assert not md5sumscript.check_md5_exist(os.path.join(new_org_rf_path_9, "sample13.bam.md5"), org_checkfilepath_9), "ERROR, md5 generated" # negative control

	md5_list = [new_org_md5_1, new_bkup_md5_1]

	# check that md5_check.txt log is created in the same directory as the md5, and that the check is recorded correctly

	for md5 in md5_list:

		try:
			f = open (org_checkfilepath_1, 'r')
			lines = f.readlines()
			for line in lines:
				if line.startswith("time of check:"):
					print ("time logged")
					return True
				elif "present" in line:
					print (md5, "OK")
					return True
				elif "missing" in line:
					print (md5, "missing status recorded")
					return True
				else:
					print ("ERROR, incorrect log recorded")
					return False

		except IOError:
			return IOError

@pytest.mark.usefixtures("scenario_7_fixture")
def test_create_md5(scenario_7_fixture):

	new_org_rf_path_7, new_bkup_rf_path_7, new_org_bam_7, new_org_md5_7, new_bkup_bam_7 = scenario_7_fixture

	"""using scenario 7, check that the program is able to generate the missing md5 files, and record the generation in a log file"""

	org_bam_list_7, org_md5_list_7 = md5sumscript.create_two_lists(new_org_rf_path_7)
	bkup_bam_list_7, bkup_md5_list_7 = md5sumscript.create_two_lists(new_bkup_rf_path_7)

	org_checkfilepath_7 = md5sumscript.create_logfile(new_org_rf_path_7)
	bkup_checkfilepath_7 = md5sumscript.create_logfile(new_bkup_rf_path_7)

	new_org_md5_list_7 = md5sumscript.create_md5(org_bam_list_7, new_org_rf_path_7, org_checkfilepath_7)
	new_bkup_md5_list_7 = md5sumscript.create_md5(bkup_bam_list_7, new_bkup_rf_path_7, bkup_checkfilepath_7)

	# check that the new_md5_lists contain md5s, and that they exist

	md5_list = new_org_md5_list_7 + new_bkup_md5_list_7

	new_md5s = []

	for md5 in md5_list:
		assert md5.endswith(".md5"), "file DO NOT end with .md5"

		if os.path.exists(md5):
			print ("{} generated OK".format(md5))
			new_md5s.append(md5)
			continue
		else:
			print ("md5 was NOT generated")

	assert md5_list == new_md5s, "some or all of the md5s in the new_md5_list are missing or NOT created"

	# check the log file md5_check.txt have recorded the new md5 files that have been generated, check that the time is also recorded

	for md5 in md5_list:

		logfile = "".join(md5.split("/")[:-1]) + "md5_check.txt"

		try:

			f = open (logfile, 'r')
			lines = f.readlines()
			for line in lines:
				if line.startswith("time of generation:"):
					print ("time of md5 generation logged")
				elif line.startswith("new md5 file") and md5 in line:
					print ("new md5 file recorded correctly in log")
				else:
					print ("ERROR, incorrect log recorded")
					continue

			f.close()
		
		except:
			None

@pytest.mark.usefixtures("scenario_2_fixture")
def test_create_logfile(scenario_2_fixture):

	"""using scenario 1, create log files that record the output of the md5sum -c program, and check the file format is correct"""
	new_org_rf_path_2, new_bkup_rf_path_2, new_org_bam_2, new_org_md5_2, new_bkup_bam_2, new_bkup_md5_2 = scenario_2_fixture
	rf_path_list = [new_org_rf_path_2, new_bkup_rf_path_2]
	
	# check the directory where the file is has write permission"""

	for rf_path in rf_path_list:
		assert os.access((rf_path), os.W_OK), "NO write permission for {}".format(rf_path)

	# check the directory where the file is has read permission"""
	
	for rf_path in rf_path_list:
		assert os.access((rf_path), os.R_OK), "NO read permission for {}".format(rf_path)

	# check the log file created exist"""

	for rf_path in rf_path_list:		
		checkfilepath = md5sumscript.create_logfile(rf_path)

		assert os.path.exists(checkfilepath), "{} DO NOT exist".format(checkfilepath)


@pytest.mark.usefixtures("scenario_3_fixture")
def test_check_md5_hash(scenario_3_fixture):

	new_org_rf_path_3, new_bkup_rf_path_3, new_org_bam_3, new_org_md5_3, new_bkup_bam_3, new_bkup_md5_3 = scenario_3_fixture
	rf_path_list = [new_org_rf_path_3, new_bkup_rf_path_3]

	"""using scenario 3, check the hash key is 32 characters long and consist alphannumeric symbols"""

	# first, create the log files
	org_checkfilepath_3 = md5sumscript.create_logfile(new_org_rf_path_3)
	bkup_checkfilepath_3 = md5sumscript.create_logfile(new_bkup_rf_path_3)

	# take the first 32 characters as checksum from the md5 files

	with open (new_org_md5_3, "r") as org_md5_file:
		org_checksum_3 = org_md5_file.readline().strip().split("  ")[0]
		print (org_checksum_3)
	with open (new_bkup_md5_3, "r") as bkup_md5_file:
		bkup_checksum_3 = bkup_md5_file.readline().strip().split("  ")[0]
		print (bkup_checksum_3)

	# check the hash using check_md5_hash function
		
	org_hash_check_3 = md5sumscript.check_md5_hash(org_checkfilepath_3, new_org_md5_3, org_checksum_3)
	bkup_hash_check_3 = md5sumscript.check_md5_hash(bkup_checkfilepath_3, new_bkup_md5_3, bkup_checksum_3)

	# 1. check the md5 hash is made of letters and numbers in both original and backup log files
	
	org_chars_3 = ""

	for char in org_checksum_3:
		if char.isalnum():
			org_chars_3 += char

	assert org_chars_3 == org_checksum_3, "characters in checksum NOT in correct format"

	if org_hash_check_3 is False:
		print ("check_md5_hash function (alphanumeric) NOT working for {}".format(new_org_md5_3))
	else:
		print ("check_md5_hash function (alphanumeric) working normally")

	bkup_chars_3 = ""

	for char in bkup_checksum_3:
		if char.isalnum():
			bkup_chars_3 += char

	assert not bkup_chars_3 == bkup_checksum_3, "characters in checksum in correct format"

	if bkup_hash_check_3 is True:
		print ("check_md5_hash function (alphanumeric) NOT working for {}".format(new_bkup_md5_3))
	else:
		print ("check_md5_hash function (alphanumeric) working normally, non-alphnumeric characters detected")


	# 2. check the md5 hash in the md5 file has correct length e.g. 32

	count = 0
	for char in org_checksum_3:
		count += 1

	org_len_3 = count

	count = 0
	for char in bkup_checksum_3:
		count += 1

	bkup_len_3 = count

	assert org_len_3 == 32, "length of md5 hash incorrect for {}".format(new_org_md5_3)

	if org_hash_check_3 is False:
		print ("check_md5_hash function (length) NOT working for {}".format(new_org_md5_3))
	else:
		print ("check_md5_hash function (length) working normally")

	assert bkup_len_3 != 32, "length of md5 hash correct for {}".format(new_bkup_md5_3)

	if bkup_hash_check_3 is True:
		print ("check_md5_hash function (length) NOT working for {}".format(new_bkup_md5_3))
	else:
		print ("check_md5_hash function (length) working normally, wrong length detected")

	# check the .chk file has recorded the error and normal logs in the correct format

	checkfilepaths = [org_checkfilepath_3, bkup_checkfilepath_3]

	for checkfilepath in checkfilepaths:

		with open (checkfilepath, "r") as checkfile:

			for line in checkfile:

				if "OK" and "have 32 characters and contain only letters or numbers" in line:
					print ("correct format recorded")

				elif "time of hash check:" in line:
					print ("hash check time logged")

				elif "ERROR" and "DO NOT have 32 characters or contain non-alphanumeric characters" in line:
					print ("correct format recorded")

				else:
					continue

@pytest.mark.usefixtures("scenario_1_fixture", "scenario_2_fixture", "scenario_3_fixture", "scenario_4_fixture", "scenario_5_fixture")
def test_check_filename(scenario_1_fixture, scenario_2_fixture, scenario_3_fixture, scenario_4_fixture, scenario_5_fixture):

	"""using scenario 1, 2, 3 and 5 as positive controls (where there is a mismatch),
	 and scenario 4 as negative control test to see whether check_filename function can identify discrepancy between bam filename prefix and the bam filename within the md5 file."""

	new_org_rf_path_1, new_bkup_rf_path_1, new_org_bam_1, new_org_md5_1, new_bkup_bam_1, new_bkup_md5_1 = scenario_1_fixture
	new_org_rf_path_2, new_bkup_rf_path_2, new_org_bam_2, new_org_md5_2, new_bkup_bam_2, new_bkup_md5_2 = scenario_2_fixture
	new_org_rf_path_3, new_bkup_rf_path_3, new_org_bam_3, new_org_md5_3, new_bkup_bam_3, new_bkup_md5_3 = scenario_3_fixture
	new_org_rf_path_4, new_bkup_rf_path_4, new_org_bam_4, new_org_md5_4, new_bkup_bam_4, new_bkup_md5_4 = scenario_4_fixture
	new_org_rf_path_5, new_bkup_rf_path_5, new_org_bam_5, new_org_md5_5, new_bkup_bam_5, new_bkup_md5_5 = scenario_5_fixture
	
	org_checkfilepath_1 = md5sumscript.create_logfile(new_org_rf_path_1)
	bkup_checkfilepath_1 = md5sumscript.create_logfile(new_bkup_rf_path_1)

	org_checkfilepath_2 = md5sumscript.create_logfile(new_org_rf_path_2)
	bkup_checkfilepath_2 = md5sumscript.create_logfile(new_bkup_rf_path_2)

	org_checkfilepath_3 = md5sumscript.create_logfile(new_org_rf_path_3)
	bkup_checkfilepath_3 = md5sumscript.create_logfile(new_bkup_rf_path_3)

	org_checkfilepath_4 = md5sumscript.create_logfile(new_org_rf_path_4)
	bkup_checkfilepath_4 = md5sumscript.create_logfile(new_bkup_rf_path_4)

	org_checkfilepath_5 = md5sumscript.create_logfile(new_org_rf_path_5)
	bkup_checkfilepath_5 = md5sumscript.create_logfile(new_bkup_rf_path_5)

	org_bam_filename_1 = new_org_bam_1.split("/")[-1]
	bkup_bam_filename_1 = new_bkup_bam_1.split("/")[-1]
	org_bam_filename_2 = new_org_bam_2.split("/")[-1]
	bkup_bam_filename_2 = new_bkup_bam_2.split("/")[-1]
	org_bam_filename_3 = new_org_bam_3.split("/")[-1]
	bkup_bam_filename_3 = new_bkup_bam_3.split("/")[-1]
	org_bam_filename_4 = new_org_bam_4.split("/")[-1]
	bkup_bam_filename_4 = new_bkup_bam_4.split("/")[-1]
	org_bam_filename_5 = new_org_bam_5.split("/")[-1]
	bkup_bam_filename_5 = new_bkup_bam_5.split("/")[-1]
	
	with open (new_org_md5_1, "r") as md5_file:
		org_bam_in_md5_1 = md5_file.readline().strip().split("  ")[-1]
	with open (new_bkup_md5_1, "r") as md5_file:
		bkup_bam_in_md5_1 = md5_file.readline().strip().split("  ")[-1]
	with open (new_org_md5_2, "r") as md5_file:
		org_bam_in_md5_2 = md5_file.readline().strip().split("  ")[-1]
	with open (new_bkup_md5_2, "r") as md5_file:
		bkup_bam_in_md5_2 = md5_file.readline().strip().split("  ")[-1]
	with open (new_org_md5_3, "r") as md5_file:
		org_bam_in_md5_3 = md5_file.readline().strip().split("  ")[-1]
	with open (new_bkup_md5_3, "r") as md5_file:
		bkup_bam_in_md5_3 = md5_file.readline().strip().split("  ")[-1]
	with open (new_org_md5_4, "r") as md5_file:
		org_bam_in_md5_4 = md5_file.readline().strip().split("  ")[-1]
	with open (new_bkup_md5_4, "r") as md5_file:
		bkup_bam_in_md5_4 = md5_file.readline().strip().split("  ")[-1]
	with open (new_org_md5_5, "r") as md5_file:
		org_bam_in_md5_5 = md5_file.readline().strip().split("  ")[-1]
	with open (new_bkup_md5_5, "r") as md5_file:
		bkup_bam_in_md5_5 = md5_file.readline().strip().split("  ")[-1]

	# for scenarios 4, the bam filenames prefix should match that which is next to the hash in the file.
	
	assert check_filename(org_checkfilepath_4, org_bam_filename_4, org_bam_in_md5_4, new_org_md5_4), "check_filename function NOT working"
	assert check_filename(bkup_checkfilepath_4, bkup_bam_filename_4, bkup_bam_in_md5_4, new_bkup_md5_4), "check_filename function NOT working"

	# for scenarios 1, 2, 3 and 5, the bam filenames prefix should NOT match that which is in the file.

	assert not check_filename(org_checkfilepath_1, org_bam_filename_1, org_bam_in_md5_1, new_org_md5_1), "check_filename function NOT working"
	assert not check_filename(bkup_checkfilepath_1, bkup_bam_filename_1, bkup_bam_in_md5_1, new_bkup_md5_1), "check_filename function NOT working"
	assert not check_filename(org_checkfilepath_2, org_bam_filename_2, org_bam_in_md5_2, new_org_md5_2), "check_filename function NOT working"
	assert not check_filename(bkup_checkfilepath_2, bkup_bam_filename_2, bkup_bam_in_md5_2, new_bkup_md5_2), "check_filename function NOT working"
	assert not check_filename(org_checkfilepath_3, org_bam_filename_3, org_bam_in_md5_3, new_org_md5_3), "check_filename function NOT working"
	assert not check_filename(bkup_checkfilepath_3, bkup_bam_filename_3, bkup_bam_in_md5_3, new_bkup_md5_3), "check_filename function NOT working"
	assert not check_filename(org_checkfilepath_5, org_bam_filename_5, org_bam_in_md5_5, new_org_md5_5), "check_filename function NOT working"
	assert not check_filename(bkup_checkfilepath_5, bkup_bam_filename_5, bkup_bam_in_md5_5, new_bkup_md5_5), "check_filename function NOT working"

	# the .chk file records the matches and mismatches in terms of filename, hence this test to see if the filename that matches have been recorded correctly

	checkfilepaths = [org_checkfilepath_1, bkup_checkfilepath_1, org_checkfilepath_2, bkup_checkfilepath_2, org_checkfilepath_3, bkup_checkfilepath_3, org_checkfilepath_4, bkup_checkfilepath_4, org_checkfilepath_5, bkup_checkfilepath_5]
	
	for checkfilepath in checkfilepaths:

		with open (checkfilepath, "r") as checkfile:
			for line in checkfile:

				if line.startswith(new_org_md5_1) and "file is being checked" in line:
					print ("{} file has been checked".format(new_org_md5_1))

				elif line.startswith(new_org_md5_2) and "file is being checked" in line:
					print ("{} file has been checked".format(new_org_md5_2))

				elif line.startswith(new_org_md5_3) and "file is being checked" in line:
					print ("{} file has been checked".format(new_org_md5_3))

				elif line.startswith(new_org_md5_4) and "file is being checked" in line:
					print ("{} file has been checked".format(new_org_md5_4))

				elif line.startswith(new_org_md5_5) and "file is being checked" in line:
					print ("{} file has been checked".format(new_org_md5_5))

				elif line.startswith(new_bkup_md5_1) and "file is being checked" in line:
					print ("{} file has been checked".format(new_bkup_md5_1))

				elif line.startswith(new_bkup_md5_2) and "file is being checked" in line:
					print ("{} file has been checked".format(new_bkup_md5_2))

				elif line.startswith(new_bkup_md5_3) and "file is being checked" in line:
					print ("{} file has been checked".format(new_bkup_md5_3))

				elif line.startswith(new_bkup_md5_4) and "file is being checked" in line:
					print ("{} file has been checked".format(new_bkup_md5_4))

				elif line.startswith(new_bkup_md5_5) and "file is being checked" in line:
					print ("{} file has been checked".format(new_bkup_md5_5))

				elif "OK" and "match" and "in" in line:
					print ("correct format recorded")

				elif "time of filename check:" in line:
					print ("filename check time logged")

				elif "ERROR" and "DO NOT match" and "in" in line:
					print ("correct format recorded")

				else:
					"Error, wrong record in log file"

@pytest.mark.usefixtures("scenario_1_fixture", "scenario_4_fixture", "scenario_7_fixture")
def test_check_md5(scenario_1_fixture, scenario_4_fixture, scenario_7_fixture):
	
	"""using scenario 1, 4 and 7, generate the.chk log files and .md5 files. test the check_md5 function will take the hash and filename from within the md5 file and add it to the check_dict, 
	assesses the md5sum -c output is recorded correctly in .chk log file"""	

	new_org_rf_path_1, new_bkup_rf_path_1, new_org_bam_1, new_org_md5_1, new_bkup_bam_1, new_bkup_md5_1 = scenario_1_fixture
	new_org_rf_path_4, new_bkup_rf_path_4, new_org_bam_4, new_org_md5_4, new_bkup_bam_4, new_bkup_md5_4 = scenario_4_fixture
	new_org_rf_path_7, new_bkup_rf_path_7, new_org_bam_7, new_org_md5_7, new_bkup_bam_7 = scenario_7_fixture

	md5_list = [new_org_md5_1, new_org_md5_4, new_org_md5_7, new_bkup_md5_1, new_bkup_md5_4]

	org_checkfilepath_1 = md5sumscript.create_logfile(new_org_rf_path_1)
	bkup_checkfilepath_1 = md5sumscript.create_logfile(new_bkup_rf_path_1)
	org_checkfilepath_4 = md5sumscript.create_logfile(new_org_rf_path_4)
	bkup_checkfilepath_4 = md5sumscript.create_logfile(new_bkup_rf_path_4)
	org_checkfilepath_7 = md5sumscript.create_logfile(new_org_rf_path_7)

	checkfilepaths = [org_checkfilepath_1, bkup_checkfilepath_1, org_checkfilepath_4, bkup_checkfilepath_4, org_checkfilepath_7]

	new_org_bam_list_1, new_org_md5_list_1 = md5sumscript.create_two_lists(new_org_rf_path_1)
	new_bkup_bam_list_1, new_bkup_md5_list_1 = md5sumscript.create_two_lists(new_bkup_rf_path_1)
	new_org_bam_list_4, new_org_md5_list_4 = md5sumscript.create_two_lists(new_org_rf_path_4)
	new_bkup_bam_list_4, new_bkup_md5_list_4 = md5sumscript.create_two_lists(new_bkup_rf_path_4)
	new_org_bam_list_7, new_org_md5_list_7 = md5sumscript.create_two_lists(new_org_rf_path_7)

	org_check_dict_1, org_md5_bam_dict_1 = md5sumscript.check_md5(org_checkfilepath_1, new_org_md5_list_1)
	bkup_check_dict_1, bkup_md5_bam_dict_1 = md5sumscript.check_md5(bkup_checkfilepath_1, new_bkup_md5_list_1)
	org_check_dict_4, org_md5_bam_dict_4 = md5sumscript.check_md5(org_checkfilepath_4, new_org_md5_list_4)
	bkup_check_dict_4, bkup_md5_bam_dict_4 = md5sumscript.check_md5(bkup_checkfilepath_4, new_bkup_md5_list_4)
	org_check_dict_7, org_md5_bam_dict_7 = md5sumscript.check_md5(org_checkfilepath_7, new_org_md5_list_7)

	# this checks the check_dict created by the check_md5 function contain the correct data for each scenario

	org_bkup_dict_1 = {"sample1.bam.md5" : "15f1891487fd6b26f752b50468adb49f"}
	org_bkup_dict_4 = {"sample1.bam.md5" : "8620b0a9852e248c88b3f7ed30e73d01"}
	org_dict_7 = {"sample11.bam.md5" : "fa0aaaeaf163a93c9e4b96d757079032"}

	for key, value in org_check_dict_1.items():
		assert key == org_bkup_dict_1.keys()[0]
		assert value == org_bkup_dict_1.values()[0]

	for key, value in bkup_check_dict_1.items():
		assert key == org_bkup_dict_1.keys()[0]
		assert value == org_bkup_dict_1.values()[0]

	for key, value in org_check_dict_4.items(): 
		assert key == org_bkup_dict_4.keys()[0]
		assert value == org_bkup_dict_4.values()[0]

	for key, value in bkup_check_dict_4.items():
		assert key == org_bkup_dict_4.keys()[0]
		assert value == org_bkup_dict_4.values()[0]

	for key, value in org_check_dict_7.items():
		assert key == org_dict_7.keys()[0]
		assert value == org_dict_7.values()[0]

	for key, value in org_check_dict_7.items():
		assert key == org_dict_7.keys()[0]
		assert value == org_dict_7.values()[0]

	# this checks the md5sum -c output for each md5 file.

	org_bam_filename_1 = new_org_bam_1.split("/")[-1]
	bkup_bam_filename_1 = new_bkup_bam_1.split("/")[-1]
	org_bam_filename_4 = new_org_bam_4.split("/")[-1]
	bkup_bam_filename_4 = new_bkup_bam_4.split("/")[-1]
	org_bam_filename_7 = new_org_bam_7.split("/")[-1]

	for checkfilepath in checkfilepaths:

		with open(checkfilepath, "r") as checkfile:

			for line in checkfile.readlines():

				if line.startswith("md5sum:") and "WARNING:" in line and "1 listed file could not be read" in line:
					print ("md5sum -c program executed correctly, with WARNING message stating file CANNOT BE READ")

				elif line.startswith("md5sum:") and "WARNING:" in line and "1 computed checksum did NOT match" in line:
					print ("md5sum -c program executed correctly, with WARNING message stating file checksum DO NOT match")

				elif line.startswith(org_bam_filename_4) and "OK" in line:
					print ("{} recorded in log file, and hash match".format(org_bam_filename_1))

				elif line.startswith(bkup_bam_filename_4) and "FAILED" in line:
					print ("{} recorded in log file, and hash DO NOT match".format(bkup_bam_filename_4))

				elif line.startswith("md5sum:") and (org_bam_filename_1) in line and "No such file or directory" in line:
					print ("md5sum -c program executed for {}, bam file stated does NOT exist". format(org_bam_filename_1))

				elif line.startswith(org_bam_filename_1) and "FAILED open or read" in line:
					print ("md5sum -c program executed for {}, bam file cannot be opened or read". format(org_bam_filename_1))

				elif line.startswith("md5sum:") and (bkup_bam_filename_1) in line and "No such file or directory" in line:
					print ("md5sum -c program executed for {}, bam file stated does NOT exist". format(bkup_bam_filename_1))

				elif line.startswith(bkup_bam_filename_1) and "FAILED open or read" in line:
					print ("md5sum -c program executed for {}, bam file cannot be opened or read". format(bkup_bam_filename_1))

				elif line.startswith(org_bam_filename_7) and "OK" in line:
					print ("{} recorded in log file, and hash match".format(org_bam_filename_7))

				elif line.startswith(org_bam_filename_7) and "OK" in line:
					print ("{} recorded in log file, and hash match".format(org_bam_filename_1))
				else:
					continue

@pytest.mark.usefixtures("scenario_7_fixture")
def test_check_filename(scenario_7_fixture):

	"""test that given filenames that should be present in dictionary are, and vice versa using scenario 7 as example"""

	new_org_rf_path_7, new_bkup_rf_path_7, new_org_bam_7, new_org_md5_7, new_bkup_bam_7 = scenario_7_fixture

	new_org_bam_list_7, new_org_md5_list_7 = md5sumscript.create_two_lists(new_org_rf_path_7)
	new_bkup_bam_list_7, new_bkup_md5_list_7 = md5sumscript.create_two_lists(new_bkup_rf_path_7)
       
	org_checkfilepath_7 = md5sumscript.create_logfile(new_org_rf_path_7)
	bkup_checkfilepath_7 = md5sumscript.create_logfile(new_bkup_rf_path_7)

	org_check_dict_7, org_md5_bam_dict_7 = md5sumscript.check_md5(org_checkfilepath_7, new_org_md5_list_7)
	bkup_check_dict_7, bkup_md5_bam_dict_7 = md5sumscript.check_md5(bkup_checkfilepath_7, new_bkup_md5_list_7)

	# for new_org_md5_7 created above, take only the md5_filename from the file path
		
	md5_filename = new_org_md5_7.split("/")[-1]

	# md5_filename should be in check_dict["storage"], but not check_dict["archive"]

	assert md5sumscript.check_md5filename(new_org_rf_path_7, md5_filename, org_check_dict_7), "md5 NOT in storage, check create_check_dict function"
	assert not md5sumscript.check_md5filename(new_bkup_rf_path_7, md5_filename, bkup_check_dict_7), "md5 in archive check create_check_dict function"

@pytest.mark.usefixtures("scenario_6_fixture", "scenario_4_fixture")
def test_compare_md5hash(scenario_6_fixture, scenario_4_fixture):

	"""using scenario 6 and 8, check if the same md5 file has matching hash (scenario 8) or not (scenario 6) between original and backup, 
	the compare_md5hash should return all matching hash as True and none matching as False, if not raise an exception"""
	
	new_org_rf_path_6, new_bkup_rf_path_6, new_org_bam_6, new_org_md5_6, new_bkup_bam_6, new_bkup_md5_6 = scenario_6_fixture
	new_org_rf_path_4, new_bkup_rf_path_4, new_org_bam_4, new_org_md5_4, new_bkup_bam_4, new_bkup_md5_4 = scenario_4_fixture

	new_org_bam_list_6, new_org_md5_list_6 = md5sumscript.create_two_lists(new_org_rf_path_6)
	new_bkup_bam_list_6, new_bkup_md5_list_6 = md5sumscript.create_two_lists(new_bkup_rf_path_6)

	new_org_bam_list_4, new_org_md5_list_4 = md5sumscript.create_two_lists(new_org_rf_path_4)
	new_bkup_bam_list_4, new_bkup_md5_list_4 = md5sumscript.create_two_lists(new_bkup_rf_path_4)


	org_checkfilepath_6 = md5sumscript.create_logfile(new_org_rf_path_6)
	bkup_checkfilepath_6 = md5sumscript.create_logfile(new_bkup_rf_path_6)

	org_checkfilepath_4 = md5sumscript.create_logfile(new_org_rf_path_4)
	bkup_checkfilepath_4 = md5sumscript.create_logfile(new_bkup_rf_path_4)

	org_check_dict_6, org_md5_bam_dict_6 = md5sumscript.check_md5(org_checkfilepath_6, new_org_md5_list_6)
	bkup_check_dict_6, bkup_md5_bam_dict_6 = md5sumscript.check_md5(bkup_checkfilepath_6, new_bkup_md5_list_6)
	org_check_dict_4, org_md5_bam_dict_4 = md5sumscript.check_md5(org_checkfilepath_4, new_org_md5_list_4)
	bkup_check_dict_4, bkup_md5_bam_dict_4 = md5sumscript.check_md5(bkup_checkfilepath_4, new_bkup_md5_list_4)

	# for matching hashes, check that the compare_md5hash function returns True

	for md5 in new_org_md5_list_4:
		md5_filename = md5.strip().split("/")[-1]
		original_hash = org_check_dict_4[md5_filename]
		backup_hash = bkup_check_dict_4[md5_filename]
		assert md5sumscript.compare_md5hash(original_hash, backup_hash), "hash matching, but output is False, compare_md5hash function NOT working"

	for md5 in new_bkup_md5_list_4:
		md5_filename = md5.strip().split("/")[-1]
		original_hash = org_check_dict_4[md5_filename]
		backup_hash = bkup_check_dict_4[md5_filename]
		assert md5sumscript.compare_md5hash(original_hash, backup_hash), "hash matching, but output is False, compare_md5hash function NOT working"

	# for mismatch hashes,  check that the compare_md5hash function returns False

	for md5 in new_org_md5_list_6:
		md5_filename = md5.strip().split("/")[-1]
		original_hash = org_check_dict_6[md5_filename]
		backup_hash = bkup_check_dict_6[md5_filename]
		assert not md5sumscript.compare_md5hash(original_hash, backup_hash), "hash matching, but output is True, compare_md5hash function NOT working"
	
	for md5 in new_bkup_md5_list_6:
		md5_filename = md5.strip().split("/")[-1]
		original_hash = org_check_dict_6[md5_filename]
		backup_hash = bkup_check_dict_6[md5_filename]
		assert not md5sumscript.compare_md5hash(original_hash, backup_hash), "hash not matching, but output is True, compare_md5hash function NOT working"

@pytest.mark.usefixtures("scenario_6_fixture", "scenario_4_fixture")
def test_check_hash_exist(scenario_6_fixture, scenario_4_fixture):

	"""check the hash in one sub-dictionary belong in another, e.g. check_dict["storage"]["X123456"] = "abcdefg123456" is also in check_dict["archive"], 
	again this should match the output from the check_hash_exist function"""

	new_org_rf_path_6, new_bkup_rf_path_6, new_org_bam_6, new_org_md5_6, new_bkup_bam_6, new_bkup_md5_6 = scenario_6_fixture
	new_org_rf_path_4, new_bkup_rf_path_4, new_org_bam_4, new_org_md5_4, new_bkup_bam_4, new_bkup_md5_4 = scenario_4_fixture

	new_org_bam_list_6, new_org_md5_list_6 = md5sumscript.create_two_lists(new_org_rf_path_6)
	new_bkup_bam_list_6, new_bkup_md5_list_6 = md5sumscript.create_two_lists(new_bkup_rf_path_6)

	new_org_bam_list_4, new_org_md5_list_4 = md5sumscript.create_two_lists(new_org_rf_path_4)
	new_bkup_bam_list_4, new_bkup_md5_list_4 = md5sumscript.create_two_lists(new_bkup_rf_path_4)

	org_checkfilepath_6 = md5sumscript.create_logfile(new_org_rf_path_6)
	bkup_checkfilepath_6 = md5sumscript.create_logfile(new_bkup_rf_path_6)

	org_checkfilepath_4 = md5sumscript.create_logfile(new_org_rf_path_4)
	bkup_checkfilepath_4 = md5sumscript.create_logfile(new_bkup_rf_path_4)

	org_check_dict_6, org_md5_bam_dict_6 = md5sumscript.check_md5(org_checkfilepath_6, new_org_md5_list_6)
	bkup_check_dict_6, bkup_md5_bam_dict_6 = md5sumscript.check_md5(bkup_checkfilepath_6, new_bkup_md5_list_6)
	org_check_dict_4, bkup_md5_bam_dict_4 = md5sumscript.check_md5(org_checkfilepath_4, new_org_md5_list_4)
	bkup_check_dict_4, bkup_md5_bam_dict_4 = md5sumscript.check_md5(bkup_checkfilepath_4, new_bkup_md5_list_4)

	# for hashes that are present in opposite runfolder

	for md5 in new_org_md5_list_4:
		md5_filename = md5.strip().split("/")[-1]
		assert md5sumscript.check_hash_exist(md5_filename, org_check_dict_4, bkup_check_dict_4), "check_hash_exist function NOT working"
	
	for md5 in new_bkup_md5_list_4:
		md5_filename = md5.strip().split("/")[-1]
		assert md5sumscript.check_hash_exist(md5_filename, bkup_check_dict_4, org_check_dict_4), "check_hash_exist function NOT working"

	# for hashes that are NOT present in opposite runfolder
	for md5 in new_org_md5_list_6:
		md5_filename = md5.strip().split("/")[-1]
		assert not md5sumscript.check_hash_exist(md5_filename, org_check_dict_6, bkup_check_dict_6), "check_hash_exist function NOT working"
	
	for md5 in new_bkup_md5_list_6:
		md5_filename = md5.strip().split("/")[-1]
		assert not md5sumscript.check_hash_exist(md5_filename, bkup_check_dict_6, org_check_dict_6), "check_hash_exist function NOT working"

@pytest.mark.usefixtures("scenario_1_fixture", "scenario_2_fixture", "scenario_6_fixture", "scenario_8_fixture")		
def test_check_org_bkup(scenario_1_fixture, scenario_2_fixture, scenario_6_fixture, scenario_8_fixture):

	"""using scenario 1, 2, 6, and 8, check that matches and mismatches of md5 filename and hashes between original and back up runfolder are recorded in
	the four dictionaries md5_matches, md5_mismatches, md5_in_bkup_not_not_org, md5_in_org_not_bkup."""

	new_org_rf_path_1, new_bkup_rf_path_1, new_org_bam_1, new_org_md5_1, new_bkup_bam_1, new_bkup_md5_1 = scenario_1_fixture
	new_org_rf_path_2, new_bkup_rf_path_2, new_org_bam_2, new_org_md5_2, new_bkup_bam_2, new_bkup_md5_2 = scenario_2_fixture
	new_org_rf_path_6, new_bkup_rf_path_6, new_org_bam_6, new_org_md5_6, new_bkup_bam_6, new_bkup_md5_6 = scenario_6_fixture
	new_org_rf_path_8, new_bkup_rf_path_8, new_org_bam_8, new_org_md5_8, new_bkup_bam_8, new_bkup_md5_8 = scenario_8_fixture

	new_org_bam_list_1, new_org_md5_list_1 = md5sumscript.create_two_lists(new_org_rf_path_1)
	new_bkup_bam_list_1, new_bkup_md5_list_1 = md5sumscript.create_two_lists(new_bkup_rf_path_1)

	new_org_bam_list_2, new_org_md5_list_2 = md5sumscript.create_two_lists(new_org_rf_path_2)
	new_bkup_bam_list_2, new_bkup_md5_list_2 = md5sumscript.create_two_lists(new_bkup_rf_path_2)

	new_org_bam_list_6, new_org_md5_list_6 = md5sumscript.create_two_lists(new_org_rf_path_6)
	new_bkup_bam_list_6, new_bkup_md5_list_6 = md5sumscript.create_two_lists(new_bkup_rf_path_6)

	new_org_bam_list_8, new_org_md5_list_8 = md5sumscript.create_two_lists(new_org_rf_path_8)
	new_bkup_bam_list_8, new_bkup_md5_list_8 = md5sumscript.create_two_lists(new_bkup_rf_path_8)

	org_checkfilepath_1 = md5sumscript.create_logfile(new_org_rf_path_1)
	bkup_checkfilepath_1 = md5sumscript.create_logfile(new_bkup_rf_path_1)

	org_checkfilepath_2 = md5sumscript.create_logfile(new_org_rf_path_2)
	bkup_checkfilepath_2 = md5sumscript.create_logfile(new_bkup_rf_path_2)

	org_checkfilepath_6 = md5sumscript.create_logfile(new_org_rf_path_6)
	bkup_checkfilepath_6 = md5sumscript.create_logfile(new_bkup_rf_path_6)

	org_checkfilepath_8 = md5sumscript.create_logfile(new_org_rf_path_8)
	bkup_checkfilepath_8 = md5sumscript.create_logfile(new_bkup_rf_path_8)

	checkfilepaths = [org_checkfilepath_1, bkup_checkfilepath_1, org_checkfilepath_2, bkup_checkfilepath_2, org_checkfilepath_6, bkup_checkfilepath_6, org_checkfilepath_8, bkup_checkfilepath_8]

	# create check_dict with original and backup check file path and md5_lists

	org_check_dict_1, org_md5_bam_dict_1 = md5sumscript.check_md5(org_checkfilepath_1, new_org_md5_list_1)
	bkup_check_dict_2, bkup_md5_bam_dict_2 = md5sumscript.check_md5(bkup_checkfilepath_2, new_bkup_md5_list_2)
	org_check_dict_6, org_md5_bam_dict_6 = md5sumscript.check_md5(org_checkfilepath_6, new_org_md5_list_6)
	bkup_check_dict_6, bkup_md5_bam_dict_6 = md5sumscript.check_md5(bkup_checkfilepath_6, new_bkup_md5_list_6)
	org_check_dict_8, org_md5_bam_dict_8 = md5sumscript.check_md5(org_checkfilepath_8, new_org_md5_list_8)
	bkup_check_dict_8, bkup_md5_bam_dict_8 = md5sumscript.check_md5(bkup_checkfilepath_8, new_bkup_md5_list_8)
	
	# create check_dict with cross referenced original and backup runfolder paths and check_dicts

	md5_matches_1_2, md5_mismatches_1_2, md5_in_org_not_bkup_1_2, md5_in_bkup_not_org_1_2 = md5sumscript.check_org_bkup(new_org_rf_path_1, new_bkup_rf_path_2, org_checkfilepath_1, bkup_checkfilepath_2, org_check_dict_1, bkup_check_dict_2)
	md5_matches_6, md5_mismatches_6, md5_in_org_not_bkup_6, md5_in_bkup_not_org_6 = md5sumscript.check_org_bkup(new_org_rf_path_6, new_bkup_rf_path_6, org_checkfilepath_6, bkup_checkfilepath_6, org_check_dict_6, bkup_check_dict_6)
	md5_matches_8, md5_mismatches_8, md5_in_org_not_bkup_8, md5_in_bkup_not_org_8 = md5sumscript.check_org_bkup(new_org_rf_path_8, new_bkup_rf_path_8, org_checkfilepath_8, bkup_checkfilepath_8, org_check_dict_8, bkup_check_dict_8)
	md5_matches_6_8, md5_mismatches_6_8, md5_in_org_not_bkup_6_8, md5_in_bkup_not_org_6_8 = md5sumscript.check_org_bkup(new_org_rf_path_6, new_bkup_rf_path_8, org_checkfilepath_6, bkup_checkfilepath_8, org_check_dict_6, bkup_check_dict_8)

	# check md5_matches dictionary content for each md5 and associated hash (same md5 filename, same hash)

	for key, values in md5_matches_8["original"].items():
		assert key == "sample12.bam.md5"
		assert values == "38f7c7e572fe7aeffe4fcf5728c31f5d"
	for key, values in md5_matches_8["backup"].items():
		assert key == "sample12.bam.md5"
		assert values == "38f7c7e572fe7aeffe4fcf5728c31f5d"

	# check md5_mismatches dictionary content (same md5 filename, different hash)

	for key, values in md5_mismatches_6["original"].items():
		assert key == "sample3.bam.md5"
		assert values == "a5119b252ab16786e583b550fab714ac"
	for key, values in md5_mismatches_6["backup"].items():
		assert key == "sample3.bam.md5"
		assert values == "a5119b252ab16786e583b550fab714ad"

	# check md5_in_bkup_not_org dictionary content (different md5 filename, same hash)

	for key, values in md5_in_bkup_not_org_1_2["hash present"].items():
		assert key == "sample2.bam.md5"
		assert values == "15f1891487fd6b26f752b50468adb49f"

	for key, values in md5_in_bkup_not_org_6_8["no hash match"].items():
		assert key == "sample12.bam.md5"
		assert values == "38f7c7e572fe7aeffe4fcf5728c31f5d"
	
	# check md5_in_org_not_org dictionary content (different md5 filename, same hash)

	for key, values in md5_in_org_not_bkup_1_2["hash present"].items():
		assert key == "sample1.bam.md5"
		assert values == "15f1891487fd6b26f752b50468adb49f"

	for key, values in md5_in_org_not_bkup_6_8["no hash match"].items():
		assert key == "sample3.bam.md5"
		assert values == "a5119b252ab16786e583b550fab714ac"

	# check output of scenario 1,2,6 and 8 in org_bkup_check.txt 

	for checkfilepath in checkfilepaths:

		with open(checkfilepath, "r") as checkfile:

			for line in checkfile.readlines():

				if line.startswith("OK") and "sample12.bam.md5" in line and "has matching hash"in line:
					print ("md5 filename and hash match")

				elif line.startswith("sample12.bam.md5") and "present in backup and original runfolder" in line:
					print ("md5 filename found in both backup and original runfolder")

				elif line.startswith("ERROR") and "a5119b252ab16786e583b550fab714ad" in line and "a5119b252ab16786e583b550fab714ac" in line and "sample3.bam.md5" in line and "DO NOT MATCH" in line:
					print ("mismatched hash recorded correctly")

				elif line.startswith("ERROR") and "not present in original runfolder" in line:
					print ("cross referenced samples recorded correctly in backup runfolder")

				elif line.startswith("{}".format("sample1.bam.md5" or "sample2.bam.md5")) and "hash in backup folder found in original folder" in line:
					print ("cross referenced samples with matching hash recorded correctly")

				elif line.startswith("sample12.bam.md5" or "sample3.bam.md5") and "hash in backup folder NOT found in original folder" in line:
					print ("cross referenced sample not matching filename or recorded correctly")

				elif line.startswith("ERROR") and "not present in backup runfolder" in line:
					print ("cross referenced samples recorded correctly in original runfolder")

				elif line.startswith("time of check:"):
					print ("time logged")

				elif "being checked" in line:
					print ("check recorded in log file")

				else:
					continue
