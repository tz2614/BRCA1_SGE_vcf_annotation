#!usr/bin/python2

import sys
import os
import subprocess
import datetime
import pprint as pp
import vcf

"""create a script that checks which CP samples contain BRCA1 variants in the BRCA1_SGE excel spreadsheet"""


def create_BRCA1_SGE_dict(csv_file):

	#create a dictionary containig containing the CHROM, POS, REF and ALT

	dir_path = os.path.dirname(csv_file)
	#BRCA1_SGE_ref = os.path.join(dir_path, "BRCA1_SGE_ref.txt")

	#if os.path.exists(BRCA1_SGE_ref):
		#return BRCA1_SGE_ref

	BRCA1_SGE_dict = {}
	
	with open (csv_file, "r") as csv_file:

		print ("creating BRCA1_SGE_ref dictionary")

		lines = [line.split(',') for line in csv_file]

		for index, fields in enumerate(lines):

			if index < 2:
				continue

			elif index == 2:
				
				for column_index, field in enumerate(fields):

					if field == "chromosome":
						chrom_num = column_index
					elif field == "position (hg19)":
						pos = column_index
					elif field == "reference":
						ref = column_index
					elif field == "alt":
						alt = column_index
					else:
						continue

			elif index > 2:

				CHROM = fields[chrom_num]
				POS = fields[pos]
				REF = fields[ref]
				ALT = fields[alt]
				if CHROM not in BRCA1_SGE_dict.keys():
					BRCA1_SGE_dict[CHROM] = {}
				if POS not in BRCA1_SGE_dict[CHROM].keys():
					BRCA1_SGE_dict[CHROM][POS] = {}
				if REF not in BRCA1_SGE_dict[CHROM][POS].keys():
					BRCA1_SGE_dict[CHROM][POS][REF] = []
				if ALT not in BRCA1_SGE_dict[CHROM][POS][REF]:
					BRCA1_SGE_dict[CHROM][POS][REF].append(ALT)

	
	print ("BRCA1_SGE_dict generated")
	return BRCA1_SGE_dict

def get_list_of_C01_vcfs(runfolder):

	vcf_list = []

	for (dirpath, dirnames, filenames) in os.walk(runfolder):
		for filename in filenames:
			if filename.endswith(".vcf") and filename.startswith("C01") and filename not in vcf_list:
				vcf_list.extend(filenames)
			else:
				continue
				
	print (sorted(vcf_list))
	return sorted(vcf_list)

def compare_BRCA1_SGE_to_vcf(vcf_list, BRCA1_SGE_dict):

	"""check which CP file contain variants from the supplementary table 1.csv file"""

	SGE_positive_vcfs = []

	for vcf_file in vcf_list:

		# check each vcf record against the BRCA1_SGE_dict
	
		vcf_data = vcf.Reader(open(vcf_file, "rb"))

		print ("analysing {}".format(vcf_file))

		for record in vcf_data:
			try:
				if BRCA1_SGE_dict[record.CHROM][record.POS][record.REF] == record.ALT[0]:
					print ("{} contain the following BRCA1_SGE_variants:".format(vcf_file))
					print (record)

			except KeyError:
				continue

		SGE_positive_vcfs.append(vcf_file)

	return SGE_positive_vcfs

def main(runfolder, csv_file):

	# check that the annotated_vcf and excel_table exist as a data source
	assert os.path.exists(csv_file), "{} DO NOT exists".format(csv_file)

	BRCA1_SGE_dict = create_BRCA1_SGE_dict(csv_file)
	vcf_list = get_list_of_C01_vcfs(runfolder)
	compare_BRCA1_SGE_to_vcf(vcf_list, BRCA1_SGE_dict)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])

	"""sys.argv1 = /path/to/CP runfolder, sys.argv2 = supplementary table 1.csv"""
