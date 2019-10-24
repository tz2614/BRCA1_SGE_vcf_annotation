#!usr/bin/python2

import sys
import os
import subprocess
import datetime
import pprint as pp
import vcf
import glob
import json

"""create a script that checks which CP samples contain BRCA1 variants in the BRCA1_SGE excel spreadsheet"""


def create_BRCA1_SGE_dict(csv_file, output=None):

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

	#pp.pprint(BRCA1_SGE_dict)
	if output:
		with open(output, "w+") as output_file:
			json.dump(BRCA1_SGE_dict, output_file)
	return BRCA1_SGE_dict

def get_list_of_C01_vcfs(runfolder):

	vcf_list = []

	filenames = glob.glob("%s*/vcfs/C01*.vcf" % runfolder)
	for filename in filenames:
		#print filename
		if ".vcf" in filename and "C01" in filename and filename not in vcf_list:
			vcf_list.append(filename)
		else:
			continue
				
	#print (sorted(vcf_list))
	return sorted(vcf_list)

def compare_BRCA1_SGE_to_vcf(vcf_list, BRCA1_SGE_dict, output=None):

	"""check which CP file contain variants from the supplementary table 1.csv file"""

	C01_sample_SGE_dict = {}

	for vcf_file in vcf_list:

		# check each vcf record against the BRCA1_SGE_dict

		vcf_records = []
	
		vcf_data = vcf.Reader(open(vcf_file, "rb"))
		
		if vcf_file not in C01_sample_SGE_dict.keys():
			C01_sample_SGE_dict[vcf_file] = vcf_records

		print ("analysing {}".format(vcf_file))

		for record in vcf_data:
			if not str(record.CHROM) in BRCA1_SGE_dict.keys():
				continue
			elif not str(record.POS) in BRCA1_SGE_dict[str(record.CHROM)].keys():
				continue
			elif not str(record.REF) in BRCA1_SGE_dict[str(record.CHROM)][str(record.POS)].keys():
				continue
			else:
				if str(record.ALT[0]) in BRCA1_SGE_dict[str(record.CHROM)][str(record.POS)][str(record.REF)]:
					vcf_records.append(str(record))
					C01_sample_SGE_dict[vcf_file] = vcf_records
					print "%s contain the following BRCA1_SGE_variants:" % vcf_file
					print str(record)
				else:
					continue

		#C01_sample_SGE_dict[vcf_file] = vcf_records
	
	if output:
		with open (output, "w+") as output_file:
			for sample in C01_sample_SGE_dict.keys():
				output_file.write("%s" % sample)

	#print C01_sample_SGE_dict.keys()
	return output_file

def main(runfolder, csv_file, SGE_dict_file, C01sample_file):

	# check that the annotated_vcf and excel_table exist as a data source
	assert os.path.exists(csv_file), "{} DO NOT exists".format(csv_file)

	BRCA1_SGE_dict = create_BRCA1_SGE_dict(csv_file, SGE_dict_file)
	vcf_list = get_list_of_C01_vcfs(runfolder)
	compare_BRCA1_SGE_to_vcf(vcf_list, BRCA1_SGE_dict, C01sample_file)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

	"""sys.argv1 = /path/to/CP runfolder, sys.argv2 = supplementary table 1.csv, sys.argv[3] = BRCA1_SGE_dict_file, sys.argv[4] = C01sample_file"""
