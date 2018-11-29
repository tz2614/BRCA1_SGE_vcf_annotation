#!usr/bin/python3

import sys
import os
import subprocess
import datetime
import pprint as pp

"""create a script that checks which CP samples contain BRCA1 variants in the BRCA1_SGE excel spreadsheet"""


def create_BRCA1_SGE_ref_txt(csv_file):

	#create a tab-delimited file - BRCA1_SGE_ref.txt containing the CHROM, POS, REF and ALT

	dir_path = os.path.dirname(csv_file)
	BRCA1_SGE_ref = os.path.join(dir_path, "BRCA1_SGE_ref.txt")

	if os.path.exists(BRCA1_SGE_ref):
		return BRCA1_SGE_ref
	
	with open (csv_file, "r") as csv_file:

		print ("generating BRCA1_SGE_ref.txt")

		lines = [line.split(',') for line in csv_file]

		for index, fields in enumerate(lines):

			if index in [0, 1]:
				continue

		# sort excel_table according to genomic position

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
					
					# write each variant information for the BRCA1_SGE_ref.txt file
					with open(BRCA1_SGE_ref, "a") as ref_file:
						ref_file.writelines(CHROM + "\t" + POS + "\t" + REF + "\t" + ALT + "\n")				
		
			else:
				break
	
	print ("BRCA1_SGE_ref.txt generated")
	return BRCA1_SGE_ref


def compare_BRCA1_SGE_2_CP(BRCA1_CP_txt, BRCA1_SGE_ref):

	"""check which CP file contain variants from the supplementary table 1.csv file"""

	with open (BRCA1_CP_txt, "r") as CP_file:

		#print ("comparing {} against {}".format(BRCA1_CP_txt, BRCA1_SGE_ref))
		
		CP_data = [line.strip().split("\t") for line in CP_file]		
				
	with open (BRCA1_SGE_ref) as SGE_file:

		SGE_data = [line.strip().split('\t') for line in SGE_file]

	for index, CP_fields in enumerate(CP_data):

		if index == 0:
			continue
		else:
			name = CP_fields[0]
			pos = CP_fields[1]
			ref = CP_fields[2]
			alt = CP_fields[3]

		for SGE_fields in SGE_data:

			if pos == SGE_fields[1] and ref == SGE_fields[2] and alt == SGE_fields[3]:

				print ("Sample", name, " contain BRCA1_SGE_variant", "\t", "chr:17", "\t", pos, "\t", ref, "\t", alt, "\n")		
			else:
				continue
	return 

def main(BRCA1_CP_txt, csv_file):

	# check that the annotated_vcf and excel_table exist as a data source
	assert os.path.exists(csv_file), "{} DO NOT exists".format(csv_file)
	assert os.path.exists(BRCA1_CP_txt), "{} DO NOT exists".format(BRCA1_CP_txt)

	BRCA1_SGE_ref = create_BRCA1_SGE_ref_txt(csv_file)
	compare_BRCA1_SGE_2_CP(BRCA1_CP_txt, BRCA1_SGE_ref)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])

	"""sys.argv1 = BRCA1_CP.txt, sys.argv2 = supplementary table 1.csv"""
