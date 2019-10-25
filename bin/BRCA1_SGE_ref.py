#!usr/bin/python3

import sys
import os
import subprocess
import datetime

"""create a vcf file containing the BRCA1 SGE class for 3893 BRCA1 SNVs in the https://www.nature.com/articles/s41586-018-0461-z#Sec9"""


def create_ref_vcf(excel_table):

	"""Create a BRCA1 reference .vcf file containing the function score (log2 scaled) and class 
	of 3893 BRCA1 SNVs in the https://www.nature.com/articles/s41586-018-0461-z#Sec9"""

	timestamp = "_".join(str(datetime.datetime.now()).split(" "))
	excel_table = os.path.abspath(excel_table)
	supplementary_table_path = os.path.dirname(excel_table)
	
	filedate = str(datetime.date.today())
	filedate = "".join(filedate.split("-"))

	ref_vcf_path = os.path.join(supplementary_table_path, "BRCA1_SGE_ref_{}.vcf".format(filedate))

	# write header for the BRCA1_SGE_ref.vcf

	with open(ref_vcf_path, "w") as vcf_file:
		vcf_file.writelines('##fileformat=VCFv4.1\n')
		vcf_file.writelines('##fileDate={}\n'.format(filedate))
		vcf_file.writelines('##reference=hg19\n')
		vcf_file.writelines('##contig=<ID=17,length=81195210>\n')
		vcf_file.writelines('##INFO=<ID=BRCA1_SGE,Number=A,Type=String,Description="BRCA1_SGE_score and BRCA1_SGE_class. Format: allele|score|class">\n')
		vcf_file.writelines('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
		
	# equate the data columns from excelsheet to the 8 fixed fields in vcf, CHROM=chromosome number, POS=position(hg19),
	# ID=".", REF=ref, ALT=alt, QUAL=".", FILTER= ".", INFO= BRCA1=allele|score|class

	with open(excel_table, "r") as SGE_table:

		print ("generating {}".format(ref_vcf_path))

		lines = [line.split(',') for line in SGE_table]

		for index, fields in enumerate(lines):

			if index in [0, 1]:
				continue

			elif index >= 2:
				
				for column_index, field in enumerate(fields):

					if field == "chromosome":
						chrom_num = column_index
					elif field == "position (hg19)":
						pos = column_index
					elif field == "reference":
						ref = column_index
					elif field == "alt":
						alt = column_index
					elif field == "function.score.mean":
						SGE_score = column_index
					elif field == "func.class":
						SGE_class = column_index
					else:
						continue

				if index > 2:
					CHROM = fields[chrom_num]
					POS = fields[pos]
					ID = "."
					REF = fields[ref]
					ALT = fields[alt]
					QUAL = "."
					FILTER = "."
					score = float(fields[SGE_score])
					func_class = fields[SGE_class]
					BRCA1_SGE = "{}|{}|{}".format(ALT, score, func_class)

					# write each variant information for the BRCA1_SGE_ref.vcf file
					with open(ref_vcf_path, "a") as vcf_file:
						vcf_file.writelines(CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + BRCA1_SGE + "\n")				
		
			else:
				break
	
	print ("{} generated".format(ref_vcf_path))
	return ref_vcf_path


def main(excel_table):

	# check that the excel_table exist as a data source
	assert os.path.exists(excel_table), "{} DO NOT exists".format(excel_table)

	# generate an artifical vcf from a tab delimited text file
	BRCA1_ref_vcf = create_ref_vcf(excel_table)
	
if __name__ == "__main__":
	main(sys.argv[1])

# navigate to where the script is located from command line, and enter

# $python3 BRCA1_SGE_ref.py 

# followed by the absolute file path of the .csv table containing the BRCA1_SGE annotations - sys.argv[1]