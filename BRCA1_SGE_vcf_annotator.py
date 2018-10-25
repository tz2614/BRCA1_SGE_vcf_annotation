#!usr/bin/python2

from __future__ import print_function
import sys
import os
import vcf
import subprocess
import datetime
import pprint as pp

"""create a vcf file containing the BRCA1 SGE class for 3893 BRCA1 SNVs in the https://www.nature.com/articles/s41586-018-0461-z#Sec9"""

def create_ref_vcf(excel_table):

	"""Create a BRCA1 reference .vcf file containing the BRCA1 SGE score (posterior probability) and class 
	of 3893 BRCA1 SNVs in the https://www.nature.com/articles/s41586-018-0461-z#Sec9"""

	excel_table = os.path.abspath(excel_table)
	ref_path = os.path.dirname(excel_table)
	vcf_path = os.path.join(ref_path, "BRCA1_SGE_ref.vcf")

	filedate = str(datetime.date.today())
	filedate = "".join(filedate.split("-"))

	# write header for the BRCA1_SGE_ref.vcf

	with open(vcf_path, "w") as vcf_file:
		vcf_file.writelines('##fileformat=VCFv4.0\n')
		vcf_file.writelines('##fileDate={}\n'.format(filedate))
		vcf_file.writelines('##reference=hg19\n')
		vcf_file.writelines('##INFO=<ID=BRCA1_SGE,Number=.,Type=String,Description="BRCA1_SGE_score and BRCA1_SGE_class. Format: allele|score|class">\n')
		vcf_file.writelines('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
		print ('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO')

	assert os.path.exists(vcf_path), "vcf reference file DO NOT exists"

	# equate the data columns from excelsheet to the 8 fixed fields in vcf, CHROM=chromosome number, POS=position(hg19), ID=".", REF=ref, ALT=alt, QUAL=".", FILTER= ".", INFO= BRCA1=allele|score|class

	with open(excel_table, "r") as SGE_table:

		lines = [line.split(',') for line in SGE_table]

		for index, fields in enumerate(lines):

			if index in [0, 1]:
				continue

		# sort excel_table according to genomic position

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
					elif field == "p.nonfunctional":
						SGE_score = column_index
					elif field == "func.class":
						SGE_class = column_index
					else:
						continue

				CHROM = fields[chrom_num]
				POS = fields[pos]
				ID = "."
				REF = fields[ref]
				ALT = fields[alt]
				QUAL = "."
				FILTER = "."
				score = fields[SGE_score]
				func_class = fields[SGE_class]
				BRCA1_SGE = "BRCA1={}|{}|{}".format(ALT, score, func_class)

				# write each variant information for the BRCA1_SGE_ref.vcf file
				with open(vcf_path, "a") as vcf_file:
					vcf_file.writelines(CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + BRCA1_SGE + "\n")				
					print(CHROM + "\t" + POS + "\t" + ID + "\t" + REF + "\t" + ALT + "\t" + QUAL + "\t" + FILTER + "\t" + BRCA1_SGE + "\n")
			else:
				break
	
	return vcf_path

def create_vcf_header(dir_path, vcf_path):

	vcf_header = os.path.join(dir_path, ".hdr")

	with open (vcf_path, "r") as vcf_file:
		for line in vcf_file:
			if line.startswith("##INFO=<ID=BRCA1_SGE_class"):
				vcf_header = line
				with open (vcf_header, "w+") as hdr:
					hdr.writelines(line)
			else:
				break

	return vcf_header

def sort_bgzip_index(vcf_path):

	vcf_path = os.path.abspath(vcf_path)
	assert os.path.exists(vcf_path), "vcf file DO NOT exist"

	dir_path = os.path.dirname(vcf_path)
	vcf_file = vcf_path.split("/")[-1]

	bigzp_path = os.path.join(vcf_path, ".gz")
	bigzp_file = bigzp_path.split("/")[-1]
	tabix_path = os.path.join(bigzp_path, ".tbi")
	tabix_file = tabix_path.split("/")[-1]
	
	files = os.listdir(dir_path)

	if not bigzp_file in files:
		print ("create {}".format(bigzp_file))
		create_bgzip = "cd {}; bgzip {}".format(dir_path, vcf_file)
		subprocess.call(create_bgzip, shell=True)

	if not tabix_file in files:
		print ("create {}".format(tabix_file, shell=True))
		create_tabix = "cd {}; tabix -p vcf {}".format(dir_path, vcf_file)
		subprocess.call(create_tabix, shell=True)

	return dir_path, bigzp_path, tabix_path

def annotate_vcf(dir_path, bgzip_file, vcf_header, BRCA1_vcf,  vcf_file):
	
	vcf_dot = vcf_file.split("/")[-1].split(".")
	for index, field in enumerate(vcf_dot):
		if field == "annotated":
			vcf_dot[index] = "BRCA1_annotated"

	annotated_vcf = ".".join(vcf_dot)
	print (annotated_vcf)

	output_file = os.path.join(dir_path, annotated_vcf)
	command = "cd {}; bcftools -a {} -h {} -c CHROM,POS,REF,ALT,INFO in.vcf.gz -o {}".format(dir_path, bgzip_file, vcf_header, output_file)
	subprocess.call(command, shell=True)

def main(excel_table, vcf_file):

	# check that the annotated_vcf and excel_table exist as a data source
	assert os.path.exists(vcf_file), "{} DO NOT exists".format(vcf_file)
	assert os.path.exists(excel_table), "{} DO NOT exists".format(excel_table)

	# generate an artifical vcf from a tab delimited text file
	BRCA1_vcf = create_ref_vcf(excel_table)
	#dir_path, bgzip_file, tabix_path = sort_bgzip_index(BRCA1_vcf)
	#vcf_header = create_vcf_header(dir_path, vcf_file)
	#annotate_vcf(dir_path, bgzip_file, vcf_header, BRCA1_vcf, vcf_file)
	
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])