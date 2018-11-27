#!usr/bin/python3

import sys
import os
import subprocess
import datetime
import BRCA1_SGE_ref

def create_vcf_header(vcf_path):

	"""create a header file containing the reference vcf header"""

	header_file = vcf_path + ".hdr"

	with open (vcf_path, "r") as vcf_file:
		for line in vcf_file:
			if line.startswith("##INFO=<ID=BRCA1_SGE"):
				print ("finding header in vcf")
				with open (header_file, "w") as hdr:
					hdr.writelines(line)
			else:
				continue
	print ("header file generated")
	return header_file


def sort_bgzip_index(ref_vcf):

	"""sort the Data lines in the vcf according to CHROM and POS (numerically), compress the file into bgzip file format, and then index the .gz file"""

	ref_vcf = os.path.abspath(ref_vcf)
	print (ref_vcf)
	assert os.path.exists(ref_vcf), "vcf file DO NOT exist"

	dir_path = os.path.dirname(ref_vcf)
	vcf_file = ref_vcf.split("/")[-1]
	print (vcf_file)

	bgzip_file = vcf_file + ".gz"
	tabix_file = bgzip_file + ".tbi"

	print ("change to current directory and sort {} according to CHROM and POS, and create {}".format(vcf_file, bgzip_file))
	sort_bgzip = "cd {}; grep -v '#' {} | sort -k1,1n -k2,2n -t$'\t' | bgzip -c > {};".format(dir_path, vcf_file, bgzip_file)
	subprocess.call(sort_bgzip, shell=True)

	print ("create {}".format(tabix_file))
	create_tabix = "cd {}; tabix -p vcf {}".format(dir_path, bgzip_file)
	subprocess.call(create_tabix, shell=True)

	return dir_path, bgzip_file, tabix_file


def annotate_vcf(dir_path, bgzip_file, vcf_header, BRCA1_ref_vcf, vcf_file):

	"""using bcftools and add BRCA1_SGE annotation to the vcf_file, using dir_path, bgzip_file, vcf_header, BRCA1_vcf and vcf_file as arguments"""

	vcf_file = str(vcf_file)

	assert ".vcf" in vcf_file, "file DO NOT end with .vcf"

	if ".annotated." in vcf_file:
		annotated_vcf_file = vcf_file.replace(".annotated.", ".BRCA1_annotated.")
	else:
		annotated_vcf_file = vcf_file.replace(".vcf", ".BRCA1_annotated.vcf")

	print ("create new vcf with additional BRCA1 annotation")
	command = "cd {}; bcftools annotate -a {} -h {} -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,'+INFO/BRCA1_SGE' -o {} {}".format(dir_path, bgzip_file, vcf_header, annotated_vcf_file, vcf_file)
	subprocess.call(command, shell=True)


def main(excel_table, vcf_file):

	# check that the annotated_vcf and excel_table exist as a data source
	assert os.path.exists(vcf_file), "{} DO NOT exists".format(vcf_file)
	assert os.path.exists(excel_table), "{} DO NOT exists".format(excel_table)

	# generate an artifical vcf from a tab delimited text file
	BRCA1_ref_vcf = BRCA1_SGE_ref.create_ref_vcf(excel_table)
	vcf_header = create_vcf_header(BRCA1_ref_vcf)
	dir_path, bgzip_file, tabix_file = sort_bgzip_index(BRCA1_ref_vcf)
	annotate_vcf(dir_path, bgzip_file, vcf_header, BRCA1_ref_vcf, vcf_file)
	
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])

# navigate to where the script is located from command line, and enter "python3", then the name of the script "BRCA1_SGE_vcf_annotator.py"

# followed by the absolute file path of the .csv table containing the BRCA1_SGE annotations - sys.argv[1]

# and the absolute file path of the vcf you wish to annotate - sys.argv[2]
