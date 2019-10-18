#!usr/bin/python3

import sys
import os
import subprocess
import datetime
import BRCA1_SGE_ref

def create_vcf_header(ref_vcf):

	"""create a header file containing the reference vcf header"""

	header_file = ref_vcf + ".hdr"

	with open (ref_vcf, "r") as vcf_file:
		for line in vcf_file:
			if line.startswith("##INFO=<ID=BRCA1_SGE"):
				#print ("finding header in vcf")
				with open (header_file, "w") as hdr:
					hdr.writelines(line)
			else:
				continue
	print ("header file generated")
	return header_file


def sort_bgzip_index(ref_vcf):

	"""sort the lines in the vcf according to CHROM and POS (numerically), compress the file into bgzip file format, and then index the .gz file"""

	ref_vcf = os.path.abspath(ref_vcf)
	print (ref_vcf)
	assert os.path.exists(ref_vcf), "reference vcf file DO NOT exist"

	bgzip_file = ref_vcf + ".gz"
	tabix_file = bgzip_file + ".tbi"

	if not os.path.exists(bgzip_file):
		print ("bgzip file DO NOT exist")
		print ("change to current directory and sort {} according to CHROM and POS, and create {}".format(ref_vcf, bgzip_file))
		sort_bgzip = "grep -v '#' {} | sort -k1V -k2n -t$'\t' | bgzip -c > {};".format(ref_vcf, bgzip_file)
		subprocess.call(sort_bgzip, shell=True)
	else:
		print("bgzip file already exists")

	if not os.path.exists(tabix_file):
		print("tabix file DO NOT exist")
		print ("create {}".format(tabix_file))
		create_tabix = "tabix -p vcf {}".format(bgzip_file)
		subprocess.call(create_tabix, shell=True)
	else:
		print("tabix file already exists")

	return bgzip_file, tabix_file

def annotate_vcf(bgzip_file, vcf_header, ref_vcf, vcf_file):

	"""using bcftools and add BRCA1_SGE annotation to the vcf_file, using bgzip_file, vcf_header, BRCA1_vcf and vcf_file as arguments"""

	# check that the ref_vcf exist as a data source
	assert ".vcf" in vcf_file, "file DO NOT end with .vcf"
	assert os.path.exists(ref_vcf), "{} DO NOT exists".format(ref_vcf)

	if ".annotated." in vcf_file:
		annotated_vcf_file = vcf_file.replace(".annotated.", ".BRCA1_annotated.")
	elif "BRCA1_annotated" in vcf_file:
		print("vcf file already annotated")
		return annotated_vcf_file
	else:
		annotated_vcf_file = vcf_file.replace(".vcf", ".BRCA1_annotated.vcf")

	print ("create new vcf with additional BRCA1 annotation")
	command = "bcftools annotate -a {} -h {} -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,'+INFO/BRCA1_SGE' -o {} {}".format(bgzip_file, vcf_header, annotated_vcf_file, vcf_file)
	subprocess.call(command, shell=True)

	return annotated_vcf_file

def main(ref_vcf, vcf_file):

	# use BRCA1_SGE_ref.vcf as reference file containing all the BRCA1_SGE variants to annotate new vcf

	ref_vcf = os.path.abspath(ref_vcf)
	vcf_file = os.path.abspath(vcf_file)
	vcf_header = create_vcf_header(ref_vcf)
	bgzip_file, tabix_file = sort_bgzip_index(ref_vcf)
	annotate_vcf( bgzip_file, vcf_header, ref_vcf, vcf_file)
	
if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])

# navigate to where the script is located from command line, and enter "python3", 

# the name of the script "BRCA1_SGE_vcf_annotator.py"

# sys.argv[1] - /path/to/BRCA1_SGE_ref.vcf
# sys.argv[2] - /path/to/vcf you wish to annotate