#!usr/bin/python3

import sys
import os
import vcf
import subprocess
import datetime
import pprint as pp

def remove_INFO(vcf):
	
	dir_path = os.path.dirname(vcf)
	new_vcf = os.path.join(dir_path, "BRCA1_test_vcf.vcf")
	
	with open(vcf, "r") as vcf_file:
		
		for line in vcf_file:
			
			if line.startswith("#"):
				with open (new_vcf, "a") as new_vcf_file:
					new_vcf_file.writelines(line)

			elif line.startswith("17"):
				fields = line.split("\t")
				line = line.replace(fields[-1], "") + "\n"
				with open (new_vcf, "a") as new_vcf_file:
					new_vcf_file.writelines(line)
			
			else:
				continue

def main(vcf):

	remove_INFO(vcf)

if __name__ == "__main__":
	main(sys.argv[1])