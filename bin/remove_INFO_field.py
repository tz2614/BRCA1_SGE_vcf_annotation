#!usr/bin/python3

import sys
import os

"""remove INFO field in vcf files"""

def remove_INFO(vcf):
	
	dir_path = os.path.dirname(vcf)
	vcf_filename = os.path.basename(vcf)
	new_vcf_filename = vcf_filename.replace(".vcf", ".INFO_removed.vcf")

	new_vcf = os.path.join(dir_path, new_vcf_filename)
	
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

# navigate to where the script is located from command line, and enter "python3", then the name of the script "remove_INFO_field.py"

# followed by the path of the vcf containing INFO fields - sys.argv[1]
