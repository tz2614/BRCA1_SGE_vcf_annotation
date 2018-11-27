#!usr/bin/python2

import sys
import os

"""check if genomic ref match transcript ref in excel table"""

def compare_genomic_transcript_variants(excel_table):

	dna_table = {"A":"T", "T":"A", "C":"G", "G":"C"}

	with open(excel_table, "r") as ref_table:
		for index, line in enumerate(ref_table):
			if index in [0, 1, 2]:
				pass
			else:
				fields = line.strip().split(",")
				if fields[3] in dna_table.keys():
					assert fields[7] == dna_table[fields[3]], "position {} contain mismatch refs between genomic and transcript coordinates"
				else:
					print ("ref not present in dna_table")

	return	ref_table			

def main(excel_table):
	assert os.path.exists(excel_table), "{} DO NOT exists".format(excel_table)

	compare_genomic_transcript_variants(excel_table)

if __name__ == "__main__":
	main(sys.argv[1])

	# provide the excel table as argument when executing the script


		
