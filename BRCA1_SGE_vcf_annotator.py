#!usr/bin/python2

from __future__ import print_function
import sys
import os
import vcf
import json
import pprint as pp

"""create a vcf file containing the BRCA1 SGE class for 3893 BRCA1 SNVs in the https://www.nature.com/articles/s41586-018-0461-z#Sec9"""

def create_vcf(excel_table):

	BRCA1_SGE_fields = {}
	gene_name = int()
	chrom_num = int()
	pos = int()
	ref = int()
	alt = int()
	transcript_id = int()
	SGE_class = int()

	with open(excel_table, "r") as SGE_table:
		for index, line in enumerate(SGE_table):
			fields = line.split(",")
			for column_index, field in enumerate(fields):
				if field == "gene":
					gene_name = column_index
				elif field == "chromosome":
					chrom_num = column_index
				elif field == "position (hg19)":
					pos = column_index
				elif field == "reference":
					ref = column_index
				elif field == "alt":
					alt = column_index
				elif field == "transcript_ID":
					transcript_id = column_index
				elif field == "func.class":
					SGE_class = column_index
				else:
					continue
			else:
				BRCA1_SGE_fields[str(index)] = []
				BRCA1_SGE_fields[str(index)].append({
					"gene": fields[gene_name],
					"chromosome": fields[chrom_num],
					"position (hg19)": fields[pos],
					"ref": fields[ref],
					"alt": fields[alt],
					"transcript_id": fields[transcript_id],
					"SGE_class": fields[SGE_class]
					})

	pp.pprint(BRCA1_SGE_fields)
	return json.dumps(BRCA1_SGE_fields, sort_keys=True, indent=4)

def main(annotated_vcf, excel_table):

	# check that the annotated_vcf and excel_table exist as a data source
	assert os.path.exists(annotated_vcf), "{} DO NOT exists".format(annotated_vcf)
	assert os.path.exists(excel_table), "{} DO NOT exist".format(excel_table)

	# create a list of vcf variants in the form of pyvcf objects
	JSON_file = create_vcf(excel_table)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])