#!usr/bin/python2

from __future__ import print_function
import sys
import os
import vcf
import json
import pprint as pp

"""create a vcf file containing the BRCA1 SGE class for 3893 BRCA1 SNVs in the https://www.nature.com/articles/s41586-018-0461-z#Sec9"""

def create_ref_json(excel_table):

	"""Create a JSON file containing gene_name, genomic coordinate, transcript_id and ref and alt allele for each variant in the excel table"""

	BRCA1_SGE_variants = []
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
				BRCA1_SGE_variants.append({
					"gene": fields[gene_name],
					"chromosome": fields[chrom_num],
					"position (hg19)": fields[pos],
					"ref": fields[ref],
					"alt": fields[alt],
					"transcript_id": fields[transcript_id],
					"SGE_class": fields[SGE_class]
					})

	#pp.pprint(BRCA1_SGE_variants)
	return json.dumps(BRCA1_SGE_variants, sort_keys=True, indent=4)

def pyvcf_objects_matching_json_ref(annotated_vcf, JSON_file):
	
	"""Query the annotated vcf file and return the pyvcf objects where a variant matches that which is found in the JSON file."""

	# iterate through the variants in the annotated vcf and return the variants that matches the profile of variants in the JSON_file

	rf_path = os.path.dirname(annotated_vcf)
	error_log = os.path.join(rf_path, "vcf_parse_error_log.txt")	
	
	BRCA1_variants_in_vcf = []

	with open(annotated_vcf, "r") as vcf_file:

		try: 
			vcf_reader = vcf.Reader(vcf_file)

		except (TypeError, RuntimeError, NameError, ValueError) as e:
			with open(error_log, "a") as err_log:
				err_log.writelines(e)

		with open(annotated_vcf, "r") as vcf_file:
			for line in vcf_file:
				if line.startswith("##INFO=<ID=CSQ"):
					fields = line.split(":")[-1].split("|")
					for index, field in enumerate(fields):
						if field == "RefSeq":
							i = index

		with open(JSON_file, "r") as ref_file:
			data = json.load(ref_file)

		for record in vcf_reader:
			for variant in data:
				if record.CHROM == variant["chromosome"] and record.POS == variant["position (hg19)"] and record.REF == variant["ref"] and record.ALT == variant["alt"] and [l.split("|")[i] for l in record.INFO["CSQ"]] == variant["transcript_id"] and record not in BRCA1_variants_in_vcf:
					print (record.CHROM, record.POS, record.REF, record.ALT)
					BRCA1_variants_in_vcf.append(record)

				else:
					continue
	
	return  (BRCA1_variants_in_vcf)

def main(annotated_vcf, excel_table):

	# check that the annotated_vcf and excel_table exist as a data source
	assert os.path.exists(annotated_vcf), "{} DO NOT exists".format(annotated_vcf)
	assert os.path.exists(excel_table), "{} DO NOT exist".format(excel_table)

	# create a list of vcf variants in the form of pyvcf objects
	JSON_file = create_ref_json(excel_table)
	#print (JSON_file)
	record_list = pyvcf_objects_matching_json_ref(annotated_vcf, JSON_file)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])