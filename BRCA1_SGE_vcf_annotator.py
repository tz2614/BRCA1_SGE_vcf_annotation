#!usr/bin/python2

from __future__ import print_function
import sys
import os
import vcf
import pprint as pp

"""create a vcf file containing the BRCA1 SGE class for 3893 BRCA1 SNVs in the https://www.nature.com/articles/s41586-018-0461-z#Sec9"""

def create_ref_tsv(excel_table):

	"""Create a reference .tsv file containing gene_name, transcript_id, chromosome number, genomic coordinate, ref allele, alt allele and SGE_class for each variant in the excel table"""

	gene_name = int()
	chrom_num = int()
	pos = int()
	ref = int()
	alt = int()
	transcript_id = int()
	SGE_class = int()
	ref_path = os.path.dirname(excel_table)
	tsv_path = os.path.join(ref_path, "BRCA1_SGE_ref.tsv")

	with open(excel_table, "r") as SGE_table:

		for index, line in enumerate(SGE_table):

			fields = line.strip().split(",")

			if line.startswith("gene"):
			
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

			elif "Variant ID" in line:
				continue

			elif "" in line:
				continue
						
			else:
				with open(tsv_path, "a") as ref_tsv:
					ref_tsv.writelines(fields[gene_name] + "\t" + fields[transcript_id] +"\t" +  fields[chrom_num] + "\t" + fields[pos] + "\t" + fields[ref] + "\t" + fields[alt] + "\t" + fields[SGE_class] + "\n")
			
	return tsv_path

def pyvcf_objects_matching_json_ref(annotated_vcf, tsv_path):
	
	"""Query the annotated vcf file and return the pyvcf objects where a variant matches that which is found in the BRCA1_SGE_variants.tsv file."""

	# iterate through the variants in the annotated vcf and return the pyvcf objects that matches the profile of variants in the BRCA1_SGE_variants.tsv file

	rf_path = os.path.dirname(annotated_vcf)
	error_log = os.path.join(rf_path, "vcf_parse_error_log.txt")	

	# open annotated vcf using vcr.Reader module

	with open(annotated_vcf, "r") as vcf_file:

		try: 
			vcf_reader = vcf.Reader(vcf_file)

		except (TypeError, ValueError, NameError, RuntimeError) as e:
			with open(error_log, "a") as err_log:
				err_log.writelines(e)

	# identifying the index within the INFO meta_information field "CSQ" that indicate the transcript_id and gene_name

	with open(annotated_vcf, "r") as vcf_file:
		
		for line in vcf_file:
			
			if line.startswith("##INFO=<ID=CSQ"):
				fields = line.split(":")[-1].split("|")

				for index, field in enumerate(fields):

					if field == "RefSeq":
						transcript_i = index
					if field == "HGNC":
						gene_i = index

	variant_dict = {}
	SGE_classes = []
	BRCA1_variants_in_vcf = []

	# parse the tsv file into a nested dictionary

	with open(tsv_path) as ref_file:

		for line in ref_file:
			
			fields = line.strip().split("\t")
			print (fields)

			if fields[-1] not in SGE_classes:
				SGE_classes.append(fields[-1])

			if fields[0] not in variant_dict.keys():
				variant_dict[fields[0]] = {}

			elif fields[0] in variant_dict.keys() and fields[1] not in variant_dict[fields[0]].keys():
				variant_dict[fields[0]][fields[1]] = {}

			elif fields[0] in variant_dict.keys() and fields[1] in variant_dict[fields[0]].keys() and fields[2] not in variant_dict[fields[0]][fields[1]].keys():
				variant_dict[fields[0]][fields[1]][fields[2]] = {}

			elif fields[0] in variant_dict.keys() and fields[1] in variant_dict[fields[0]].keys() and fields[2] in variant_dict[fields[0]][fields[1]].keys() and fields[3] not in variant_dict[fields[0]][fields[1]][fields[2]].keys():
				variant_dict[fields[0]][fields[1]][fields[2]][fields[3]] = {}

			elif fields[0] in variant_dict.keys() and fields[1] in variant_dict[fields[0]].keys() and fields[2] in variant_dict[fields[0]][fields[1]].keys() and fields[3] in variant_dict[fields[0]][fields[1]][fields[2]].keys() and fields[4] not in variant_dict[fields[0]][fields[1]][fields[2]][fields[3]].keys():
				variant_dict[fields[0]][fields[1]][fields[2]][fields[3]][fields[4]] = {}

			elif fields[0] in variant_dict.keys() and fields[1] in variant_dict[fields[0]].keys() and fields[2] in variant_dict[fields[0]][fields[1]].keys() and fields[3] in variant_dict[fields[0]][fields[1]][fields[2]].keys() and fields[4] in variant_dict[fields[0]][fields[1]][fields[2]][fields[3]].keys() and fields[5] not in variant_dict[fields[0]][fields[1]][fields[2]][fields[3]][fields[4]].keys():
				variant_dict[fields[0]][fields[1]][fields[2]][fields[3]][fields[4]][fields[5]] = fields[6]
			else:
				continue

	pp.pprint(variant_dict)

	for record in vcf_reader:

		for line in record.INFO["CSQ"]:
			gene_name = line.split("|")[gene_i]
			transcript_id = line.split("|")[transcript_i]

			if variant_dict[gene_name][transcript_id][record.CHROM][record.POS][record.REF][record.ALT] in SGE_classes and record not in BRCA1_variants_in_vcf:
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
	tsv_path = create_ref_tsv(excel_table)
	
	record_list = pyvcf_objects_matching_json_ref(annotated_vcf, tsv_path)

if __name__ == "__main__":
	main(sys.argv[1], sys.argv[2])