#!usr/bin/python3

import BRCA1_SGE_vcf_annotator
import pytest
import os
import shutil

# VCF variant records that should be present in the annotated vcf for pos_variants.vcf
pos_var1_annot = "ABHet=0.52;AC=1;AF=0.5;AN=2;BaseQRankSum=10.362;DP=245;Dels=0;FS=2.017;HRun=0;HaplotypeScore=17.5355;MQ=55.8;MQ0=1;MQRankSum=0.294;OND=0.01;QD=10.34;ReadPosRankSum=2.074;SB=-708.03;BRCA1_SGE=T|-0.68844702|FUNC"
pos_var2_annot = "ABHet=0.51;AC=1;AF=0.5;AN=2;BaseQRankSum=-8.311;DP=435;Dels=0;FS=6.464;HRun=0;HaplotypeScore=21.1309;MQ=57.79;MQ0=1;MQRankSum=-1.237;QD=8.48;ReadPosRankSum=-2.61;SB=-1201.33;BRCA1_SGE=C|0.408757217|FUNC"
pos_var3_annot = "ABHet=0.51;AC=1;AF=0.5;AN=2;BaseQRankSum=9.959;DP=940;Dels=0;FS=14.99;HRun=0;HaplotypeScore=43.1044;MQ=57.16;MQ0=1;MQRankSum=6.032;QD=9.4;ReadPosRankSum=-1.558;SB=-4820.63;BRCA1_SGE=T|0.143488826|FUNC"
pos_var4_annot = "ABHet=0.48;AC=1;AF=0.5;AN=2;BaseQRankSum=12.957;DP=999;DS;Dels=0;FS=0.769;HRun=0;HaplotypeScore=66.4204;MQ=58.1;MQ0=1;MQRankSum=1.872;OND=0;QD=10.26;ReadPosRankSum=1.395;SB=-5821.58;BRCA1_SGE=A|0.203667345|FUNC"

# VCF variant records that should be present in the annotated vcf for neg_variants.vcf
neg1_var_annot = "ABHet=0.83;AC=1;AF=0.5;AN=2;BaseQRankSum=4.709;DP=337;Dels=0;FS=17.206;HRun=0;HaplotypeScore=50.2036;MQ=51.85;MQ0=0;MQRankSum=-10.412;OND=0;QD=0.05;ReadPosRankSum=0.346;SB=-29.84"
neg2_var_annot = "ABHet=0.87;AC=1;AF=0.5;AN=2;BaseQRankSum=-1.144;DP=664;Dels=0;FS=54.801;HRun=0;HaplotypeScore=366.431;MQ=45.89;MQ0=0;MQRankSum=-4.984;OND=0.34;QD=0.1;ReadPosRankSum=-2.911;SB=-0.01"
neg3_var_annot = "ABHom=1;AC=2;AF=1;AN=2;BaseQRankSum=1.678;DP=957;DS;Dels=0;FS=0;HRun=0;HaplotypeScore=58.3949;MQ=55.33;MQ0=0;MQRankSum=1.016;OND=0.01;QD=25.49;ReadPosRankSum=0.903;SB=-12358"
neg4_var_annot = "ABHom=1;AC=2;AF=1;AN=2;DP=3;Dels=0;FS=0;HRun=0;HaplotypeScore=2.9664;MQ=25.63;MQ0=0;QD=6.38;SB=-0.01"

# reference header information
ref_header = '##INFO=<ID=BRCA1_SGE,Number=A,Type=String,Description="BRCA1_SGE_score and BRCA1_SGE_class. Format: allele|score|class">'

def test_vcf_header(create_test_work_dir):

	# function to test if header file begins with ##INFO

	test_work_dir, new_ref_hdr, new_ref_vcf, new_ref_vcf_gz, new_ref_vcf_gz_tbi, new_merged_vcf_gz, new_merged_vcf_gz_tbi = create_test_work_dir

	hdr_file = BRCA1_SGE_vcf_annotator.create_vcf_header(new_ref_vcf)

	with open(hdr_file, "r") as hdr:
		for index, line in enumerate(hdr):
			if index == 0:
				assert ref_header == line.strip(), "header do not contain right INFO"
			else:
				break


def test_annotate_ref(create_test_vcfs_dir, create_test_work_dir):

	# function to test if correct functions score and class for BRCA1_SGE postive and negative variants have been added to reference.vcf

	# create the temporary test files in the test work directory

	test_vcf_dir, pos_vcf, neg_vcf = create_test_vcfs_dir
	print (test_vcf_dir,pos_vcf, neg_vcf)

	new_ref_dir, new_ref_hdr, new_ref_vcf, new_ref_vcf_gz, new_ref_vcf_gz_tbi, new_merged_vcf_gz, new_merged_vcf_gz_tbi = create_test_work_dir
	print (new_ref_dir, new_ref_hdr, new_ref_vcf, new_ref_vcf_gz, new_ref_vcf_gz_tbi, new_merged_vcf_gz, new_merged_vcf_gz_tbi)

	# generate the annotated reference vcf
	pos_annotated_vcf = BRCA1_SGE_vcf_annotator.annotate_ref(new_ref_vcf_gz, new_ref_hdr, pos_vcf)

	neg_annotated_vcf = BRCA1_SGE_vcf_annotator.annotate_ref(new_ref_vcf_gz, new_ref_hdr, neg_vcf)


	# assess annotation of BRCA1_SGE positive variants

	with open(pos_annotated_vcf, "r") as pos_a_vcf:
		for i_pos, line in enumerate(pos_a_vcf):
			if line.startswith("#CHROM"):
				var1 = i_pos + 1

				if i_pos == var1:
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41197829" \
					and fields[3] == "G" and fields[4] == "T" \
					and fields[8] == pos_var1_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_pos == (var1 + 1):
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41199671" \
					and fields[3] == "T" and fields[4] == "C" \
					and fields[8] == pos_var2_annot, \
					"#variant INFO wrong  \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_pos == (var1 + 2):
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41222975" \
					and fields[3] == "C" and fields[4] == "T" \
					and fields[8] == pos_var2_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_pos == (var1 + 3):
					fields = line.split("\t")

					assert fields[0] == "17"	and fields[1] == "41267810" \
					and fields[3] == "G"and fields[4] == "A" \
					and fields[8] == pos_var2_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				else:
					continue

	# assess annotation of BRCA1_SGE negative variants

	with open(neg_annotated_vcf, "r") as neg_a_vcf:
		for i_neg, line in enumerate(neg_a_vcf):
			if line.startswith("#CHROM"):
				var1 = i_neg + 1

				if i_neg == var1:
					fields = line.split("\t")


					assert fields[0] == "17" and fields[1] == "41199935" \
					and fields[3] == "C" and fields[4] == "T" \
					and fields[8] == neg_var1_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_neg == (var1 + 1):
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41202400" \
					and fields[3] == "T" and fields[4] == "C" \
					and fields[8] == neg_var2_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_neg == (var1 + 2):
					fields = line.split("\t")

					assert fields[0] == "17"	and fields[1] == "41226675" \
					and fields[3] == "A" and fields[4] == "T" \
					and fields[8] == neg_var3_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_neg == (var1 + 3):
					fields = line.split("\t")

					assert fields[0] == "17"	and fields[1] == "41299122" \
					and fields[3] == "A" and fields[4] == "C" \
					and fields[8] == neg_var4_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])


def test_annotate_vcf(create_test_vcfs_dir, create_test_work_dir):

	# function to test if correct functions score and class for BRCA1_SGE postive and negative variants have been added to <sample>.vcf

	# create the temporary test files in the test work directory

	test_vcf_dir, pos_vcf, neg_vcf = create_test_vcfs_dir
	print (test_vcf_dir,pos_vcf, neg_vcf)
	new_ref_dir, new_ref_hdr, new_ref_vcf, new_ref_vcf_gz, new_ref_vcf_gz_tbi, new_merged_vcf_gz, new_merged_vcf_gz_tbi = create_test_work_dir
	print (new_ref_dir, new_ref_hdr, new_ref_vcf, new_ref_vcf_gz, new_ref_vcf_gz_tbi, new_merged_vcf_gz, new_merged_vcf_gz_tbi)

	# generate the annotated vcf
	pos_annotated_vcf = BRCA1_SGE_vcf_annotator.annotate_vcf(new_ref_vcf_gz, new_ref_hdr, pos_vcf, new_merged_vcf_gz, new_merged_vcf_gz_tbi)

	neg_annotated_vcf = BRCA1_SGE_vcf_annotator.annotate_vcf(new_ref_vcf_gz, new_ref_hdr, neg_vcf, new_merged_vcf_gz, new_merged_vcf_gz_tbi)


	# assess annotation of BRCA1_SGE postive variants

	with open(pos_annotated_vcf, "r") as pos_a_vcf:
		for i_pos, line in enumerate(pos_a_vcf):
			if line.startswith("#CHROM"):
				var1 = i_pos + 1

				if i_pos == var1:
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41197829" \
					and fields[3] == "G" and fields[4] == "T" \
					and fields[8] == pos_var1_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_pos == (var1 + 1):
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41199671" \
					and fields[3] == "T" and fields[4] == "C" \
					and fields[8] == pos_var2_annot, \
					"#variant INFO wrong  \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_pos == (var1 + 2):
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41222975" \
					and fields[3] == "C" and fields[4] == "T" \
					and fields[8] == pos_var3_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_pos == (var1 + 3):
					fields = line.split("\t")

					assert fields[0] == "17"	and fields[1] == "41267810" \
					and fields[3] == "G"and fields[4] == "A" \
					and fields[8] == pos_var4_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				else:
					continue

	# assess annotation of BRCA1_SGE negative variants

	with open(neg_annotated_vcf, "r") as neg_a_vcf:
		for i_neg, line in enumerate(neg_a_vcf):
			if line.startswith("#CHROM"):
				var1 = i_neg + 1

				if i_neg == var1:
					fields = line.split("\t")


					assert fields[0] == "17" and fields[1] == "41199935" \
					and fields[3] == "C" and fields[4] == "T" \
					and fields[8] == neg_var1_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_neg == (var1 + 1):
					fields = line.split("\t")

					assert fields[0] == "17" and fields[1] == "41202400" \
					and fields[3] == "T" and fields[4] == "C" \
					and fields[8] == neg_var2_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_neg == (var1 + 2):
					fields = line.split("\t")

					assert fields[0] == "17"	and fields[1] == "41226675" \
					and fields[3] == "A" and fields[4] == "T" \
					and fields[8] == neg_var3_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])

				elif i_neg == (var1 + 3):
					fields = line.split("\t")

					assert fields[0] == "17"	and fields[1] == "41299122" \
					and fields[3] == "A" and fields[4] == "C" \
					and fields[8] == neg_var4_annot, \
					"#variant INFO wrong \n{}\t{}\t{}\t{}\t{}".format(fields[0], fields[1], fields[3], fields[4], fields[8])
