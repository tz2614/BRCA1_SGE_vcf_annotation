#!usr/bin/python3

""" create_config_NILES.py

- create the config file for snakemake pipeline for a given sample name, or runfolder (batch mode)
- create an initial logging file for each sample name

Input:
- runfolder (default: current folder)
- sample name
Output:
- {runfolder/current_folder}/logs/{sample}_config.yaml
- {runfolder/current_folder}/logs/{sample}.log

Example:
CP0358/
|-- C011022.vcf

python3.5 BRCA1_SGE_vcf_annotator.py -v C014255b.vcf
python3.5 BRCA1_SGE_vcf_annotator.py -f .

CP0358/
|-- C011022.vcf
|-- C011022.vcf.gz
|-- C011022.vcf.gz.tbi
|-- C011022.BRCA1_annotated.vcf
|-- C011022.BRCA1_annotated_20191024.log
"""

import sys
import os
import subprocess
import datetime
import BRCA1_SGE_ref
import remove_INFO_field
import argparse


def create_log_file(vcf_file, log_info):

    """Create a log file for each vcf file that has been annotated,
    including the name of vcf, assciated files generated and timestamp"""

    timestamp = "_".join(str(datetime.datetime.now()).split(" "))
    vcf_file = os.path.abspath(vcf_file)
    
    filedate = str(datetime.date.today())
    filedate = "".join(filedate.split("-"))
    log_file = vcf_file.replace(".vcf", ".BRCA1_annotated.log")

    with open (log_file, "a+") as log:
        log.write(timestamp + "\n")
        log.writelines("{}\n".format(str(log_info)))

    return log_file

def create_vcf_header(ref_vcf):

    """create a header file containing the reference vcf header"""

    header_file = ref_vcf + ".hdr"

    if not os.path.exists(header_file):
        print ("reference vcf header file DO NOT exist")
        print ("generating reference vcf header file")
        with open (ref_vcf, "r") as vcf_file:
            for line in vcf_file:
                if line.startswith("##INFO=<ID=BRCA1_SGE"):
                    #print ("finding header in vcf")
                    with open (header_file, "w") as hdr:
                        hdr.writelines(line)
                else:
                    continue
        print ("reference vcf header file generated")

    else:
        print ("reference vcf header file already exists")

    return header_file

def sort_vcf_file(ref_vcf):

    """sort the lines in the vcf according to CHROM and POS (numerically)"""

    ref_vcf = os.path.abspath(ref_vcf)
    #print (ref_vcf)
    assert os.path.exists(ref_vcf), "reference vcf file DO NOT exist"

    if "sorted" not in ref_vcf:
        sorted_ref_vcf = ref_vcf.replace(".vcf", ".sorted.vcf")
    else:
        sorted_ref_vcf = ref_vcf

    if not os.path.exists(sorted_ref_vcf):
        print ("sort {} according to CHROM and POS".format(sorted_ref_vcf))
        sort_vcf = "cat {ref_vcf} | awk '$1 ~ /^#/ {{print $0;next}}' > {sorted_ref_vcf}".format(ref_vcf=ref_vcf, sorted_ref_vcf=sorted_ref_vcf)
        subprocess.call(sort_vcf, shell=True)
        sort_vcf2 = "grep -v '#' {ref_vcf} | sort -k1,1V -k2,2n >> {sorted_ref_vcf}".format(ref_vcf=ref_vcf, sorted_ref_vcf=sorted_ref_vcf)
        subprocess.call(sort_vcf2, shell=True)
        #print (sort_ref_vcf)
    else:
        print ("sorted reference vcf already exists")

    return sorted_ref_vcf

def merge_vcf_variants(sorted_ref_vcf, bgzip_file, tabix_file):

    assert os.path.exists(bgzip_file), "{} DO NOT exists".format(bgzip_file)
    assert os.path.exists(tabix_file), "{} DO NOT exists".format(tabix_file)

    merged_bgzip_file = sorted_ref_vcf + ".merged.gz"
    merged_tabix_file = merged_bgzip_file + ".tbi"

    if not os.path.exists(merged_bgzip_file):
        print ("merge multi-allelic variants to create {}".format(merged_bgzip_file))
        merge_variants = "bcftools norm -m +any {bgzip_file} -Oz -o {merged_bgzip_file}".format(bgzip_file=bgzip_file, merged_bgzip_file=merged_bgzip_file)
        subprocess.call(merge_variants, shell=True)
    else:
        print("merged_bgzip file already exists")
    
    if not os.path.exists(merged_tabix_file):
        print ("create {}".format(merged_tabix_file))
        create_tabix = "bcftools index -t {merged_bgzip_file}".format(merged_bgzip_file=merged_bgzip_file)
        subprocess.call(create_tabix, shell=True)
    else:
        print("merged tabix file already exists")

    return sorted_ref_vcf, merged_bgzip_file, merged_tabix_file

def compress_file(vcf_file):

    # compress the file into bgzip file format

    if vcf_file.endswith(".vcf"):
        bgzip_file = vcf_file + ".gz"
    else:
        print ("{} NOT a vcf file".format(vcf_file))
        exit()

    if not os.path.exists(bgzip_file):
        create_bgzip = "cat {vcf_file} | bgzip -c > {bgzip_file}".format(vcf_file=vcf_file, bgzip_file=bgzip_file)
        print("bgzipping {}".format(vcf_file))
        subprocess.call(create_bgzip, shell=True)
        return bgzip_file
    else:
        print ("bgzip file {} already exists".format(bgzip_file))
        return bgzip_file

def index_file(bgzip_file):

    #index the .gz file using tabix

    if bgzip_file.endswith(".gz"):
        tabix_file = bgzip_file + ".tbi"
    else:
        print ("{} NOT a gz file".format(bgzip_file))
        exit()

    if not os.path.exists(tabix_file):
        create_tabix = "tabix -p vcf {bgzip_file}".format(bgzip_file=bgzip_file)
        subprocess.call(create_tabix, shell=True)
        print ("indexing {}".format(bgzip_file))
        return tabix_file
    else:
        print ("tabix file {} already exists".format(tabix_file))
        return tabix_file

def annotate_ref(bgzip_file, vcf_header, ref_vcf):

	# check that the reference vcf.gz exist as a data source and that the vcf file to be annotated has.vcf.gz and vcf.gz.tbi files
    assert bgzip_file.endswith(".gz"), "reference.vcf.gz file DO NOT end with .gz"
    assert vcf_header.endswith(".hdr"), "header file DO NOT end with .hdr"
    assert ref_vcf.endswith(".vcf"), "ref vcf DO NOT end with .vcf"

    assert os.path.exists(bgzip_file), "{} DO NOT exists".format(bgzip_file)
    assert os.path.exists(vcf_header), "{} DO NOT exists".format(vcf_header)
    assert os.path.exists(ref_vcf), "{} DO NOT exists".format(ref_vcf)

    if "BRCA1_annotated" in ref_vcf:
        print("ref vcf already annotated")
        return annotated_vcf_file
    else:
        annotated_vcf_file = ref_vcf.replace(".vcf", ".BRCA1_annotated.vcf")

    print ("create ref vcf with additional BRCA1 annotation")
    command = "bcftools annotate -a {} -h {} -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,+INFO/BRCA1_SGE -o {} {}".format(bgzip_file, vcf_header, annotated_vcf_file, ref_vcf)
    subprocess.call(command, shell=True)

    print ("annotated vcf: {}".format(annotated_vcf_file))
    return annotated_vcf_file

def annotate_vcf(bgzip_file, vcf_header, vcf_file, vcf_bgzip_file, vcf_tabix_file):

    """using bcftools and add BRCA1_SGE annotation to the vcf_file"""

    # check that the annotation vcf.gz exist as a data source and that the vcf file to be annotated has.vcf.gz and vcf.gz.tbi files
    assert bgzip_file.endswith(".gz"), "annotation.gz file DO NOT end with .gz"
    assert vcf_header.endswith(".hdr"), "header file DO NOT end with .hdr"
    assert vcf_file.endswith(".vcf"), "vcf file DO NOT end with .vcf"
    assert vcf_bgzip_file.endswith(".gz"), "vcf.gz file DO NOT end with .gz"
    assert vcf_tabix_file.endswith(".tbi"), " tabix file DO NOT end with .tbi"
    
    assert os.path.exists(bgzip_file), "{} DO NOT exists".format(bgzip_file)
    assert os.path.exists(vcf_header), "{} DO NOT exists".format(vcf_header)
    assert os.path.exists(vcf_file), "{} DO NOT exists".format(vcf_file)
    assert os.path.exists(vcf_bgzip_file), "{} DO NOT exists".format(vcf_bgzip_file)
    assert os.path.exists(vcf_tabix_file), "{} DO NOT exists".format(vcf_tabix_file)

    if "BRCA1_annotated" in vcf_file:
        print("vcf file already annotated")
        return annotated_vcf_file
    else:
        annotated_vcf_file = vcf_file.replace(".vcf", ".BRCA1_annotated.vcf")

    print ("create new vcf with additional BRCA1 annotation")
    command = "bcftools annotate -a {} -h {} -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,+INFO/BRCA1_SGE -Ov -o {} {}".format(bgzip_file, vcf_header, annotated_vcf_file, vcf_bgzip_file)
    subprocess.call(command, shell=True)

    print ("annotated vcf: {}".format(annotated_vcf_file))
    return annotated_vcf_file

def sample_annotation(ref_csv, vcf_file, runfolder=None):

    # generate BRCA1_SGE_ref.vcf using 41586_2018_461_MOESM3_ESM.csv
    # sort, bgzip, index and merge BRCA1_SGE_ref.vcf as reference file containing all the BRCA1_SGE variants to annotate new vcf

    if runfolder == ".":
        rf_path = os.getcwd()

    if not runfolder:
        rf_path = os.path.dirname(vcf_file)

    vcf_file = os.path.basename(vcf_file)

    try:
        os.path.getsize("{}/{}".format(rf_path, vcf_file))

    except FileNotFoundError:
        print ("{}/{} doesn't exist".format(rf_path, vcf_file))
        return -1
        
    if os.path.getsize("{}/{}".format(rf_path, vcf_file)) == 0:
        print ("{}/{} is empty".format(rf_path, vcf_file))
        return -1

    if vcf_file.startswith("C01"):

        vcf_file = os.path.abspath(vcf_file)
        log = create_log_file(vcf_file, vcf_file)

        ref_csv = os.path.abspath(ref_csv)
        log = create_log_file(vcf_file, ref_csv)

        ref_vcf = BRCA1_SGE_ref.create_ref_vcf(ref_csv)
        log = create_log_file(vcf_file, ref_vcf)

        bgzip_file, tabix_file = sort_bgzip_index(ref_vcf)
        log = create_log_file(vcf_file, bgzip_file)
        log = create_log_file(vcf_file, tabix_file)

        vcf_header = create_vcf_header(ref_vcf)
        log = create_log_file(vcf_file, vcf_header)
       
        ref_vcf_remove_info = remove_INFO_field.remove_INFO(ref_vcf)
        log = create_log_file(vcf_file, ref_vcf_remove_info)

        annotated_ref_vcf = annotate_ref(bgzip_file, vcf_header, ref_vcf_remove_info)
        log = create_log_file(vcf_file, annotated_ref_vcf)

        annotated_sorted_ref_vcf = sort_vcf_file(annotated_ref_vcf)
        log = create_log_file(vcf_file, annotated_sorted_ref_vcf)
       
        annotated_ref_bgzip_file = compress_file(annotated_sorted_ref_vcf)
        log = create_log_file(vcf_file, annotated_ref_bgzip_file)

        annotated_ref_tabix_file = index_file(annotated_ref_bgzip_file)       
        log = create_log_file(vcf_file, annotated_ref_tabix_file)

        annotated_sorted_ref_vcf, merged_bgzip_file, merged_tabix_file = merge_vcf_variants(annotated_sorted_ref_vcf, annotated_ref_bgzip_file, annotated_ref_tabix_file)
        log = create_log_file(vcf_file, merged_bgzip_file)
        log = create_log_file(vcf_file, merged_tabix_file)

        vcf_bgzip_file = compress_file(vcf_file)
        log = create_log_file(vcf_file, vcf_bgzip_file)

        vcf_tabix_file = index_file(vcf_bgzip_file)
        log = create_log_file(vcf_file, vcf_tabix_file)

        annotated_final_vcf = annotate_vcf(merged_bgzip_file, vcf_header, vcf_file, vcf_bgzip_file, vcf_tabix_file)
        log = create_log_file(vcf_file, annotated_final_vcf)
    
def batch_annotation(ref_csv, runfolder="."):

    C01_samples = []
    
    # if no runfolder specified 
    # get runfolder for batch annotation

    if runfolder == ".":
        rf_path = os.getcwd()
    else:
        # make sure path is absolute_path
        rf_path = os.path.abspath(runfolder)
    
    #check runfolder exists
    if os.path.isdir(rf_path):
         #get the sample names from each vcf files
        for root, dirname, filenames in os.walk(rf_path):
            for filename in filenames:
                #find all paths to C01####.vcf in the runfolder
                if filename.endswith(".vcf") and "BRCA1_annotated" not in filename and filename.startswith("C01"):
                    C01_samples.append(os.path.join(root, filename))           
    else:
        print ("{} is not a directory, check the runfolder given is correct".format(rf_path))

    if C01_samples == []:
        print ("No unannotated C01 samples with C01#####.vcf found in this runfolder")
        return -1

    else:
        # check both vcf and associated output files exist and contain data or not
        C01_samples = sorted(C01_samples)

        # annotate each sample in the vcf file
        for sample in C01_samples:
            sample_annotation(ref_csv, sample)

#if __name__ == "__main__":
    #main(sys.argv[1], sys.argv[2])

# navigate to where the script is located from command line, and enter "python3", 

# the name of the script "BRCA1_SGE_vcf_annotator.py"

# sys.argv[1] - /path/to/BRCA1_SGE_ref.vcf
# sys.argv[2] - /path/to/vcf you wish to annotate

def main():
    
    parser = argparse.ArgumentParser()

    parser.add_argument("-v", "--vcf", help = "run BRCA1_SGE annotation for this vcf")
    parser.add_argument("-f", "--runfolder", help = "run BRCA1_SGE annotation for all vcf files in this runfolder")
    parser.add_argument("-r", "--BRCA1_SGE_reference", help = "Reference csv file containing BRCA1_SGE annotations e.g. 41586_2018_461_MOESM3_ESM.csv for use in the analysis")
    args = parser.parse_args()

    if args.vcf and args.runfolder:
        sample_annotation(args.BRCA1_SGE_reference, args.vcf, runfolder=args.runfolder)
    elif args.runfolder:
        batch_annotation(args.BRCA1_SGE_reference, runfolder=args.runfolder)
    elif args.vcf:
        sample_annotation(args.BRCA1_SGE_reference, args.vcf, runfolder=".")

    else:
        parser.print_help()
        exit(-1)

if __name__ == "__main__":
    main()