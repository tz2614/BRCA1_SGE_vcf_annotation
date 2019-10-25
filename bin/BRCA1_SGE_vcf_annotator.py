#!usr/bin/python3

import sys
import os
import subprocess
import datetime
import BRCA1_SGE_ref
import remove_INFO_field

def create_vcf_header(ref_vcf):

    """create a header file containing the reference vcf header"""

    header_file = ref_vcf + ".hdr"

    if os.path.exists(header_file):
        print ("header file already exists")
        return header_file
    else:
        print ("creating header file")
        with open (ref_vcf, "r") as vcf_file:
            for line in vcf_file:
                if line.startswith("##INFO=<ID=BRCA1_SGE"):
                    #print ("finding header in vcf")
                    with open (header_file, "w") as hdr:
                        hdr.writelines(line)
                else:
                    continue

    print ("header file {} generated".format(header_file))
    return header_file


def sort_bgzip_index(ref_vcf):

    """sort the lines in the vcf according to CHROM and POS (numerically), compress the file into bgzip file format, and then index the .gz file"""

    ref_vcf = os.path.abspath(ref_vcf)
    assert os.path.exists(ref_vcf), "reference vcf file DO NOT exist"

    bgzip_file = ref_vcf + ".gz"
    tabix_file = bgzip_file + ".tbi"

    if not os.path.exists(bgzip_file):
        print ("sort {} according to CHROM and POS, and create {}".format(ref_vcf, bgzip_file))
        sort_bgzip = "grep -v '#' {} | sort -k1V -k2n -t$'\t' | bgzip -c > {};".format(ref_vcf, bgzip_file)
        subprocess.call(sort_bgzip, shell=True)
    else:
        print("reference vcf.gz file already exists".format(bgzip_file))

    if not os.path.exists(tabix_file):
        print ("create {}".format(tabix_file))
        create_tabix = "tabix -p vcf {}".format(bgzip_file)
        subprocess.call(create_tabix, shell=True)
    else:
        print("reference.vcf.gz.tbi file already exists".format(tabix_file))

    return bgzip_file, tabix_file

def annotate_ref(bgzip_file, vcf_header, ref_vcf):

    """using bcftools and add BRCA1_SGE annotation to the reference vcf, using bgzip_file, vcf_header, BRCA1_vcf and vcf_file as arguments"""

    # check that the reference vcf.gz exist as a data source
    assert bgzip_file.endswith(".gz"), "reference.vcf.gz file DO NOT end with .gz"
    assert vcf_header.endswith(".hdr"), "header file DO NOT end with .hdr"
    assert ref_vcf.endswith(".vcf"), "ref vcf DO NOT end with .vcf"

    assert os.path.exists(bgzip_file), "reference vcf.gz file DO NOT exists".format(bgzip_file)
    assert os.path.exists(vcf_header), "header file DO NOT exists".format(vcf_header)
    assert os.path.exists(ref_vcf), "ref vcf DO NOT exists".format(ref_vcf)

    annotated_ref_vcf = ref_vcf.replace(".vcf", ".BRCA1_annotated.vcf")

    if os.path.exists(annotated_ref_vcf):
        print ("annotated ref vcf already exists")
    else:
	    print ("create ref vcf with additional BRCA1 annotation")
	    command = "bcftools annotate -a {} -h {} -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,+INFO/BRCA1_SGE -o {} {}".format(bgzip_file, vcf_header, annotated_ref_vcf, ref_vcf)
	    subprocess.call(command, shell=True)

    print ("annotated ref vcf: {}".format(annotated_ref_vcf))
    return annotated_ref_vcf

def annotate_vcf(bgzip_file, vcf_header, vcf_file, vcf_bgzip_file, vcf_tabix_file):
    
    """using bcftools and add BRCA1_SGE annotation to the vcf_file"""

    # check that the annotation vcf.gz exist as a data source and that the vcf file to be annotated has.vcf.gz and vcf.gz.tbi files
    assert bgzip_file.endswith(".gz"), "merged.gz file DO NOT end with .gz"
    assert vcf_header.endswith(".hdr"), "header file DO NOT end with .hdr"
    assert vcf_file.endswith(".vcf"), "vcf file DO NOT end with .vcf"
    assert vcf_bgzip_file.endswith(".gz"), "vcf.gz file DO NOT end with .gz"
    assert vcf_tabix_file.endswith(".tbi"), " tabix file DO NOT end with .tbi"
    
    assert os.path.exists(bgzip_file), "{} DO NOT exists".format(bgzip_file)
    assert os.path.exists(vcf_header), "{} DO NOT exists".format(vcf_header)
    assert os.path.exists(vcf_file), "{} DO NOT exists".format(vcf_file)
    assert os.path.exists(vcf_bgzip_file), "{} DO NOT exists".format(vcf_bgzip_file)
    assert os.path.exists(vcf_tabix_file), "{} DO NOT exists".format(vcf_tabix_file)

    annotated_vcf = vcf_file.replace(".vcf", ".BRCA1_annotated.vcf")

    if os.path.exists(annotated_vcf):
        print ("annotated vcf already exists")
    else:
        print ("create new vcf with additional BRCA1 annotation")
        command = "bcftools annotate -a {} -h {} -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,+INFO/BRCA1_SGE -Ov -o {} {}".format(bgzip_file, vcf_header, annotated_vcf, vcf_bgzip_file)
        subprocess.call(command, shell=True)
        print ("annotated vcf: {}".format(annotated_vcf))

    return annotated_vcf

def sort_vcf_file(ref_vcf):

    """sort the lines in the vcf according to CHROM and POS (numerically)"""

    assert os.path.exists(ref_vcf), "reference vcf file {} DO NOT exist for sorting".format(ref_vcf)

    if "sorted" not in ref_vcf:
        sorted_ref_vcf = ref_vcf.replace(".vcf", ".sorted.vcf")
    else:
        sorted_ref_vcf = ref_vcf

    if not os.path.exists(sorted_ref_vcf):
        print ("sort {} according to CHROM and POS".format(sorted_ref_vcf))
        sort_vcf = "cat {ref_vcf} | awk '$1 ~ /^#/ {{print $0;next}}' > {sorted_ref_vcf}".format(ref_vcf=ref_vcf, sorted_ref_vcf=sorted_ref_vcf)
        subprocess.call(sort_vcf, shell=True)
        sort_vcf2 = "grep -v '#' {ref_vcf} | sort -k1V -k2n >> {sorted_ref_vcf}".format(ref_vcf=ref_vcf, sorted_ref_vcf=sorted_ref_vcf)
        subprocess.call(sort_vcf2, shell=True)
        #print (sort_ref_vcf)
    else:
        print ("sorted reference vcf already exists")

    return sorted_ref_vcf

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
        print ("{} already exists".format(bgzip_file))
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
        print ("{} already exists".format(tabix_file))
        return tabix_file

def merge_vcf_variants(ref_vcf, bgzip_file, tabix_file):

    # merge vcf records at the same position, ready for annotating mult-allelic variants

    assert os.path.exists(bgzip_file), "{} DO NOT exists".format(bgzip_file)
    assert os.path.exists(tabix_file), "{} DO NOT exists".format(tabix_file)
    assert os.path.exists(ref_vcf), "{} DO NOT exists".format(ref_vcf)

    merged_bgzip_file = ref_vcf + ".merged.gz"
    merged_tabix_file = merged_bgzip_file + ".tbi"

    if not os.path.exists(merged_bgzip_file):
        # merge multi-allelic variants in reference vcf.gz
        print ("merge multi-allelic variants")
        merge_variants = "bcftools norm -m +any {bgzip_file} -Oz -o {merged_bgzip_file}".format(bgzip_file=bgzip_file, merged_bgzip_file=merged_bgzip_file)
        subprocess.call(merge_variants, shell=True)
        print ("{} created".format(merged_bgzip_file))
    else:
        print("merged_bgzip file already exists")
    
    if not os.path.exists(merged_tabix_file):
        # indexing merged vcf.gz
        print ("indexing merged gz")
        create_tabix = "bcftools index -t {merged_bgzip_file}".format(merged_bgzip_file=merged_bgzip_file)
        subprocess.call(create_tabix, shell=True)
        print ("{} created".format(merged_tabix_file))
    else:
        print("merged tabix file already exists")

    return merged_bgzip_file, merged_tabix_file

def main(ref_vcf, vcf_file):

    # use BRCA1_SGE_ref.vcf as reference file containing all the BRCA1_SGE variants to annotate new vcf

    ref_vcf = os.path.abspath(ref_vcf)
    vcf_file = os.path.abspath(vcf_file)
    vcf_header = create_vcf_header(ref_vcf)
    bgzip_file, tabix_file = sort_bgzip_index(ref_vcf)
    #remove INFO field in ref vcf to facilitate enable BRCA1_SGE annotations to be added to ref vcf
    removed_info_ref_vcf = remove_INFO_field.remove_INFO(ref_vcf)
    annotated_ref_vcf = annotate_ref(bgzip_file, vcf_header, removed_info_ref_vcf)
    sorted_annotated_ref_vcf = sort_vcf_file(annotated_ref_vcf)
    ref_bgzip_file = compress_file(sorted_annotated_ref_vcf)
    ref_tabix_file = index_file(ref_bgzip_file)
    merged_bgzip_file, merged_tabix_file = merge_vcf_variants(sorted_annotated_ref_vcf, ref_bgzip_file, ref_tabix_file)

    # compress and index existing vcf before annotating vcf file
    bgzip_file = compress_file(vcf_file)
    tabix_file = index_file(bgzip_file)

    # annotating vcf with BRCA1 annotations
    annotated_vcf = annotate_vcf(merged_bgzip_file, vcf_header, vcf_file, bgzip_file, tabix_file)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

# navigate to where the script is located from command line, and enter "python3", 

# the name of the script "BRCA1_SGE_vcf_annotator.py"

# sys.argv[1] - /path/to/BRCA1_SGE_ref.vcf
# sys.argv[2] - /path/to/vcf you wish to annotate