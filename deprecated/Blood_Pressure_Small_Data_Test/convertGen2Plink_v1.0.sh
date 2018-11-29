#!/bin/bash/

#----------------------------------------------------------------------------------------------#
# --- convertGen2Plink_v1.0.sh
# --- Stephen Newhouse
# --- stephen.newhouse@kcl.ac.uk
#----------------------------------------------------------------------------------------------#


# Convert Oxford *.gen to PLINK Binary using https://www.cog-genomics.org/plink2/
#
# Assumes
#
# gen (Oxford genotype file format)
#
# Native text genotype file format for Oxford statistical genetics tools, such as IMPUTE2 and SNPTEST. Should always be accompanied by a .sample file. Loaded with --data/--gen, and produced by '--recode oxford'.
# A text file with no header line, and one line per variant with 3N+5 fields where N is the number of individuals. Each line stores information for a single SNP. The first five fields are:
#
# Chromosome code (that's what PLINK assumes, anyway; this field can technically be repurposed)
# RS ID of SNP
# Position in base-pair units
# Allele 1 (usually minor)
# Allele 2 (usually major)
#
# Each subsequent triplet of values then indicate likelihoods of homozygote A1, heterozygote, and homozygote A2 genotypes at this SNP, respectively, for one individual. If they add up to less than one, the remainder is a no-call probability weight.
# PLINK2
# URL :- https://www.cog-genomics.org/plink2/
# USAGE:- sh convertGen2Plink_v1.0.sh [GENFILE_PREFIX] [SAMPLE_FILE] [OUTFILE_PREFIX]

#----------------------------------------------------------------------------------------------#
# Set PLINK exe                                                                                #
#----------------------------------------------------------------------------------------------#
PLINK_EXE="~/plink2"

#----------------------------------------------------------------------------------------------#
# get options                                                                                  #
#----------------------------------------------------------------------------------------------#
gen_file=${1}
sample_file=${2}
out_file=${3}
#----------------------------------------------------------------------------------------------#
# BGEN FORMAT:-                                                                                #
#----------------------------------------------------------------------------------------------#
# IF BGEN FORMAT THEN USE THIS :-
# To specify that chromosome codes should be read from the .bgen 'SNP ID' field, use the 'snpid-chr' modifier. (Without this modifier, the SNP ID field is ignored.)
# See :- https://www.cog-genomics.org/plink2/input#oxford

# Remove the HASH's below and HASH out the STANDARD GEN Section

 echo ""
 echo "Converting to PLINK Binary"
 echo "calling:- plink2 --bgen ${1} --sample ${2} --make-bed --out ${3} --hard-call-threshold 0.0"
 echo ""
 ~/plink2 --bgen ${1} --sample ${2} --missing-code 'NA,NaN' --make-bed --out ${3} --hard-call-threshold 0.0;
#----------------------------------------------------------------------------------------------#
# STANDARD GEN FORMAT:-                                                                        #
#----------------------------------------------------------------------------------------------#

# echo ""
# echo "Converting to PLINK Binary"
# echo "calling:- plink2 --gen ${1} --sample ${2} --make-bed --out ${3} --hard-call-threshold 0.0"
# echo ""

# ${PLINK_EXE} --gen ${1} --sample ${2} --make-bed --out ${3} --hard-call-threshold 0.0;






