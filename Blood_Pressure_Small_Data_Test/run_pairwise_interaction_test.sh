FILE_PREFIX=${1}
TARGET_SNP=${2}
#Target gene which mimics a drug
rm -f $FILE_PREFIX".set"
echo -e "TARGET_SNP\n"$TARGET_SNP"\nEND\n" >> $FILE_PREFIX".set"
#All other SNPs
echo -e "OTHER_SNPS" >> $FILE_PREFIX".set"
awk '$2!=$TARGET_SNP {print$2}' $FILE_PREFIX".bim" >> $FILE_PREFIX".set"
echo -e "END" >> $FILE_PREFIX".set"
#Run PLINK test for interaction
plink --file $FILE_PREFIX  --epistasis --allow-no-sex --epi1 1 --set-test --set $FILE_PREFIX".set"
