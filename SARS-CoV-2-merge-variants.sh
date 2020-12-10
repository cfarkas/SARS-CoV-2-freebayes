#!/bin/bash
set -e 

{

Reference=${1}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Reference]"
  echo ""
  echo "This script will merge variants using jacquard and the calculate viral frequencies"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Reference]"
  echo ""
  echo "This script will merge variants using jacquard and the calculate viral frequencies"
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Reference]"
  echo ""
  echo "This script will merge variants using jacquard and the calculate viral frequencies"
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [Reference]"
  echo ""
  echo "This script will merge variants using jacquard and the calculate viral frequencies"
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [Reference]"; exit 1; }

if [ $# -ne 1 ]; then
  echo 1>&2 "Usage: ./`basename $0` [Reference]"
  exit 3
fi

### Calculating Number of Variants per genome
echo ""
echo "Calculating Number of Variants per genome"
echo ""
ulimit -s 299999
{
vcf= ls -1 *.vcf
for vcf in *.vcf; do grep -P 'NC_045512.2\t' ${vcf} -c
done
#
} | tee logfile_variants_GISAID_freebayes
#
grep "EPI_ISL_" logfile_variants_GISAID_freebayes > vcf_files
grep -v "EPI_ISL_" logfile_variants_GISAID_freebayes > variants_per_sample
paste vcf_files variants_per_sample > logfile_variants_GISAID
rm vcf_files variants_per_sample
sed -i 's/.fa.left.vcf//'g logfile_variants_GISAID
rm logfile_variants_GISAID_freebayes

# Merge VCFs using jacquard
echo "Merge VCFs using jacquard"
echo ""
ulimit -n 1000000 && jacquard merge --include_all ./ merged.GISAID.vcf
echo ""
# Left only genotypes in merged VCF
echo "Fixing genotypes in merged VCF"
echo ""
vcfkeepgeno merged.GISAID.vcf GT > merged.GISAID.GT.vcf

# Split variants and header from merged.GT.vcf
echo "Splitting variants and header from merged.GT.vcf"
echo ""
grep "#" merged.GISAID.GT.vcf > header
grep -v "#" merged.GISAID.GT.vcf > variants.vcf

sed -i 's|1/1|1|'g variants.vcf   
sed -i 's|0/1|1|'g variants.vcf   
sed -i 's|1/0|1|'g variants.vcf  
sed -i 's/[.]/0/'g variants.vcf   # convert point to zeros 

# Reconstitute vcf file
echo "Reconstituting vcf file"
echo ""
cat header variants.vcf > merged.GISAID.fixed.vcf
rm header variants.vcf
sed -i 's/NC_04551202/NC_045512.2/'g merged.GISAID.fixed.vcf

# left-align vcf file and fix names
echo "left-aligning vcf file and fix names"
echo ""
vcfleftalign -r ${1} merged.GISAID.fixed.vcf > merged.GISAID.left.vcf
sed -i 's/|unknown//'g merged.GISAID.left.vcf

# calculate AF
echo "calculating viral frequency with vcflib"
echo ""
vcffixup merged.GISAID.left.vcf > merged.GISAID.AF.vcf
rm merged.GISAID.fixed.vcf merged.GISAID.left.vcf
gzip merged.GISAID.vcf
ulimit -s 299999
gzip *.fasta 

# Filter variants by Viral Frequency: 0.0099 (1%)
echo "Filtering variants by Viral Frequency: 0.0099 (1%)"
echo ""
vcffilter -f "AF > 0.0099" merged.GISAID.AF.vcf > merged.GISAID.AF_1%.vcf
grep -v "##" merged.GISAID.AF_1%.vcf > merged.GISAID.AF_1%.table

echo "#######"
echo "Summary:"
echo "#######"
echo ""
echo "merged.GISAID.AF.vcf contain merged variants"
echo ""
echo "merged.GISAID.AF_1%.vcf contain merged variants (Viral Frequency >=1%)"
echo ""
echo "merged.GISAID.AF_1%.table contain merged variants (Viral Frequency >=1%), without VCF header and suitable for plotting"
echo ""
echo "All done."

###############################################################
#
} | tee logfile
#
