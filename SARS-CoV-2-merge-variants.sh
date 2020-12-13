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

begin=`date +%s`

### Merge VCFs using jacquard
echo "Merge VCFs using jacquard"
echo ""
jacquard merge --include_all ./ merged.GISAID.vcf
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
vcffixup merged.GISAID.left.vcf > merged.GISAID.AF.raw.vcf
wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf
sed -i 's/MN908947.3/NC_045512.2/'g problematic_sites_sarsCov2.vcf
vcfintersect -i problematic_sites_sarsCov2.vcf merged.GISAID.AF.raw.vcf -r ${2} --invert > merged.GISAID.AF.vcf
rm merged.GISAID.fixed.vcf merged.GISAID.left.vcf
gzip merged.GISAID.vcf merged.GISAID.AF.raw.vcf

# Filter variants by Viral Frequency: 0.0099 (1%)
echo "Filtering variants by Viral Frequency: 0.0099 (1%)"
echo ""
vcffilter -f "AF > 0.0099" merged.GISAID.AF.vcf > merged.GISAID.AF_1%.vcf
grep -v "##" merged.GISAID.AF_1%.vcf > merged.GISAID.AF_1%.table

echo "#######"
echo "Summary:"
echo "#######"
echo ""
echo "merged.GISAID.AF.raw.vcf.gz contain all merged variants"
echo ""
echo "merged.GISAID.AF.vcf contain merged variants, problematic sites excluded. See https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473."
echo ""
echo "merged.GISAID.AF_1%.vcf contain merged variants with viral Frequency >=1%."
echo ""
echo "merged.GISAID.AF_1%.table contain merged variants (Viral Frequency >=1%), without VCF header, suitable for plotting"
echo ""
echo "All done."

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed

mdate=`date +'%d/%m/%Y %H:%M:%S'`
mcpu=$[100-$(vmstat 1 2|tail -1|awk '{print $15}')]%
mmem=`free | grep Mem | awk '{print $3/$2 * 100.0}'`
echo "$mdate | $mcpu | $mmem" >> ./stats-cpu

###############################################################
#
} | tee logfile
#
