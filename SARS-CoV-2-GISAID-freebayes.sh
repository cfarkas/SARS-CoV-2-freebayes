#!/bin/bash
set -e 

{

GISAID_fasta=${1}
Reference=${2}
Threads=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  exit 3
fi

begin=`date +%s`

### Split fasta files
echo "fixing names in FASTA file"
echo ""
sed -i 's/ /-/'g ${1}
sed -i "s|hCoV-19/.*./2020||"g ${1}
sed -i "s|hCoV-19/.*./2019||"g ${1}
sed -i 's/|/\t/'g ${1}
sed -i 's/>\t/>/'g ${1}
seqkit fx2tab ${1} > merged.GISAID.tabular
awk '{print $1"\t"$3}' merged.GISAID.tabular > merged.GISAID.tab && rm merged.GISAID.tabular
seqkit tab2fx merged.GISAID.tab > merged.GISAID.fasta && rm merged.GISAID.tab
echo "Splitting fasta files with seqkit"
echo ""
seqkit split --by-id merged.GISAID.fasta
cd merged.GISAID.fasta.split/
for filename in *.fasta; do mv "./$filename" "./$(echo "$filename" | sed -e 's/merged.GISAID.id_//g')";  done
cd ..
cp ./merged.GISAID.fasta.split/EPI*.fasta ./
rm -r -f merged.GISAID.fasta.split merged.GISAID.fasta ${1}
echo "Split is done. Continue with FASTA alignments"
echo ""

### Align fasta files to reference (covid19-refseq.fasta, provided in this repository) and call variants with freebayes (option C 1)
echo "Aligning fasta files to reference and call variants with freebayes (option C 1)"
echo ""
samtools faidx ${2}
fasta= ls -1 *.fasta
for fasta in *.fasta; do
minimap2 -ax asm5 -t ${3} ${2} ${fasta} > ${fasta}.sam
samtools view -bS ${fasta}.sam > ${fasta}.bam
samtools sort -o ${fasta}.sorted.bam ${fasta}.bam
freebayes -f ${2} -C 1 ${fasta}.sorted.bam > ${fasta}.vcf
vcfleftalign -r ${2} ${fasta}.vcf > ${fasta}.left.vcf
rm ${fasta}.sam ${fasta}.bam ${fasta}.sorted.bam ${fasta}.vcf
done
echo "Done. Continue with variant aggregation"
echo ""

### Merge of Variants
echo "Merging variants with Jacquard"
echo ""

# fixing VCF files for merge
echo "fixing VCF files for merge"
echo ""
vcf= ls -1 *.fasta.left.vcf; for vcf in *.fasta.left.vcf; do sed -i "s|0/0|1/1|"g ${vcf}; done

# Renaming files in bash
echo "Renaming files in bash"
echo ""
for filename in *.fasta.left.vcf; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fasta.left.vcf/.vcf/g')";  done

# Merge VCFs using jacquard
echo "Merge VCFs using jacquard"
echo ""
ulimit -n 1000000 && ulimit -s 299999 && jacquard merge --include_all ./ merged.GISAID.vcf
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
vcfleftalign -r ${2} merged.GISAID.fixed.vcf > merged.GISAID.left.vcf
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

###############################################################
#
} | tee logfile
#
