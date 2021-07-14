#!/bin/bash
set -e
{
usage="$(basename "$0") [-h] [-f <sequences.fasta>] [-g <reference_genome.fasta>] [-t <threads>]
This program will call variants using freebayes in given  fasta sequences files to obtain major viral variants.
Arguments:
    -h  show this help text
    -f  genomes in FASTA format
    -g  PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name.
    -t  Number of CPU processors"
options=':hf:g:t:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    f) f=$OPTARG;;
    g) g=$OPTARG;;
    t) t=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$f" ] || [ ! "$g" ] || [ ! "$t" ]; then
  echo "arguments -f, -g and -t must be provided"
  echo "$usage" >&2; exit 1
fi

begin=`date +%s`


### Split fasta files
echo "fixing names in FASTA file"
echo ""
fasta_name=$(echo ${f} | sed "s/.fasta//")
cp ${f} ./merged.sequence.fasta
#
echo "Splitting fasta files with faidx (python)"
echo ""
mkdir ${fasta_name} && mv merged.sequence.fasta ./${fasta_name}
cd ${fasta_name}
faidx --split-files merged.sequence.fasta
rm merged.sequence.fasta merged.sequence.fasta.fai
echo "Split is done. Continue with FASTA alignments"
echo ""

### Align fasta files to reference (covid19-refseq.fasta, provided in this repository) and call variants with freebayes (option C 1)
echo "Aligning fasta files to reference and call variants with freebayes (option C 1)"
echo ""
samtools faidx ${g}
fasta= ls -1 *.fasta
for fasta in *.fasta; do
minimap2 -ax asm5 -t ${t} ${g} ${fasta} > ${fasta}.sam
samtools view -bS ${fasta}.sam > ${fasta}.bam
samtools sort -o ${fasta}.sorted.bam ${fasta}.bam
freebayes -f ${g} -C 1 ${fasta}.sorted.bam > ${fasta}.vcf
vcfleftalign -r ${g} ${fasta}.vcf > ${fasta}.left.vcf
rm ${fasta}.sam ${fasta}.bam ${fasta}.sorted.bam ${fasta}.vcf
done
echo "Done. Continue with variant aggregation"
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
jacquard merge --include_all ./ merged.vcf
echo ""
# Left only genotypes in merged VCF
echo "Fixing genotypes in merged VCF"
echo ""
vcfkeepgeno merged.vcf GT > merged.GT.vcf

# Split variants and header from merged.GT.vcf
echo "Splitting variants and header from merged.GT.vcf"
echo ""
grep "#" merged.GT.vcf > header
grep -v "#" merged.GT.vcf > variants.vcf

sed -i 's|1/1|1|'g variants.vcf
sed -i 's|0/1|1|'g variants.vcf
sed -i 's|1/0|1|'g variants.vcf
sed -i 's/[.]/0/'g variants.vcf   # convert point to zeros

# Reconstitute vcf file
echo "Reconstituting vcf file"
echo ""
cat header variants.vcf > merged.fixed.vcf
rm header variants.vcf
sed -i 's/NC_04551202/NC_045512.2/'g merged.fixed.vcf

# left-align vcf file and fix names
echo "left-aligning vcf file and fixing names"
echo ""
vcfleftalign -r ${g} merged.fixed.vcf > merged.left.vcf
sed -i 's/|unknown//'g merged.left.vcf

# calculate AF
echo "calculating viral frequency with vcflib"
echo ""
vcffixup merged.left.vcf > merged.AF.raw.vcf
wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf
sed -i 's/MN908947.3/NC_045512.2/'g problematic_sites_sarsCov2.vcf
vcfintersect -i problematic_sites_sarsCov2.vcf merged.AF.raw.vcf -r ${g} --invert > merged.AF.vcf
rm merged.fixed.vcf merged.left.vcf
gzip merged.vcf merged.AF.raw.vcf
gzip *.fasta

# Filter variants by Viral Frequency: 0.0099 (1%)
echo "Filtering variants by Viral Frequency: 0.0099 (1%)"
echo ""
vcffilter -f "AF > 0.0099" merged.AF.vcf > merged.AF_1%.vcf
grep -v "##" merged.AF_1%.vcf > merged.AF_1%.table

echo "#######"
echo "Summary:"
echo "#######"
echo "in folder:"
echo ""
echo "./"${fasta_name}
echo ""
echo "The outputs:"
echo ""
echo "merged.AF.raw.vcf.gz contain all merged variants as genotype array, including all samples"
echo ""
echo "merged.AF.vcf contain merged variants, problematic sites excluded. See https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473."
echo ""
echo "merged.AF_1%.vcf contain merged variants with viral Frequency >=1%."
echo ""
echo "merged.AF_1%.table contain merged variants (Viral Frequency >=1%), without VCF header, suitable for plotting"
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
