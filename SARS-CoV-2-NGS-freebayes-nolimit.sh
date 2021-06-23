#!/bin/bash

{
usage="$(basename "$0") [-h] [-l <SRA_list>] [-g <reference_genome.fasta>] [-a <Frequency>] [-t <threads>]
This program will call variants using freebayes in given SRA NGS sequences files to obtain major viral variants.
    -h  show this help text
    -l  File or path to SRA accession list in tabular format
    -g  PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name.
    -a  A number between (0-1) indicating Viral Frequency for Freebayes variant calling (i.e.: 0.5 = 50% viral frequency).
    -t  Number of CPU processors"
options=':hl:g:a:t:'
while getopts $options option; do
  case "$option" in
    h) echo "$usage"; exit;;
    l) l=$OPTARG;;
    g) g=$OPTARG;;
    a) a=$OPTARG;;
    t) t=$OPTARG;;
    :) printf "missing argument for -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
   \?) printf "illegal option: -%s\n" "$OPTARG" >&2; echo "$usage" >&2; exit 1;;
  esac
done

# mandatory arguments
if [ ! "$l" ] || [ ! "$g" ] || [ ! "$a" ] || [ ! "$t" ]; then
  echo "arguments -l, -g, -a and -t must be provided"
  echo "$usage" >&2; exit 1
fi

begin=`date +%s`

echo "Downloading SRA files from the given list of accessions"
prefetch --max-size 800G -O ./ --option-file ${l}
echo "SRA files were downloaded in current directory"
echo ""
echo "Done"
echo ""
echo "Converting SRA files to fastq.gz"
SRA= ls -1 *.sra
for SRA in *.sra; do fastq-dump --gzip ${SRA}
done

##################################################################################
# Trimming downloaded Illumina datasets with fastp, using 16 threads (-w option) #
##################################################################################

echo "Trimming downloaded Illumina datasets with fastp."
echo ""

z= ls -1 *.fastq.gz
for z in *.fastq.gz; do fastp -w ${t} -i ${z} -o ${z}.fastp
gzip ${z}.fastp
done

###########################################################################################
# Aligning illumina datasets againts reference with minimap, using 20 threads (-t option) #
###########################################################################################

echo "Aligning illumina datasets againts reference with minimap, using n threads."
echo ""
b= ls -1 *.fastq.gz.fastp.gz
for b in *.fastq.gz.fastp.gz; do minimap2 -ax sr ${g} ${b} > ${b}.sam -t ${t}
samtools sort ${b}.sam > ${b}.sam.sorted.bam -@ ${t}
rm ${b}.sam
rm ${b}
done

######################
# Renaming BAM files #
######################

echo "Renaming files in bash"
echo ""
for filename in *.bam; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fastq.gz.fastp.gz.sam.sorted//g')";  done

######################
# Indexing BAM files #
######################

echo "Indexing BAM files."
echo ""

f= ls -1 *.bam
for f in *.bam; do samtools index ${f}; done

#################################################
### Performing Variant Calling with freebayes ###
#################################################

echo "Performing Variant Calling with freebayes:"
echo ""

x= ls -1 *.bam
for x in *.bam; do freebayes-parallel <(fasta_generate_regions.py ${g}.fai 2000) ${t} -f ${g} -F ${a} -b ${x} > ${x}.freebayes.vcf
done

#######################################
### Merging variants using jacquard ###
#######################################
echo "Merging variants using jacquard"
echo ""
echo "for information, please see: https://jacquard.readthedocs.io/en/v0.42/overview.html#why-would-i-use-jacquard"
echo ""
# Removing 0-byte files in folder
find . -size 0 -delete
echo "Merge VCFs using jacquard"
echo ""
jacquard merge --include_all ./ merged.vcf

echo "Left only genotypes in merged VCF"
echo ""
vcfkeepgeno merged.vcf GT > merged.GT.vcf

echo "Splitting variants and header from merged.GT.vcf"
echo ""
grep "#" merged.GT.vcf > header
grep -v "#" merged.GT.vcf > variants.vcf
echo ""
sed -i 's|1/1|1|'g variants.vcf   # convert diploid to haploid
sed -i 's|0/1|1|'g variants.vcf   # convert diploid to haploid
sed -i 's|1/0|1|'g variants.vcf   # convert diploid to haploid
sed -i 's/[.]/0/'g variants.vcf   # convert points to zeros
echo ""
# Reconstitute vcf file
echo "Reconstituting vcf file"
echo ""
cat header variants.vcf > merged.fixed.vcf
rm header variants.vcf
sed -i 's/NC_04551202/NC_045512.2/'g merged.fixed.vcf

echo "left-aligning vcf file and fixing names"
echo ""
vcfleftalign -r ${g} merged.fixed.vcf > merged.left.vcf
sed -i 's/|unknown//'g merged.left.vcf

# Calculating Viral Frequencies
echo "Calculating viral frequencies with vcflib"
echo ""
vcffixup merged.left.vcf > merged.AF.raw.vcf
wget https://raw.githubusercontent.com/W-L/ProblematicSites_SARS-CoV2/master/problematic_sites_sarsCov2.vcf
sed -i 's/MN908947.3/NC_045512.2/'g problematic_sites_sarsCov2.vcf
vcfintersect -i problematic_sites_sarsCov2.vcf merged.AF.raw.vcf -r ${g} --invert > merged.AF.vcf
rm merged.fixed.vcf merged.left.vcf
gzip merged.vcf merged.AF.raw.vcf

# Filter variants by Viral Frequency: 0.0099 (1%)
echo "Filtering variants by Viral Frequency: 0.0099 (1%)"
echo ""
vcffilter -f "AF > 0.0099" merged.AF.vcf > merged.AF_1%.vcf
grep -v "##" merged.AF_1%.vcf > merged.AF_1%.table

echo "#######"
echo "Summary:"
echo "#######"
echo ""
echo "merged.AF.raw.vcf.gz contain all merged variants."
echo ""
echo "merged.AF.vcf contain merged variants, problematic sites excluded. See https://virological.org/t/issues-with-sars-cov-2-sequencing-data/473."
echo ""
echo "merged.AF_1%.vcf contain merged variants with viral frequency >=1%."
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
