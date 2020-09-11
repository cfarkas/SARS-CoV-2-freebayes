#!/bin/bash

{

SRA_list=${1}
Reference=${2}
Threads=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  echo ""
  echo "This program will call variants using Strelka in given SRA NGS sequences files to obtain viral founder variants."
  echo ""
  echo "[SRA_list]: File of path to SRA accession list in tabular format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"; exit 1; }

if [ $# -ne 3 ]; then
  echo 1>&2 "Usage: ./`basename $0` [SRA_list] [Reference] [Threads]"
  exit 3
fi
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)

echo "Downloading SRA files from the given list of accessions"
prefetch --max-size 800G -O ./ --option-file ${1}
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

a= ls -1 *.fastq.gz
for a in *.fastq.gz; do fastp -w ${3} -i ${a} -o ${a}.fastp
gzip ${a}.fastp
done

###########################################################################################
# Aligning illumina datasets againts reference with minimap, using 20 threads (-t option) #
###########################################################################################

echo "Aligning illumina datasets againts reference with minimap, using n threads."
echo ""

b= ls -1 *.fastq.gz.fastp.gz
for b in *.fastq.gz.fastp.gz; do minimap2 -ax sr ${2} ${b} > ${b}.sam -t ${3}
samtools sort ${b}.sam > ${b}.sam.sorted.bam -@ ${3}
rm ${b}.sam
rm ${b}
done

######################
# Renaming BAM files #
######################

echo "Renaming files in bash"
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

a= ls -1 *.bam
for a in *.bam; do freebayes-parallel <(fasta_generate_regions.py ${2}.fai 2000) ${3} -f covid19-refseq.fasta -F 0.49 -b ${a} > ${a}.freebayes.vcf
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
ulimit -n 1000000
jacquard merge --include_all ./ merged.vcf

echo "Left only genotypes in merged VCF"
echo ""
vcfkeepgeno merged.vcf GT > merged.GT.vcf

echo "Split variants and header from merged.GT.vcf"
echo ""
grep "#" merged.GT.vcf > header
grep -v "#" merged.GT.vcf > variants.vcf

sed -i 's|1/1|1|'g variants.vcf   # convert diploid to haploid
sed -i 's|0/1|1|'g variants.vcf   # convert diploid to haploid
sed -i 's|1/0|1|'g variants.vcf   # convert diploid to haploid
sed -i 's/[.]/0/'g variants.vcf   # convert points to zeros

# Reconstitute vcf file
cat header variants.vcf > merged.fixed.vcf
rm header variants.vcf
sed -i 's/NC_04551202/NC_045512.2/'g merged.fixed.vcf

echo "left-align vcf file and fix names"
echo ""
vcfleftalign -r /home/user/MITACS/July_12_2020/SARS-CoV-2_illumina_analysis/1/covid19-refseq.fasta merged.fixed.vcf > merged.left.vcf
sed -i 's/|unknown//'g merged.left.vcf

echo "calculate AF with vcflib"
echo ""
vcffixup merged.left.vcf > merged.AF.vcf
rm merged.fixed.vcf merged.left.vcf

echo "Filter variants by AF: 0.0099 (1%, founders)"
echo ""
vcffilter -f "AF > 0.0099"  merged.AF.vcf > merged.AF_0.01.vcf

echo "Filter variants by AF: 0.0049 (0.5%)"
echo ""
vcffilter -f "AF > 0.00499"  merged.AF.vcf > merged.AF_0.005.vcf

echo "All done. Merged vcf files are called merged.AF.vcf, merged.AF_0.01.vcf and merged.AF_0.005.vcf and are located in current directory"

###############################################################
#
} | tee logfile
#
