#!/bin/bash

set -e 
usage="$(basename "$0") [-h] [-f <GISAID.fasta>] [-g <reference_genome.fasta>] [-t <threads>]
This program will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants.
Arguments:
    -h  show this help text
    -f  GISAID genomes in FASTA format
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

### Split fasta files
echo "fixing names in FASTA file"
echo ""
fasta_name=$(echo ${f} | sed "s/.fasta//")
sed -i 's/ /-/'g ${f}
sed -i "s|hCoV-19/.*./2021||"g ${f}
sed -i "s|hCoV-19/.*./2020||"g ${f}
sed -i "s|hCoV-19/.*./2019||"g ${f}
sed -i 's/|/\t/'g ${f}
sed -i 's/>\t/>/'g ${f}
seqkit fx2tab ${f} > merged.GISAID.tabular
awk '{print $1"\t"$NF}' merged.GISAID.tabular > merged.GISAID.tab && rm merged.GISAID.tabular
grep "EPI_" merged.GISAID.tab > merged.GISAID.tabular && rm merged.GISAID.tab
seqkit tab2fx merged.GISAID.tabular > merged.GISAID.fasta && rm merged.GISAID.tabular
echo "Splitting fasta files with faidx (python)"
echo ""
mkdir ${fasta_name} && mv merged.GISAID.fasta ./${fasta_name}
cd ${fasta_name}
faidx --split-files merged.GISAID.fasta
echo "Split is done. Continue with FASTA alignments"
echo ""

### Align fasta files to reference (covid19-refseq.fasta, provided in this repository) and call variants with freebayes (option C 1)
echo "Aligning fasta files to reference and call variants with freebayes (option C 1)"
echo ""
ulimit -n 1000000 && ulimit -s 1000000
samtools faidx ${g}
fasta= ls -1 EPI_ISL_*.fasta
for fasta in EPI_ISL_*.fasta; do
minimap2 -ax asm5 -t ${t} ${g} ${fasta} > ${fasta}.sam
samtools view -bS ${fasta}.sam > ${fasta}.bam
samtools sort -o ${fasta}.sorted.bam ${fasta}.bam
freebayes -f ${g} -C 1 ${fasta}.sorted.bam > ${fasta}.vcf
vcfleftalign -r ${g} ${fasta}.vcf > ${fasta}.left.vcf
rm ${fasta}.sam ${fasta}.bam ${fasta}.sorted.bam ${fasta}.vcf
done
rm merged.GISAID.fasta*
### fixing VCF files for merge
echo "fixing VCF files for merge"
echo ""
vcf= ls -1 *.fasta.left.vcf; for vcf in *.fasta.left.vcf; do sed -i "s|0/0|1/1|"g ${vcf}; done

### Renaming files in bash
echo "Renaming files in bash"
echo ""
for filename in *.fasta.left.vcf; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fasta.left.vcf/.vcf/g')";  done
gzip *.fasta 
cd ..
echo ""
echo "#######################################################"
echo "All done. Variants were called from aligned FASTA files"
echo "#######################################################"
echo ""
