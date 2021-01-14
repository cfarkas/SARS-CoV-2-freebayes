#!/bin/bash
set -e 

GISAID_fasta=${1}
Reference=${2}
Threads=${3}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This script will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This script will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This script will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [GISAID_fasta] [Reference] [Threads]"
  echo ""
  echo "This script will call variants using freebayes in given GISAID fasta sequences files to obtain major viral variants."
  echo ""
  echo "[GISAID_fasta] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: Full PATH where the SARS-CoV-2 reference genome (in fasta format) is located."
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

### Split fasta files
echo "fixing names in FASTA file"
echo ""
fasta_name=$(echo ${1} | sed "s/.fasta//")
sed -i 's/ /-/'g ${1}
sed -i "s|hCoV-19/.*./2020||"g ${1}
sed -i "s|hCoV-19/.*./2019||"g ${1}
sed -i 's/|/\t/'g ${1}
sed -i 's/>\t/>/'g ${1}
seqkit fx2tab ${1} > merged.GISAID.tabular
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
samtools faidx ${2}
fasta= ls -1 EPI_ISL_*.fasta
for fasta in EPI_ISL_*.fasta; do
minimap2 -ax asm5 -t ${3} ${2} ${fasta} > ${fasta}.sam
samtools view -bS ${fasta}.sam > ${fasta}.bam
samtools sort -o ${fasta}.sorted.bam ${fasta}.bam
freebayes -f ${2} -C 1 ${fasta}.sorted.bam > ${fasta}.vcf
vcfleftalign -r ${2} ${fasta}.vcf > ${fasta}.left.vcf
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