#!/bin/bash
set -e 

SnpEff_vcf=${1}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SnpEff_vcf]"
  echo ""
  echo "This script will parse SnpEff annotated variants per protein and also output synonymous, non-synonymous, frameshift and stop-gained variants"
  echo ""
  echo "[SnpEff_vcf] : GISAID genomes in FASTA format"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SnpEff_vcf]"
  echo ""
  echo "This script will parse SnpEff annotated variants per protein and also output synonymous, non-synonymous, frameshift and stop-ga"
  echo ""
  echo "[SnpEff_vcf] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SnpEff_vcf]"
  echo ""
  echo "This script will parse SnpEff annotated variants per protein and also output synonymous, non-synonymous, frameshift and stop-ga"
  echo ""
  echo "[SnpEff_vcf] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [SnpEff_vcf]"
  echo ""
  echo "This script will parse SnpEff annotated variants per protein and also output synonymous, non-synonymous, frameshift and stop-ga"
  echo ""
  echo "[SnpEff_vcf] : GISAID genomes in FASTA format"
  echo ""
  echo "[Reference]: PATH where the SARS-CoV-2 reference genome (in fasta format) is located. If the genome is located in the working folder, just specify the name."
  echo ""
  echo "[Threads]: Number of CPUs for the task (integer)"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [SnpEff_vcf]"; exit 1; }

if [ $# -ne 1 ]; then
  echo 1>&2 "Usage: ./`basename $0` [SnpEff_vcf]"
  exit 3
fi

echo "Parsing SnpEff annotated file by protein"
echo ""
grep "#" -v ${1} > variants.vcf
awk '{print $1"\t"$2"\t"$4"\t"$5"\t"$8}' variants.vcf > SnpEff.sites
sed -i 's/|/\t/'g SnpEff.sites
sed -i 's/Ala/A/'g SnpEff.sites
sed -i 's/Arg/R/'g SnpEff.sites
sed -i 's/Asn/N/'g SnpEff.sites
sed -i 's/Asp/D/'g SnpEff.sites
sed -i 's/Cys/C/'g SnpEff.sites
sed -i 's/Glu/E/'g SnpEff.sites
sed -i 's/Gln/Q/'g SnpEff.sites
sed -i 's/Gly/G/'g SnpEff.sites
sed -i 's/His/H/'g SnpEff.sites
sed -i 's/His/H/'g SnpEff.sites
sed -i 's/Ile/I/'g SnpEff.sites
sed -i 's/Leu/L/'g SnpEff.sites
sed -i 's/Lys/K/'g SnpEff.sites
sed -i 's/Met/M/'g SnpEff.sites
sed -i 's/Phe/F/'g SnpEff.sites
sed -i 's/Pro/P/'g SnpEff.sites
sed -i 's/Ser/S/'g SnpEff.sites
sed -i 's/Thr/T/'g SnpEff.sites
sed -i 's/Trp/W/'g SnpEff.sites
sed -i 's/Tyr/Y/'g SnpEff.sites
sed -i 's/Val/V/'g SnpEff.sites

rm variants.vcf
grep "protein_coding" SnpEff.sites > SnpEff.coding.sites
mkdir variants_per_protein
grep "GU280_gp04" SnpEff.coding.sites > ./variants_per_protein/E.variants
grep "GU280_gp05" SnpEff.coding.sites > ./variants_per_protein/M.variants
grep "GU280_gp10" SnpEff.coding.sites > ./variants_per_protein/N.variants
grep "GU280_gp11" SnpEff.coding.sites > ./variants_per_protein/ORF10.variants
grep "YP_009725297.1" SnpEff.coding.sites > ./variants_per_protein/leader_protein.variants
grep "YP_009725298.1" SnpEff.coding.sites > ./variants_per_protein/nsp2.variants
grep "YP_009725299.1" SnpEff.coding.sites > ./variants_per_protein/nsp3.variants
grep "YP_009725300.1" SnpEff.coding.sites > ./variants_per_protein/nsp4.variants
grep "YP_009725301.1" SnpEff.coding.sites > ./variants_per_protein/3C-like-proteinase.variants
grep "YP_009725302.1" SnpEff.coding.sites > ./variants_per_protein/nsp6.variants
grep "YP_009725303.1" SnpEff.coding.sites > ./variants_per_protein/nsp7.variants
grep "YP_009725304.1" SnpEff.coding.sites > ./variants_per_protein/nsp8.variants
grep "YP_009725305.1" SnpEff.coding.sites > ./variants_per_protein/nsp9.variants
grep "YP_009725306.1" SnpEff.coding.sites > ./variants_per_protein/nsp10.variants
grep "YP_009725307.1" SnpEff.coding.sites > ./variants_per_protein/RNA-dependent-RNA-polymerase.variants
grep "YP_009725308.1" SnpEff.coding.sites > ./variants_per_protein/helicase.variants
grep "YP_009725309.1" SnpEff.coding.sites > ./variants_per_protein/3-to-5-exonuclease.variants
grep "YP_009725310.1" SnpEff.coding.sites > ./variants_per_protein/endoRNAse.variants
grep "YP_009725311.1" SnpEff.coding.sites > ./variants_per_protein/2-O-ribose-methyltransferase.variants
grep "YP_009725312.1" SnpEff.coding.sites > ./variants_per_protein/nsp11.variants
grep "GU280_gp03" SnpEff.coding.sites > ./variants_per_protein/ORF3a.variants
grep "GU280_gp06" SnpEff.coding.sites > ./variants_per_protein/ORF6.variants
grep "GU280_gp07" SnpEff.coding.sites > ./variants_per_protein/ORF7a.variants
grep "GU280_gp08" SnpEff.coding.sites > ./variants_per_protein/ORF7b.variants
grep "GU280_gp09" SnpEff.coding.sites > ./variants_per_protein/ORF8.variants
grep "GU280_gp02" SnpEff.coding.sites > ./variants_per_protein/S.variants
grep "intergenic_region" SnpEff.sites > ./variants_per_protein/intergenic_region.variants

cd variants_per_protein

{
file= ls -1 *.variants
for file in *.variants; do cat ${file} | wc -l
done
#
} | tee logfile_variants
#
grep ".variants" logfile_variants > features
grep -v ".variants" logfile_variants > variants_per_feature
paste features variants_per_feature > logfile_variants
rm features variants_per_feature
cd ..
echo ""
echo "Computing the frequencies of synonymous, missense, nonsense and frameshift variants (Overall)"
echo ""
grep "missense_variant" ${1} > missense_variant.SnpEff
grep "stop_gained" ${1} > stop_gained.SnpEff
grep "synonymous_variant" ${1} > synonymous_variant.SnpEff
grep "frameshift_variant" ${1} > frameshift_variant.SnpEff
sed -i 's/;/\t/'g missense_variant.SnpEff
sed -i 's/;/\t/'g stop_gained.SnpEff
sed -i 's/;/\t/'g synonymous_variant.SnpEff
sed -i 's/;/\t/'g frameshift_variant.SnpEff
awk '{print $8}' missense_variant.SnpEff > missense_variant.counts
awk '{print $8}' stop_gained.SnpEff > stop_gained.counts
awk '{print $8}' synonymous_variant.SnpEff > synonymous_variant.counts
awk '{print $8}' frameshift_variant.SnpEff > frameshift_variant.counts
sed -i 's/AC=//'g missense_variant.counts
sed -i 's/AC=//'g stop_gained.counts
sed -i 's/AC=//'g synonymous_variant.counts
sed -i 's/AC=//'g frameshift_variant.counts
echo ""
echo "Processing AA changes"
echo ""
grep -v "synonymous_variant" SnpEff.coding.sites > SnpEff.coding.sites.nonsynonymous
awk '{print $5"\t"$8"\t"$15}' SnpEff.coding.sites.nonsynonymous > SnpEff.coding.sites.subset
sed -i s'/;AF=/\t/'g SnpEff.coding.sites.subset
sed -i s'/;AN=/\t/'g SnpEff.coding.sites.subset
awk '{print $2"\t"$4"\t"$5}' SnpEff.coding.sites.subset > SnpEff.AAchanges && rm SnpEff.coding.sites.subset
sed -i s'/p.//'g SnpEff.AAchanges
awk '{print $1*100"\t"$2"\t"$3}' SnpEff.AAchanges > SnpEff.AAchanges.percentages
sort -k1,1nr -k2,2 SnpEff.AAchanges.percentages > SnpEff.AAchanges.percentages.sorted
echo "#########"
echo "All done"
echo "Check variants_per_protein folder, containing the provided VCF file, splitted by SARS-CoV-2 protein"
echo ".SnpEff files contain parsed synonymous, missense, nonsense and frameshift variants"
echo ".counts files contain parsed synonymous, missense, nonsense and frameshift frequency counts"
echo "SnpEff.AAchanges.percentages contains viral frequencies (in %) of missense, framshift and stop-gained changes ocurring in SARS-CoV-2 proteins"
echo "SnpEff.AAchanges.percentages.sorted contains ordered viral frequencies (in %) of missense, framshift and stop-gained changes ocurring in SARS-CoV-2 proteins"
echo "SnpEff.AAchanges contains viral frequencies of missense, framshift and stop-gained changes ocurring in SARS-CoV-2 proteins"
echo "#########"
