#!/bin/bash
set -e
{
pi_tajima=${1}
merged_vcf_file=${2}

if [ "$1" == "-h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [pi_tajima] [merged_vcf_file]"
  echo ""
  echo "This script will invoke R, several R libraries and vcfstats to process merged pi and Tajima's D values per bin."
  echo ""
  echo "[pi_tajima]: Provide binned pi and Tajima's D values, in tabular format. As example: Europe.50.pi.D"
  echo ""
  echo "[merged_vcf_file]: Provide merged VCF file with viral frequency values per variant. As example: merged.GISAID.AF.vcf"
  echo ""
  exit 0
fi

if [ "$1" == "-help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [pi_tajima] [merged_vcf_file]"
  echo ""
  echo "This script will invoke R, several R libraries and vcfstats to process merged pi and Tajima's D values per bin."
  echo ""
  echo "[pi_tajima]: Provide binned pi and Tajima's D values, in tabular format. As example: Europe.50.pi.D"
  echo ""
  echo "[merged_vcf_file]: Provide merged VCF file with viral frequency values per variant. As example: merged.GISAID.AF.vcf"
  echo ""
  exit 0
fi
if [ "$1" == "--h" ]; then
  echo ""
  echo "Usage: ./`basename $0` [pi_tajima] [merged_vcf_file]"
  echo ""
  echo "This script will invoke R, several R libraries and vcfstats to process merged pi and Tajima's D values per bin."
  echo ""
  echo "[pi_tajima]: Provide binned pi and Tajima's D values, in tabular format. As example: Europe.50.pi.D"
  echo ""
  echo "[merged_vcf_file]: Provide merged VCF file with viral frequency values per variant. As example: merged.GISAID.AF.vcf"
  echo ""
  exit 0
fi

if [ "$1" == "--help" ]; then
  echo ""
  echo "Usage: ./`basename $0` [pi_tajima] [merged_vcf_file]"
  echo ""
  echo "This script will invoke R, several R libraries and vcfstats to process merged pi and Tajima's D values per bin."
  echo ""
  echo "[pi_tajima]: Provide binned pi and Tajima's D values, in tabular format. As example: Europe.50.pi.D"
  echo ""
  echo "[merged_vcf_file]: Provide merged VCF file with viral frequency values per variant. As example: merged.GISAID.AF.vcf"
  echo ""
  exit 0
fi

[ $# -eq 0 ] && { echo "Usage: ./`basename $0` [pi_tajima] [merged_vcf_file]"; exit 1; }

if [ $# -ne 2 ]; then
  echo 1>&2 "Usage: ./`basename $0` [pi_tajima] [merged_vcf_file]"
  exit 3
fi

# Obtaining 95% confidence intervals in Tajima's D and parsing provided VCF file
echo "====> Obtaining 95% confidence intervals in Tajima's D and parsing provided VCF file"
echo ""
dir1=$(cd -P -- "$(dirname -- "$0")" && pwd -P)
echo "====> The working directory is the following: ${dir1}/"
echo ""
wget https://raw.githubusercontent.com/cfarkas/SARS-CoV-2-freebayes/master/confidence_interval.R
cp ${1} input.pi.D
grep -v "nan" input.pi.D > input.pi.D.clean && rm input.pi.D
mv input.pi.D.clean input.pi.D
Rscript confidence_interval.R input.pi.D && rm input.pi.D
sed -i 's/V1/bins/'g bins_2.5%_confidence.tab
sed -i 's/V1/bins/'g bins_97.5%_confidence.tab
awk '{print "NC_045512.2""\t"($1-50)"\t"$1}' bins_2.5%_confidence.tab > bins_2.5%_confidence.bed
awk '{print "NC_045512.2""\t"($1-50)"\t"$1}' bins_97.5%_confidence.tab > bins_97.5%_confidence.bed
tail -n +2 bins_2.5%_confidence.bed > bins_2.5%_conf_interval.bed && rm bins_2.5%_confidence.tab bins_2.5%_confidence.bed
tail -n +2 bins_97.5%_confidence.bed > bins_97.5%_conf_interval.bed && rm bins_97.5%_confidence.tab bins_97.5%_confidence.bed
mv bins_2.5%_conf_interval.bed bins_2.5%_confidence.bed
mv bins_97.5%_conf_interval.bed bins_97.5%_confidence.bed
vcftools --vcf ${2} --bed bins_2.5%_confidence.bed --out 2.5_CI_confidence --recode --keep-INFO-all
vcftools --vcf ${2} --bed bins_97.5%_confidence.bed --out 97.5_CI_confidence --recode --keep-INFO-all
echo ""
# making temporal directory and parsing vcf file by protein or feature
echo "====> Making temporal directory and parsing vcf file by protein or feature: outside 2.5% CI"
echo ""
mkdir vcfsplit-temporal-directory
cp 2.5_CI_confidence.recode.vcf ./vcfsplit-temporal-directory/
cp 97.5_CI_confidence.recode.vcf ./vcfsplit-temporal-directory/
cd vcfsplit-temporal-directory
# Working with bins outside 2.5% CI
echo "====> Working with bins outside 2.5% CI"
echo ""
grep "#" -v 2.5_CI_confidence.recode.vcf > variants.vcf
grep "#" 2.5_CI_confidence.recode.vcf > vcfheader
awk '{ if ($2>=1 && $2<=265) { print } }' variants.vcf > five_prime_utr.vcfsplit
awk '{ if ($2>=266 && $2<=805) { print } }' variants.vcf > leader_protein.vcfsplit
awk '{ if ($2>=806 && $2<=2719) { print } }' variants.vcf > nsp2.vcfsplit
awk '{ if ($2>=2720 && $2<=8554) { print } }' variants.vcf > nsp3.vcfsplit
awk '{ if ($2>=8555 && $2<=10054) { print } }' variants.vcf > nsp4.vcfsplit
awk '{ if ($2>=10055 && $2<=10972) { print } }' variants.vcf > 3C_like_proteinase.vcfsplit
awk '{ if ($2>=10973 && $2<=11842) { print } }' variants.vcf > nsp6.vcfsplit
awk '{ if ($2>=11843 && $2<=12091) { print } }' variants.vcf > nsp7.vcfsplit
awk '{ if ($2>=12092 && $2<=12685) { print } }' variants.vcf > nsp8.vcfsplit
awk '{ if ($2>=12686 && $2<=13024) { print } }' variants.vcf > nsp9.vcfsplit
awk '{ if ($2>=13025 && $2<=13441) { print } }' variants.vcf > nsp10.vcfsplit
awk '{ if ($2>=13442 && $2<=16236) { print } }' variants.vcf > RNA_dependent_RNA_polymerase.vcfsplit
awk '{ if ($2>=16237 && $2<=18039) { print } }' variants.vcf > helicase.vcfsplit
awk '{ if ($2>=18040 && $2<=19620) { print } }' variants.vcf > 3_prime-to-5_prime_exonuclease.vcfsplit
awk '{ if ($2>=19621 && $2<=20658) { print } }' variants.vcf > endoRNAse.vcfsplit
awk '{ if ($2>=20659 && $2<=21552) { print } }' variants.vcf > 2_prime-O-ribose_methyltransferase.vcfsplit
awk '{ if ($2>=21563 && $2<=25384) { print } }' variants.vcf > S.vcfsplit
awk '{ if ($2>=25393 && $2<=26220) { print } }' variants.vcf > ORF3a.vcfsplit
awk '{ if ($2>=26245 && $2<=26472) { print } }' variants.vcf > E.vcfsplit
awk '{ if ($2>=26523 && $2<=27191) { print } }' variants.vcf > M.vcfsplit
awk '{ if ($2>=27202 && $2<=27387) { print } }' variants.vcf > ORF6.vcfsplit
awk '{ if ($2>=27394 && $2<=27759) { print } }' variants.vcf > ORF7a.vcfsplit
awk '{ if ($2>=27756 && $2<=27887) { print } }' variants.vcf > ORF7b.vcfsplit
awk '{ if ($2>=27894 && $2<=28259) { print } }' variants.vcf > ORF8.vcfsplit
awk '{ if ($2>=28274 && $2<=29533) { print } }' variants.vcf > N.vcfsplit
awk '{ if ($2>=29558 && $2<=29674) { print } }' variants.vcf > ORF10.vcfsplit
awk '{ if ($2>=29675 && $2<=29903) { print } }' variants.vcf > three_prime_utr.vcfsplit
# generating stats
echo "====> Generating stats: outside 2.5% CI"
echo ""
echo "five_prime_utr" >> variants_per_feature.txt
echo "leader_protein" >> variants_per_feature.txt
echo "nsp2" >> variants_per_feature.txt
echo "nsp3" >> variants_per_feature.txt
echo "nsp4" >> variants_per_feature.txt
echo "3C_like_proteinase" >> variants_per_feature.txt
echo "nsp6" >> variants_per_feature.txt
echo "nsp7" >> variants_per_feature.txt
echo "nsp8" >> variants_per_feature.txt
echo "nsp9" >> variants_per_feature.txt
echo "nsp10" >> variants_per_feature.txt
echo "RNA_dependent_RNA_polymerase" >> variants_per_feature.txt
echo "helicase" >> variants_per_feature.txt
echo "3_prime-to-5_prime_exonuclease" >> variants_per_feature.txt
echo "endoRNAse" >> variants_per_feature.txt
echo "2_prime-O-ribose_methyltransferase" >> variants_per_feature.txt
echo "S" >> variants_per_feature.txt
echo "ORF3a" >> variants_per_feature.txt
echo "E" >> variants_per_feature.txt
echo "M" >> variants_per_feature.txt
echo "ORF6" >> variants_per_feature.txt
echo "ORF7a" >> variants_per_feature.txt
echo "ORF7b" >> variants_per_feature.txt
echo "ORF8" >> variants_per_feature.txt
echo "N" >> variants_per_feature.txt
echo "ORF10" >> variants_per_feature.txt
echo "three_prime_utr" >> variants_per_feature.txt
cat five_prime_utr.vcfsplit | wc -l >> variants_per_feature.txt
cat leader_protein.vcfsplit | wc -l >> variants_per_feature.txt
cat nsp2.vcfsplit | wc -l >> variants_per_feature.txt
cat nsp3.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp4.vcfsplit | wc -l  >> variants_per_feature.txt
cat 3C_like_proteinase.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp6.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp7.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp8.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp9.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp10.vcfsplit | wc -l  >> variants_per_feature.txt
cat RNA_dependent_RNA_polymerase.vcfsplit | wc -l  >> variants_per_feature.txt
cat helicase.vcfsplit | wc -l  >> variants_per_feature.txt
cat 3_prime-to-5_prime_exonuclease.vcfsplit | wc -l  >> variants_per_feature.txt
cat endoRNAse.vcfsplit | wc -l  >> variants_per_feature.txt
cat 2_prime-O-ribose_methyltransferase.vcfsplit | wc -l  >> variants_per_feature.txt
cat S.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF3a.vcfsplit | wc -l  >> variants_per_feature.txt
cat E.vcfsplit | wc -l  >> variants_per_feature.txt
cat M.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF6.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF7a.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF7b.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF8.vcfsplit | wc -l  >> variants_per_feature.txt
cat N.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF10.vcfsplit | wc -l  >> variants_per_feature.txt
cat three_prime_utr.vcfsplit | wc -l  >> variants_per_feature.txt
grep '[Aa-Zz]' variants_per_feature.txt > SARS-CoV-2_proteins
grep -v '[Aa-Zz]' variants_per_feature.txt > variants_per_sample
paste SARS-CoV-2_proteins variants_per_sample > variants_per_feature
rm SARS-CoV-2_proteins variants_per_sample variants_per_feature.txt
mv variants_per_feature variants_per_feature_2.5-CI.tabular
# Replacing chromosome name with gene names 
echo "====> Replacing chromosome name with gene names : outside 2.5% CI"
echo ""
sed -i 's/NC_045512.2/five_prime_utr/'g five_prime_utr.vcfsplit
sed -i 's/NC_045512.2/leader_protein/'g leader_protein.vcfsplit
sed -i 's/NC_045512.2/nsp2/'g nsp2.vcfsplit
sed -i 's/NC_045512.2/nsp3/'g nsp3.vcfsplit
sed -i 's/NC_045512.2/nsp4/'g nsp4.vcfsplit
sed -i 's/NC_045512.2/3C_like_proteinase/'g 3C_like_proteinase.vcfsplit
sed -i 's/NC_045512.2/nsp6/'g nsp6.vcfsplit
sed -i 's/NC_045512.2/nsp7/'g nsp7.vcfsplit
sed -i 's/NC_045512.2/nsp8/'g nsp8.vcfsplit
sed -i 's/NC_045512.2/nsp9/'g nsp9.vcfsplit
sed -i 's/NC_045512.2/nsp10/'g nsp10.vcfsplit
sed -i 's/NC_045512.2/RNA_dependent_RNA_polymerase/'g RNA_dependent_RNA_polymerase.vcfsplit
sed -i 's/NC_045512.2/helicase/'g helicase.vcfsplit
sed -i 's/NC_045512.2/3_prime-to-5_prime_exonuclease/'g 3_prime-to-5_prime_exonuclease.vcfsplit
sed -i 's/NC_045512.2/endoRNAse/'g endoRNAse.vcfsplit
sed -i 's/NC_045512.2/2_prime-O-ribose_methyltransferase/'g 2_prime-O-ribose_methyltransferase.vcfsplit
sed -i 's/NC_045512.2/S/'g S.vcfsplit
sed -i 's/NC_045512.2/ORF3a/'g ORF3a.vcfsplit
sed -i 's/NC_045512.2/E/'g E.vcfsplit
sed -i 's/NC_045512.2/M/'g M.vcfsplit
sed -i 's/NC_045512.2/ORF6/'g ORF6.vcfsplit
sed -i 's/NC_045512.2/ORF7a/'g ORF7a.vcfsplit
sed -i 's/NC_045512.2/ORF7b/'g ORF7b.vcfsplit
sed -i 's/NC_045512.2/ORF8/'g ORF8.vcfsplit
sed -i 's/NC_045512.2/N/'g N.vcfsplit
sed -i 's/NC_045512.2/ORF10/'g ORF10.vcfsplit
sed -i 's/NC_045512.2/three_prime_utr/'g three_prime_utr.vcfsplit
# tweaking 2.5_CI_confidence.recode.vcf file replacing chromosome with proteins-feature names
echo "====> Tweaking 2.5_CI_confidence.recode.vcf file replacing chromosome with proteins-feature names"
echo ""
cat vcfheader five_prime_utr.vcfsplit leader_protein.vcfsplit nsp2.vcfsplit nsp3.vcfsplit nsp4.vcfsplit 3C_like_proteinase.vcfsplit nsp6.vcfsplit nsp7.vcfsplit nsp8.vcfsplit nsp9.vcfsplit nsp10.vcfsplit RNA_dependent_RNA_polymerase.vcfsplit helicase.vcfsplit 3_prime-to-5_prime_exonuclease.vcfsplit endoRNAse.vcfsplit 2_prime-O-ribose_methyltransferase.vcfsplit S.vcfsplit ORF3a.vcfsplit E.vcfsplit M.vcfsplit ORF6.vcfsplit ORF7a.vcfsplit ORF7b.vcfsplit ORF8.vcfsplit N.vcfsplit ORF10.vcfsplit three_prime_utr.vcfsplit > 2.5_CI_confidence.tweaked.vcf   
rm *vcfsplit variants.vcf vcfheader
# Working with bins outside 97.5% CI
echo "====> Working with bins outside 97.5% CI"
echo ""
grep "#" -v 97.5_CI_confidence.recode.vcf > variants.vcf
grep "#" 97.5_CI_confidence.recode.vcf > vcfheader
awk '{ if ($2>=1 && $2<=265) { print } }' variants.vcf > five_prime_utr.vcfsplit
awk '{ if ($2>=266 && $2<=805) { print } }' variants.vcf > leader_protein.vcfsplit
awk '{ if ($2>=806 && $2<=2719) { print } }' variants.vcf > nsp2.vcfsplit
awk '{ if ($2>=2720 && $2<=8554) { print } }' variants.vcf > nsp3.vcfsplit
awk '{ if ($2>=8555 && $2<=10054) { print } }' variants.vcf > nsp4.vcfsplit
awk '{ if ($2>=10055 && $2<=10972) { print } }' variants.vcf > 3C_like_proteinase.vcfsplit
awk '{ if ($2>=10973 && $2<=11842) { print } }' variants.vcf > nsp6.vcfsplit
awk '{ if ($2>=11843 && $2<=12091) { print } }' variants.vcf > nsp7.vcfsplit
awk '{ if ($2>=12092 && $2<=12685) { print } }' variants.vcf > nsp8.vcfsplit
awk '{ if ($2>=12686 && $2<=13024) { print } }' variants.vcf > nsp9.vcfsplit
awk '{ if ($2>=13025 && $2<=13441) { print } }' variants.vcf > nsp10.vcfsplit
awk '{ if ($2>=13442 && $2<=16236) { print } }' variants.vcf > RNA_dependent_RNA_polymerase.vcfsplit
awk '{ if ($2>=16237 && $2<=18039) { print } }' variants.vcf > helicase.vcfsplit
awk '{ if ($2>=18040 && $2<=19620) { print } }' variants.vcf > 3_prime-to-5_prime_exonuclease.vcfsplit
awk '{ if ($2>=19621 && $2<=20658) { print } }' variants.vcf > endoRNAse.vcfsplit
awk '{ if ($2>=20659 && $2<=21552) { print } }' variants.vcf > 2_prime-O-ribose_methyltransferase.vcfsplit
awk '{ if ($2>=21563 && $2<=25384) { print } }' variants.vcf > S.vcfsplit
awk '{ if ($2>=25393 && $2<=26220) { print } }' variants.vcf > ORF3a.vcfsplit
awk '{ if ($2>=26245 && $2<=26472) { print } }' variants.vcf > E.vcfsplit
awk '{ if ($2>=26523 && $2<=27191) { print } }' variants.vcf > M.vcfsplit
awk '{ if ($2>=27202 && $2<=27387) { print } }' variants.vcf > ORF6.vcfsplit
awk '{ if ($2>=27394 && $2<=27759) { print } }' variants.vcf > ORF7a.vcfsplit
awk '{ if ($2>=27756 && $2<=27887) { print } }' variants.vcf > ORF7b.vcfsplit
awk '{ if ($2>=27894 && $2<=28259) { print } }' variants.vcf > ORF8.vcfsplit
awk '{ if ($2>=28274 && $2<=29533) { print } }' variants.vcf > N.vcfsplit
awk '{ if ($2>=29558 && $2<=29674) { print } }' variants.vcf > ORF10.vcfsplit
awk '{ if ($2>=29675 && $2<=29903) { print } }' variants.vcf > three_prime_utr.vcfsplit
# generating stats
echo "====> Generating stats: outside 97.5% CI"
echo ""
echo "five_prime_utr" >> variants_per_feature.txt
echo "leader_protein" >> variants_per_feature.txt
echo "nsp2" >> variants_per_feature.txt
echo "nsp3" >> variants_per_feature.txt
echo "nsp4" >> variants_per_feature.txt
echo "3C_like_proteinase" >> variants_per_feature.txt
echo "nsp6" >> variants_per_feature.txt
echo "nsp7" >> variants_per_feature.txt
echo "nsp8" >> variants_per_feature.txt
echo "nsp9" >> variants_per_feature.txt
echo "nsp10" >> variants_per_feature.txt
echo "RNA_dependent_RNA_polymerase" >> variants_per_feature.txt
echo "helicase" >> variants_per_feature.txt
echo "3_prime-to-5_prime_exonuclease" >> variants_per_feature.txt
echo "endoRNAse" >> variants_per_feature.txt
echo "2_prime-O-ribose_methyltransferase" >> variants_per_feature.txt
echo "S" >> variants_per_feature.txt
echo "ORF3a" >> variants_per_feature.txt
echo "E" >> variants_per_feature.txt
echo "M" >> variants_per_feature.txt
echo "ORF6" >> variants_per_feature.txt
echo "ORF7a" >> variants_per_feature.txt
echo "ORF7b" >> variants_per_feature.txt
echo "ORF8" >> variants_per_feature.txt
echo "N" >> variants_per_feature.txt
echo "ORF10" >> variants_per_feature.txt
echo "three_prime_utr" >> variants_per_feature.txt
cat five_prime_utr.vcfsplit | wc -l >> variants_per_feature.txt
cat leader_protein.vcfsplit | wc -l >> variants_per_feature.txt
cat nsp2.vcfsplit | wc -l >> variants_per_feature.txt
cat nsp3.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp4.vcfsplit | wc -l  >> variants_per_feature.txt
cat 3C_like_proteinase.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp6.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp7.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp8.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp9.vcfsplit | wc -l  >> variants_per_feature.txt
cat nsp10.vcfsplit | wc -l  >> variants_per_feature.txt
cat RNA_dependent_RNA_polymerase.vcfsplit | wc -l  >> variants_per_feature.txt
cat helicase.vcfsplit | wc -l  >> variants_per_feature.txt
cat 3_prime-to-5_prime_exonuclease.vcfsplit | wc -l  >> variants_per_feature.txt
cat endoRNAse.vcfsplit | wc -l  >> variants_per_feature.txt
cat 2_prime-O-ribose_methyltransferase.vcfsplit | wc -l  >> variants_per_feature.txt
cat S.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF3a.vcfsplit | wc -l  >> variants_per_feature.txt
cat E.vcfsplit | wc -l  >> variants_per_feature.txt
cat M.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF6.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF7a.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF7b.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF8.vcfsplit | wc -l  >> variants_per_feature.txt
cat N.vcfsplit | wc -l  >> variants_per_feature.txt
cat ORF10.vcfsplit | wc -l  >> variants_per_feature.txt
cat three_prime_utr.vcfsplit | wc -l  >> variants_per_feature.txt
grep '[Aa-Zz]' variants_per_feature.txt > SARS-CoV-2_proteins
grep -v '[Aa-Zz]' variants_per_feature.txt > variants_per_sample
paste SARS-CoV-2_proteins variants_per_sample > variants_per_feature
rm SARS-CoV-2_proteins variants_per_sample variants_per_feature.txt
mv variants_per_feature variants_per_feature_97.5-CI.tabular
# Replacing chromosome name with gene names 
echo "====> Replacing chromosome name with gene names: outside 97.5% CI"
echo ""
sed -i 's/NC_045512.2/five_prime_utr/'g five_prime_utr.vcfsplit
sed -i 's/NC_045512.2/leader_protein/'g leader_protein.vcfsplit
sed -i 's/NC_045512.2/nsp2/'g nsp2.vcfsplit
sed -i 's/NC_045512.2/nsp3/'g nsp3.vcfsplit
sed -i 's/NC_045512.2/nsp4/'g nsp4.vcfsplit
sed -i 's/NC_045512.2/3C_like_proteinase/'g 3C_like_proteinase.vcfsplit
sed -i 's/NC_045512.2/nsp6/'g nsp6.vcfsplit
sed -i 's/NC_045512.2/nsp7/'g nsp7.vcfsplit
sed -i 's/NC_045512.2/nsp8/'g nsp8.vcfsplit
sed -i 's/NC_045512.2/nsp9/'g nsp9.vcfsplit
sed -i 's/NC_045512.2/nsp10/'g nsp10.vcfsplit
sed -i 's/NC_045512.2/RNA_dependent_RNA_polymerase/'g RNA_dependent_RNA_polymerase.vcfsplit
sed -i 's/NC_045512.2/helicase/'g helicase.vcfsplit
sed -i 's/NC_045512.2/3_prime-to-5_prime_exonuclease/'g 3_prime-to-5_prime_exonuclease.vcfsplit
sed -i 's/NC_045512.2/endoRNAse/'g endoRNAse.vcfsplit
sed -i 's/NC_045512.2/2_prime-O-ribose_methyltransferase/'g 2_prime-O-ribose_methyltransferase.vcfsplit
sed -i 's/NC_045512.2/S/'g S.vcfsplit
sed -i 's/NC_045512.2/ORF3a/'g ORF3a.vcfsplit
sed -i 's/NC_045512.2/E/'g E.vcfsplit
sed -i 's/NC_045512.2/M/'g M.vcfsplit
sed -i 's/NC_045512.2/ORF6/'g ORF6.vcfsplit
sed -i 's/NC_045512.2/ORF7a/'g ORF7a.vcfsplit
sed -i 's/NC_045512.2/ORF7b/'g ORF7b.vcfsplit
sed -i 's/NC_045512.2/ORF8/'g ORF8.vcfsplit
sed -i 's/NC_045512.2/N/'g N.vcfsplit
sed -i 's/NC_045512.2/ORF10/'g ORF10.vcfsplit
sed -i 's/NC_045512.2/three_prime_utr/'g three_prime_utr.vcfsplit
# tweaking 97.5_CI_confidence.recode.vcf file replacing chromosome with proteins-feature names
echo "====> Tweaking 97.5_CI_confidence.recode.vcf file replacing chromosome with proteins-feature names"
echo ""
cat vcfheader five_prime_utr.vcfsplit leader_protein.vcfsplit nsp2.vcfsplit nsp3.vcfsplit nsp4.vcfsplit 3C_like_proteinase.vcfsplit nsp6.vcfsplit nsp7.vcfsplit nsp8.vcfsplit nsp9.vcfsplit nsp10.vcfsplit RNA_dependent_RNA_polymerase.vcfsplit helicase.vcfsplit 3_prime-to-5_prime_exonuclease.vcfsplit endoRNAse.vcfsplit 2_prime-O-ribose_methyltransferase.vcfsplit S.vcfsplit ORF3a.vcfsplit E.vcfsplit M.vcfsplit ORF6.vcfsplit ORF7a.vcfsplit ORF7b.vcfsplit ORF8.vcfsplit N.vcfsplit ORF10.vcfsplit three_prime_utr.vcfsplit > 97.5_CI_confidence.tweaked.vcf   
rm *vcfsplit variants.vcf vcfheader
cd ..
cp ./vcfsplit-temporal-directory/2.5_CI_confidence.tweaked.vcf ./
cp ./vcfsplit-temporal-directory/97.5_CI_confidence.tweaked.vcf ./
cp ./vcfsplit-temporal-directory/variants_per_feature_2.5-CI.tabular ./
cp ./vcfsplit-temporal-directory/variants_per_feature_97.5-CI.tabular ./
rm -rf vcfsplit-temporal-directory
# invoking vcfstats on tweaked vcf files
echo "====> Invoking vcfstats on tweaked vcf files"
echo ""
echo "====> Working in postprocessing_pi_D_output_files"
echo ""
mkdir postprocessing_pi_D_output_files
vcfstats --vcf 2.5_CI_confidence.tweaked.vcf \
	 --outdir postprocessing_pi_D_output_files/ \
	 --formula 'COUNT(1) ~ CONTIG' \
	 --title 'Number of variants per protein or feature 2.5 CI'

vcfstats --vcf 2.5_CI_confidence.tweaked.vcf \
	 --outdir postprocessing_pi_D_output_files/ \
	 --formula 'AAF ~ CONTIG' \
	 --title 'Viral frequency of variants per protein or feature 2.5 CI' \
	 --figtype boxplot

vcfstats --vcf 97.5_CI_confidence.tweaked.vcf \
	 --outdir postprocessing_pi_D_output_files/ \
	 --formula 'COUNT(1) ~ CONTIG' \
	 --title 'Number of variants per protein or feature 97.5 CI'

vcfstats --vcf 97.5_CI_confidence.tweaked.vcf \
	 --outdir postprocessing_pi_D_output_files/ \
	 --formula 'AAF ~ CONTIG' \
	 --title 'Viral frequency of variants per protein or feature 97.5 CI' \
	 --figtype boxplot
# cleaning up
echo ""
echo "====> Cleaning up .."
echo ""
gzip 2.5_CI_confidence.tweaked.vcf 2.5_CI_confidence.recode.vcf 97.5_CI_confidence.tweaked.vcf 97.5_CI_confidence.recode.vcf
mv 2.5_CI_confidence.tweaked.vcf.gz 2.5_CI_confidence.recode.vcf.gz 97.5_CI_confidence.tweaked.vcf.gz 97.5_CI_confidence.recode.vcf.gz variants_per_feature_2.5-CI.tabular variants_per_feature_97.5-CI.tabular bins_2.5%_confidence.bed bins_97.5%_confidence.bed pi_tajima.pdf ./postprocessing_pi_D_output_files/
#
} | tee logfile_tajima_d_postprocessing
#
echo "##################"
echo "##################"
echo ""
echo "All done. Check postprocessing_pi_D_output_files folder containing output files"
echo ""
echo "bins_2.5%_confidence.bed contains bins with Tajima's D values out of 2.5% CI. Check pi_tajima.pdf plot"
echo ""
echo "bins_97.5%_confidence.bed contains bins with Tajima's D values out of 97.5% CI. Check pi_tajima.pdf plot"
echo ""
echo "2.5_CI_confidence.recode.vcf.gz contain sites in Tajima's D values out of 2.5% CI"
echo ""
echo "97.5_CI_confidence.recode.vcf.gz contain sites in Tajima's D values out of 97.5% CI"
echo ""
echo "A tweaked version of the latter files conatining protein names instead of chromosome name are called 2.5 and 97.5_CI_confidence.tweaked.vcf.gz, respectively"
echo ""
echo "variants_per_feature_2.5-CI.tabular contain variants falling in bins with Tajima's D values out of 2.5% CI"
echo ""
echo "variants_per_feature_97.5-CI.tabular contain variants falling in bins with Tajima's D values out of 97.5% CI"
echo ""
echo "##################"
echo "##################"
