# SARS-CoV-2-freebayes
Analysis of SARS-CoV-2 genome variants collected with freebayes variant caller.

# Requirements: 

### 1) Installing minimap2 aligner (for install details, please see: https://github.com/lh3/minimap2)
```
### Installing minimap2
git clone https://github.com/lh3/minimap2
# build minimap2
cd minimap2 && make
# with sudo privileges
sudo cp minimap2 /usr/local/bin/
```

### 2) Installing fastp: An ultra-fast all-in-one FASTQ preprocessor (for details, please see: https://github.com/OpenGene/fastp)
```
git clone https://github.com/OpenGene/fastp.git
# build fastp
cd fastp
make
# with sudo privileges
sudo cp fastp /usr/local/bin/
```

### 3) Obtaining and installing Freebayes:
Clone Freebayes in a specific directory: 
```
git clone --recursive git://github.com/ekg/freebayes.git

# If this line not work, try:

git config --global url.https://github.com/.insteadOf git://github.com/
git clone --recursive git://github.com/ekg/freebayes.git

#Enter Freebayes directory and make:
cd freebayes
make
sudo make install
sudo cp scripts/* /usr/local/bin/

#To check installation, type in terminal:
freebayes
bamleftalign
```

### 4) Installing vcflib
```
git config --global url.https://github.com/.insteadOf git://github.com/
git clone --recursive git://github.com/vcflib/vcflib.git

#Enter vcflib directory and make
cd vcflib
make   # Needs CMake compiler, with sudo privileges do: sudo apt-get install cmake
cp scripts/* /usr/local/bin/
cp bin/* /usr/local/bin/
```

### 5) Obtaining and installing up-to-date SAMtools with htslib (version >= 1.9)
(Old samtools version can also work). Users need to install version up to date of these three packages. Users can first install htslib v1.9 and then samtools with bcftools v1.9, respectively. For downloading these packages, see http://www.htslib.org/download/). The latter can be accomplished by downloading the three packages, decompressing it, and doing the following:
```
wget https://github.com/samtools/htslib/releases/download/1.10.2/htslib-1.10.2.tar.bz2
bzip2 -d htslib-1.10.2.tar.bz2
tar -xvf htslib-1.10.2.tar
rm htslib-1.10.2.tar
cd htslib-1.10.2    # and similarly for samtools
sudo ./configure --prefix=/usr/local/bin
sudo make
sudo make install
# this step is only for samtools
sudo cp samtools /usr/local/bin/

# Similarly as htslib, samtools and bcftools can be downloaded as follows:

wget https://github.com/samtools/samtools/releases/download/1.10/samtools-1.10.tar.bz2
wget https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2
```

Then in a terminal type
>samtools

to check 1.10 version (using htslib v1.10)

### 6) Obtaining SRA toolkit from ncbi (for downloading reads from SRA archive).
```
### Installing SRA toolkit from ncbi
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.9.6/sratoolkit.2.9.6-ubuntu64.tar.gz
gunzip sratoolkit.2.9.6-ubuntu64.tar.gz
tar -xvf sratoolkit.2.9.6-ubuntu64.tar
cp sratoolkit.2.9.6-ubuntu64/bin/fastq-dump /usr/local/bin/
cp sratoolkit.2.9.6-ubuntu64/bin/prefetch /usr/local/bin/
```

### 7) Installing jaqcuard
For information, please see :https://jacquard.readthedocs.io/en/v0.42/installation.html
```
pip install jacquard
```

### (8) Installing seqkit
SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation (https://bioinf.shenwei.me/seqkit/) can be installed from repository as follows:
```
wget https://github.com/shenwei356/seqkit/releases/download/v0.12.1/seqkit_linux_386.tar.gz
gunzip seqkit_linux_386.tar.gz
tar -xvf seqkit_linux_386.tar
sudo cp seqkit /usr/local/bin/
```

# Colecting Variants (Sequence Read Archive datasets)

In order to obtain SARS-CoV-2 variants (viral frequency >= 0.5) users need to provide:  

- Sequence read archive accessions of each datasets (SRR prefix list, in tabular format or txt format. As example, see: SARS-CoV-2_curated_list_17_07_2020.tabular, provided in this repository)
- SARS-CoV-2 reference in fasta format (covid19-refseq.fasta, provided in this repository)
- number of threads for calculations 

Execution (from scratch): 
```
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git
cd SARS-CoV-2-freebayes
samtools faidx covid19-refseq.fasta
chmod 755 SARS-CoV-2* covid19-refseq.fasta*
./SARS-CoV-2-freebayes.sh SRA_list covid19-refseq.fasta Threads
```
This execution will:

- Download SRA datasets from the provided list, convert to fastq, trim adaptors and gzip reads, for each line of the provided list
- Align trimmed reads against SARS-CoV-2 reference genome (NC_045512.2) by using minimap2
- Call variants (viral frequency >= 0.5) by using freebayes as frequency-based pooled caller
- Merge all variants in a single VCF file by using jacquard. This VCF file also contains viral frequencies in the AF field (AF=AO/AO+RO). See AO and RO fields for alternative and reference allele counts.

### Example 

We provided SARS-CoV-2_curated_list_17_07_2020.tabular, containing a curated list of 16586 SARS-CoV-2 worldwide datasets until July 17, 2020. We also provided curated lists in txt format by continent (see July_28_2020_*.txt files). As an example, to collect variants from July_28_2020_North_America.txt datasets using 30 threads:

```
./SARS-CoV-2-freebayes.sh July_28_2020_North_America.txt covid19-refseq.fasta 30
```
will collect variants (VF>=0.5) in each Sample. To change VF, edit F value in line 138 of SARS-CoV-2-freebayes.sh script.


# Collecting variants per genome: GISAID datasets
- This operation will split every fasta file from a given collection of GISAID SARS-CoV-2 genomes, call variants per genome and will merge variants in a single VCF containing all sample names. The fasta collections from GISAID are the following:
-  gisaid_Africa_08_03_2020.fasta
-  gisaid_Asia_08_03_2020.fasta
-  gisaid_Europe_08_03_2020.fasta
-  gisaid_North_America_08_03_2020.fasta
-  gisaid_Oceania_08_03_2020.fasta
-  gisaid_South_America_08_03_2020.fasta
-  gisaid.merge.fasta (by doing cat gisaid_*_08_03_2020.fasta > gisaid.merge.fasta)

all these fasta files are available here: 

As an example for gisaid.merge.fasta (containing all genomes per geographical region) place gisaid.merge.fasta and covid19-refseq.fasta in a folder and do:

```

### Split fasta files

ulimit -s 80000
sed -i 's/ /-/'g gisaid.merge.fasta
sed -i "s|hCoV-19/.*./2020||"g gisaid.merge.fasta
sed -i "s|hCoV-19/.*./2019||"g gisaid.merge.fasta
sed -i 's/|/\t/'g gisaid.merge.fasta
sed -i 's/>\t/>/'g gisaid.merge.fasta
seqkit fx2tab gisaid.merge.fasta > gisaid.merge.tabular
awk '{print $1"\t"$3}' gisaid.merge.tabular > gisaid.merge.tab && rm gisaid.merge.tabular
seqkit tab2fx gisaid.merge.tab > gisaid.merge.fasta && rm gisaid.merge.tab
seqkit split --by-id gisaid.merge.fasta
cd gisaid.merge.fasta.split/
for filename in *.fasta; do mv "./$filename" "./$(echo "$filename" | sed -e 's/gisaid.merge.id_//g')";  done
cd ..
cp ./gisaid.merge.fasta.split/EPI*.fasta ./
rm -r -f gisaid.merge.fasta.split gisaid.merge.fasta


### Align fasta files to reference (covid19-refseq.fasta, provided in this repository) and call variants with freebayes (option C 1)

samtools faidx covid19-refseq.fasta
fasta= ls -1 *.fasta
for fasta in *.fasta; do
minimap2 -ax asm5 -t 50 covid19-refseq.fasta ${fasta} > ${fasta}.sam
samtools view -bS ${fasta}.sam > ${fasta}.bam
samtools sort -o ${fasta}.sorted.bam ${fasta}.bam
freebayes -f covid19-refseq.fasta -C 1 ${fasta}.sorted.bam > ${fasta}.vcf
vcfleftalign -r /home/user/MITACS/July_12_2020/SARS-CoV-2_illumina_analysis/1/covid19-refseq.fasta ${fasta}.vcf > ${fasta}.left.vcf
rm ${fasta}.sam ${fasta}.bam ${fasta}.sorted.bam ${fasta}.vcf
done


### Merge of Variants

# fixing VCF files for merge
ulimit -s 80000 && vcf= ls -1 *.fasta.left.vcf; for vcf in *.fasta.left.vcf; do sed -i "s|0/0|1/1|"g ${vcf}; done

# Renaming files in bash
for filename in *.fasta.left.vcf; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.fasta.left.vcf/.vcf/g')";  done

# Calculating Number of Variants per genome
ulimit -s 80000
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

# Merge VCFs using jacquard
ulimit -n 1000000 && jacquard merge --include_all ./ merged.GISAID.vcf

# Left only genotypes in merged VCF
vcfkeepgeno merged.GISAID.vcf GT > merged.GISAID.GT.vcf

# Split variants and header from merged.GT.vcf
grep "#" merged.GISAID.GT.vcf > header
grep -v "#" merged.GISAID.GT.vcf > variants.vcf

sed -i 's|1/1|1|'g variants.vcf   
sed -i 's|0/1|1|'g variants.vcf   
sed -i 's|1/0|1|'g variants.vcf  
sed -i 's/[.]/0/'g variants.vcf   # convert point to zeros 

# Reconstitute vcf file
cat header variants.vcf > merged.GISAID.fixed.vcf
rm header variants.vcf
sed -i 's/NC_04551202/NC_045512.2/'g merged.GISAID.fixed.vcf

# left-align vcf file and fix names
vcfleftalign -r /home/user/MITACS/July_12_2020/SARS-CoV-2_illumina_analysis/1/covid19-refseq.fasta merged.GISAID.fixed.vcf > merged.GISAID.left.vcf
sed -i 's/|unknown//'g merged.GISAID.left.vcf

# calculate AF
vcffixup merged.GISAID.left.vcf > merged.GISAID.AF.vcf
rm merged.GISAID.fixed.vcf merged.GISAID.left.vcf
gzip merged.GISAID.vcf
ulimit -s 80000
gzip *.fasta

# Filter variants by Viral Frequency: 0.0499 (5%)
vcffilter -f "AF > 0.0499" merged.GISAID.AF.vcf > merged.GISAID.AF_5%.vcf
grep -v "##" merged.GISAID.AF_5%.vcf > merged.GISAID.AF_5%.table
```

Users can do the same for each fasta collection file to collect aggregated variants per region (merged.GISAID.AF.vcf) and aggregated variants filtered with Viral Frequency > 5% (merged.GISAID.AF_5%.vcf). 


# Collecting variants per Protein (SnpEff-classified GISAID merged variants: working with SnpEff-eff_merged.GISAID.vcf file)
In a given folder, place SnpEff-eff_merged.GISAID.vcf (available for download here: https://usegalaxy.org/u/carlosfarkas/h/sars-cov-2-variants-gisaid-august-03-2020) and do the following: 
```
grep "#" -v SnpEff-eff_merged.GISAID.vcf > variants.vcf
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
cd ..

### Compute the frequencies of synonymous, missense, nonsense and frameshift variants (Overall)  ###
grep "missense_variant" SnpEff-eff_merged.GISAID.vcf > missense_variant.GISAID.SnpEff
grep "stop_gained" SnpEff-eff_merged.GISAID.vcf > stop_gained.GISAID.SnpEff
grep "synonymous_variant" SnpEff-eff_merged.GISAID.vcf > synonymous_variant.GISAID.SnpEff
grep "frameshift_variant" SnpEff-eff_merged.GISAID.vcf > frameshift_variant.GISAID.SnpEff
sed -i 's/;/\t/'g missense_variant.GISAID.SnpEff
sed -i 's/;/\t/'g stop_gained.GISAID.SnpEff
sed -i 's/;/\t/'g synonymous_variant.GISAID.SnpEff
sed -i 's/;/\t/'g frameshift_variant.GISAID.SnpEff
awk '{print $8}' missense_variant.GISAID.SnpEff > missense_variant.GISAID.counts
awk '{print $8}' stop_gained.GISAID.SnpEff > stop_gained.GISAID.counts
awk '{print $8}' synonymous_variant.GISAID.SnpEff > synonymous_variant.GISAID.counts
awk '{print $8}' frameshift_variant.GISAID.SnpEff > frameshift_variant.GISAID.counts
sed -i 's/AC=//'g missense_variant.GISAID.counts
sed -i 's/AC=//'g stop_gained.GISAID.counts
sed -i 's/AC=//'g synonymous_variant.GISAID.counts
sed -i 's/AC=//'g frameshift_variant.GISAID.counts
gzip SnpEff-eff_merged.GISAID.vcf
```


# Ka-Ks calculation per geographical region
In a folder (i.e. Ka-Ks) place covid19-refseq.fasta and the fasta collections per geographical region. Then do:
```
mkdir Ka-Ks
cd ./Ka-Ks/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz
gunzip *
gffread GCF_009858895.2_ASM985889v3_genomic.gff -T -o SARS-CoV-2.gtf
git clone https://github.com/chasewnelson/SNPGenie.git


################
### Analysis ###
################

### Africa
minimap2 -ax asm5 -t 50 covid19-refseq.fasta gisaid_Africa_08_03_2020.fasta > Africa_alignment.sam
samtools view -bS Africa_alignment.sam > Africa_alignment.bam
samtools sort -o Africa_alignment.sorted.bam Africa_alignment.bam
samtools index Africa_alignment.sorted.bam
freebayes -f covid19-refseq.fasta -F 0.001 Africa_alignment.sorted.bam > Africa_alignment.vcf
vcfleftalign -r covid19-refseq.fasta Africa_alignment.vcf > Africa.left.vcf
perl ./SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=Africa.left.vcf --fastafile=covid19-refseq.fasta --gtffile=SARS-CoV-2.gtf --outdir=Africa


### Asia
minimap2 -ax asm5 -t 50 covid19-refseq.fasta gisaid_Asia_08_03_2020.fasta > Asia_alignment.sam
samtools view -bS Asia_alignment.sam > Asia_alignment.bam
samtools sort -o Asia_alignment.sorted.bam Asia_alignment.bam
samtools index Asia_alignment.sorted.bam
freebayes -f covid19-refseq.fasta -F 0.001 Asia_alignment.sorted.bam > Asia_alignment.vcf
vcfleftalign -r covid19-refseq.fasta Asia_alignment.vcf > Asia.left.vcf
perl ./SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=Asia.left.vcf --fastafile=covid19-refseq.fasta --gtffile=SARS-CoV-2.gtf --outdir=Asia


### Europe
minimap2 -ax asm5 -t 50 covid19-refseq.fasta gisaid_Europe_08_03_2020.fasta > Europe_alignment.sam
samtools view -bS Europe_alignment.sam > Europe_alignment.bam
samtools sort -o Europe_alignment.sorted.bam Europe_alignment.bam
samtools index Europe_alignment.sorted.bam
# Split Variant calling for large datasets (memory issues)
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:1-5000 Europe_alignment.sorted.bam > Europe_alignment_1_5000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:5001-6000 Europe_alignment.sorted.bam > Europe_alignment_5001_6000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:6001-7000 Europe_alignment.sorted.bam > Europe_alignment_6001_7000.vcf #
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:7001-8000 Europe_alignment.sorted.bam > Europe_alignment_7001_8000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:8001-9000 Europe_alignment.sorted.bam > Europe_alignment_8001_9000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:9001-10000 Europe_alignment.sorted.bam > Europe_alignment_9001_10000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:10001-11000 Europe_alignment.sorted.bam > Europe_alignment_10001_11000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:11001-12000 Europe_alignment.sorted.bam > Europe_alignment_11001_12000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:13001-14000 Europe_alignment.sorted.bam > Europe_alignment_13001_14000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:14001-15000 Europe_alignment.sorted.bam > Europe_alignment_14001_15000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:15001-16000 Europe_alignment.sorted.bam > Europe_alignment_15001_16000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:16001-17000 Europe_alignment.sorted.bam > Europe_alignment_16001_17000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:17001-18000 Europe_alignment.sorted.bam > Europe_alignment_17001_18000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:18001-19000 Europe_alignment.sorted.bam > Europe_alignment_18001_19000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:19001-20000 Europe_alignment.sorted.bam > Europe_alignment_19001_20000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:20001-21000 Europe_alignment.sorted.bam > Europe_alignment_20001_21000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:21001-22000 Europe_alignment.sorted.bam > Europe_alignment_21001_22000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:22001-23000 Europe_alignment.sorted.bam > Europe_alignment_22001_23000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:23001-24000 Europe_alignment.sorted.bam > Europe_alignment_23001_24000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:24001-25000 Europe_alignment.sorted.bam > Europe_alignment_24001_25000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:25001-26000 Europe_alignment.sorted.bam > Europe_alignment_25001_26000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:26001-27000 Europe_alignment.sorted.bam > Europe_alignment_26001_27000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:27001-28000 Europe_alignment.sorted.bam > Europe_alignment_27001_28000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:28001-29903 Europe_alignment.sorted.bam > Europe_alignment_28001_29903.vcf

vcfcombine Europe_alignment_* > Europe_alignment.vcf
rm Europe_alignment_*.vcf

vcfleftalign -r covid19-refseq.fasta Europe_alignment.vcf > Europe.left.vcf
perl ./SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=Europe.left.vcf --fastafile=covid19-refseq.fasta --gtffile=SARS-CoV-2.gtf --outdir=Europe


### North_America
minimap2 -ax asm5 -t 50 covid19-refseq.fasta gisaid_North_America_08_03_2020.fasta > North_America_alignment.sam
samtools view -bS North_America_alignment.sam > North_America_alignment.bam
samtools sort -o North_America_alignment.sorted.bam North_America_alignment.bam
samtools index North_America_alignment.sorted.bam
# Split Variant calling for large datasets (memory issues)
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:1-5000 North_America_alignment.sorted.bam > North_America_alignment_1_5000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:5001-6000 North_America_alignment.sorted.bam > North_America_alignment_5001_6000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:6001-7000 North_America_alignment.sorted.bam > North_America_alignment_6001_7000.vcf #
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:7001-8000 North_America_alignment.sorted.bam > North_America_alignment_7001_8000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:8001-9000 North_America_alignment.sorted.bam > North_America_alignment_8001_9000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:9001-10000 North_America_alignment.sorted.bam > North_America_alignment_9001_10000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:10001-11000 North_America_alignment.sorted.bam > North_America_alignment_10001_11000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:11001-12000 North_America_alignment.sorted.bam > North_America_alignment_11001_12000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:13001-14000 North_America_alignment.sorted.bam > North_America_alignment_13001_14000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:14001-15000 North_America_alignment.sorted.bam > North_America_alignment_14001_15000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:15001-16000 North_America_alignment.sorted.bam > North_America_alignment_15001_16000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:16001-17000 North_America_alignment.sorted.bam > North_America_alignment_16001_17000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:17001-18000 North_America_alignment.sorted.bam > North_America_alignment_17001_18000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:18001-19000 North_America_alignment.sorted.bam > North_America_alignment_18001_19000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:19001-20000 North_America_alignment.sorted.bam > North_America_alignment_19001_20000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:20001-21000 North_America_alignment.sorted.bam > North_America_alignment_20001_21000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:21001-22000 North_America_alignment.sorted.bam > North_America_alignment_21001_22000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:22001-23000 North_America_alignment.sorted.bam > North_America_alignment_22001_23000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:23001-24000 North_America_alignment.sorted.bam > North_America_alignment_23001_24000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:24001-25000 North_America_alignment.sorted.bam > North_America_alignment_24001_25000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:25001-26000 North_America_alignment.sorted.bam > North_America_alignment_25001_26000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:26001-27000 North_America_alignment.sorted.bam > North_America_alignment_26001_27000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:27001-28000 North_America_alignment.sorted.bam > North_America_alignment_27001_28000.vcf
freebayes -f covid19-refseq.fasta -F 0.001 --region NC_045512.2:28001-29903 North_America_alignment.sorted.bam > North_America_alignment_28001_29903.vcf

vcfcombine North_America_alignment_* > North_America_alignment.vcf
rm North_America_alignment_*.vcf

vcfleftalign -r covid19-refseq.fasta North_America_alignment.vcf > North_America.left.vcf
perl ./SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=North_America.left.vcf --fastafile=covid19-refseq.fasta --gtffile=SARS-CoV-2.gtf --outdir=North_America


### South_America
minimap2 -ax asm5 -t 50 covid19-refseq.fasta gisaid_South_America_08_03_2020.fasta > South_America_alignment.sam
samtools view -bS South_America_alignment.sam > South_America_alignment.bam
samtools sort -o South_America_alignment.sorted.bam South_America_alignment.bam
samtools index South_America_alignment.sorted.bam
freebayes -f covid19-refseq.fasta -F 0.001 South_America_alignment.sorted.bam > South_America_alignment.vcf
vcfleftalign -r covid19-refseq.fasta South_America_alignment.vcf > South_America.left.vcf
perl ./SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=South_America.left.vcf --fastafile=covid19-refseq.fasta --gtffile=SARS-CoV-2.gtf --outdir=South_America


### Oceania
minimap2 -ax asm5 -t 50 covid19-refseq.fasta gisaid_Oceania_08_03_2020.fasta > Oceania_alignment.sam
samtools view -bS Oceania_alignment.sam > Oceania_alignment.bam
samtools sort -o Oceania_alignment.sorted.bam Oceania_alignment.bam
samtools index Oceania_alignment.sorted.bam
freebayes -f covid19-refseq.fasta -F 0.001 Oceania_alignment.sorted.bam > Oceania_alignment.vcf
vcfleftalign -r covid19-refseq.fasta Oceania_alignment.vcf > Oceania.left.vcf
perl ./SNPGenie/snpgenie.pl --vcfformat=4 --snpreport=Oceania.left.vcf --fastafile=covid19-refseq.fasta --gtffile=SARS-CoV-2.gtf --outdir=Oceania
```
