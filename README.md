# SARS-CoV-2-freebayes
Analysis of SARS-CoV-2 genome variants collected with freebayes variant caller.

## Requirements: 

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

### 8) Installing pyfasta
For information, please see : https://pypi.org/project/pyfasta/
```
pip install pyfasta
```

### 9) Installing seqkit
SeqKit - a cross-platform and ultrafast toolkit for FASTA/Q file manipulation (https://bioinf.shenwei.me/seqkit/) can be installed from repository as follows:
```
wget https://github.com/shenwei356/seqkit/releases/download/v0.12.1/seqkit_linux_386.tar.gz
gunzip seqkit_linux_386.tar.gz
tar -xvf seqkit_linux_386.tar
sudo cp seqkit /usr/local/bin/
```

### 10) Installing BEDOPS
For information, please see: https://bedops.readthedocs.io/en/latest/content/installation.html#linux
```
wget https://github.com/bedops/bedops/releases/download/v2.4.39/bedops_linux_x86_64-v2.4.39.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2
cp bin/* /usr/local/bin/
```

### 11) Installing inStrain 
Tool for highly accurate genome comparisons, analysis of coverage, microdiversity, linkage and sensitive SNP detection. For information, please see: https://instrain.readthedocs.io/en/latest/

```
pip install instrain
```

### 12) Obtaining and Installing VCFtools
We employed the version of Julien Y. Dutheil https://github.com/jydu/vcftools that includes the --haploid flag. Safe install can be achieved with root (sudo -i) as follows: 

```
### Installing vcftools 
sudo -i
cd /home/user/ # go to home directory or another directory
sudo apt-get install libz-dev zlib1g-dev  # zlib requirements in Ubuntu
git clone https://github.com/jydu/vcftools.git
cd vcftools/
./autogen.sh
export PERL5LIB=/path/to/your/vcftools-directory/src/perl/ 
./configure
make
make install
exit
```

## I) Colecting Variants (Sequence Read Archive datasets)

In order to obtain SARS-CoV-2 variants (viral frequency >= 0.5) users need to provide:  

- Sequence read archive accessions of each datasets (SRR prefix list, in tabular format or txt format). As example, see SRA_Accessions_Aug_03_2020.tabular, for curated SARS-CoV-2 accessions from SRA until August 03, 2020, provided in this repository.
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

We provided SRA_Accessions_Aug_03_2020.tabular, containing a curated list of 17560 SARS-CoV-2 worldwide datasets until August 03, 2020. We also provided curated lists in txt format by continent (see July_28_2020_*.txt files). As an example, we will collect variants from July_28_2020_North_America.txt datasets using 30 threads. In Ubuntu it is recommended to increase open file limit and stack size accordingly at the number of genomes to process, otherwise these steps will crush. see README_ulimit for details in Ubuntu. 

```
ulimit -n 1000000 && ulimit -s 299999  # increase stack size and open file limit, see README_ulimit for details.
./SARS-CoV-2-NGS-freebayes.sh July_28_2020_North_America.txt covid19-refseq.fasta 30
```
will collect variants (VF>=0.5) in each Sample. To change VF, edit F value in line 139 of SARS-CoV-2-NGS-freebayes.sh script.


## II) Collecting Variants (GISAID FASTA genomes)
- By using freebayes we can call variants from GISAID SARS-CoV-2 FASTA genomes and merge these variants in a single VCF containing all sample names. The fasta collections from GISAID are the following:
-  gisaid_Africa_08_03_2020.fasta
-  gisaid_Asia_08_03_2020.fasta
-  gisaid_Europe_08_03_2020.fasta
-  gisaid_North_America_08_03_2020.fasta
-  gisaid_Oceania_08_03_2020.fasta
-  gisaid_South_America_08_03_2020.fasta
-  merged.GISAID.fasta.gz (merge)

all these fasta files are available for download here: https://usegalaxy.org/u/carlosfarkas/h/sars-cov-2-variants-gisaid-august-03-2020 and can be downloaded with wget: 

```
wget -O gisaid_Africa_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5df5a9de556b60745/display?to_ext=fasta.gz
wget -O gisaid_Asia_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5c7bff6a669e318dc/display?to_ext=fasta.gz
wget -O gisaid_Europe_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b507b027e055bf2df9/display?to_ext=fasta.gz
wget -O gisaid_North_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ee919645a4a97d76/display?to_ext=fasta.gz
wget -O gisaid_Oceania_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ad62fc70fed0a55b/display?to_ext=fasta.gz
wget -O gisaid_South_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5134c7103a63c1db1/display?to_ext=fasta.gz
wget -O merged.GISAID.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b50b3becb49899ed42/display?to_ext=fasta.gz
```

### Execution

As an example for merged.GISAID.fasta.gz (containing worldwide GISAID genomes) we can obtain aggregated variants from merged.GISAID.fasta.gz dataset in a folder called "GISAID_merge" in Ubuntu as follows:

```
# Clone repository anywhere
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git
samtools faidx ./SARS-CoV-2-freebayes/covid19-refseq.fasta && chmod 755 ./SARS-CoV-2-freebayes/SARS-CoV-2* ./SARS-CoV-2-freebayes/covid19-refseq.fasta*

# In the previous directory, download merged.GISAID.fasta.gz inside "GISAID_merge" folder and decompress.
mkdir GISAID_merge && cd GISAID_merge
wget -O merged.GISAID.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b50b3becb49899ed42/display?to_ext=fasta.gz && gunzip merged.GISAID.fasta.gz

# Execute the pipeline using 10 threads
ulimit -n 1000000 && ulimit -s 299999  # increase stack size and open file limit, see README_ulimit for details.
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh merged.GISAID.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
```

-This operation will obtain aggregated variants per region (merged.GISAID.AF.vcf) and aggregated variants filtered with Viral Frequency > 1% (merged.GISAID.AF_1%.vcf) inside the folder "GISAID_merge". Users can change the name of the folder (i.e.: GISAID_North_America for North America GISAID genomes). Execution takes 98000 cpu seconds (~27.2 hrs) with a peak of ~ 4GB of RAM.   

-NOTE: It is recommended to process larger FASTA collections by chunks (i.e. chunks of 100000 genomes). We provide up to date analysis of GISAID genomes here: https://github.com/cfarkas/SARS-CoV-2-freebayes/wiki (November 2020) using ~230000 genomes. Ubuntu users need to change ulimit -s and -n parameters, see README_ulimit file in this repository for details. 

### Number of variants per genome

Users can compute the number of variants per genome as follows. Inside all_variants folder, do:

```
ulimit -n 1000000 && ulimit -s 299999
{
vcf= ls -1 EPI_*.vcf
for vcf in EPI_*.vcf; do grep -P 'NC_045512.2\t' ${vcf} -c
done
#
} | tee logfile_variants_GISAID_freebayes
#
grep "EPI_ISL_" logfile_variants_GISAID_freebayes > vcf_files
grep -v "EPI_ISL_" logfile_variants_GISAID_freebayes > variants_per_sample
paste vcf_files variants_per_sample > logfile_variants_GISAID
rm vcf_files variants_per_sample
sed -i 's/.fa.left.vcf//'g logfile_variants_GISAID
rm logfile_variants_GISAID_freebayes
```
logfile_variants_GISAID file contains the GISAID accession along with the number of detected variants. 


## III) Collecting variants per protein on SnpEff-classified GISAID variants
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


## IV) Nucleotide diversity and Tajima's D calculation per geographical region
To estimate nucleotide diversity (π) and Tajima's D parameters we will employ vcftools program version from Julien Y. Dutheil (accepting --haploid flag) (https://github.com/jydu/vcftools). We will download with wget FASTA genomes from each continent submitted to GISAID until August 03, 2020 and we will execute from scratch variant calling analysis including the vcftools analysis. In a folder (i.e. diversity). From scratch, do: 

```
mkdir diversity
cd ./diversity/
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.gff.gz
gunzip *
gffread GCF_009858895.2_ASM985889v3_genomic.gff -T -o SARS-CoV-2.gtf
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git
samtools faidx ./SARS-CoV-2-freebayes/covid19-refseq.fasta && chmod 755 ./SARS-CoV-2-freebayes/SARS-CoV-2* ./SARS-CoV-2-freebayes/covid19-refseq.fasta*

### Full analysis, using 10 threads 

### Africa                                 
mkdir GISAID_Africa && cd GISAID_Africa
wget -O gisaid_Africa_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5df5a9de556b60745/display?to_ext=fasta.gz && gunzip gisaid_Africa_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 299999
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Africa_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 100 --haploid --out Africa
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 100 --haploid --out Africa

### Asia                                 
mkdir GISAID_Asia && cd GISAID_Asia
wget -O gisaid_Asia_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5c7bff6a669e318dc/display?to_ext=fasta.gz && gunzip gisaid_Asia_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 299999
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Asia_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 100 --haploid --out Asia
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 100 --haploid --out Asia

### Europe                                 
mkdir GISAID_Europe && cd GISAID_Europe
wget -O gisaid_Europe_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b507b027e055bf2df9/display?to_ext=fasta.gz && gunzip gisaid_Europe_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 299999
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Europe_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 100 --haploid --out Europe
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 100 --haploid --out Europe

### North_America                                
mkdir GISAID_North_America && cd GISAID_North_America  
wget -O gisaid_North_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ee919645a4a97d76/display?to_ext=fasta.gz && gunzip gisaid_North_America_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 299999
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_North_America_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 100 --haploid --out North_America
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 100 --haploid --out North_America

### South_America                                
mkdir GISAID_South_America && cd GISAID_South_America
wget -O gisaid_South_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5134c7103a63c1db1/display?to_ext=fasta.gz && gunzip gisaid_South_America_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 299999
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_South_America_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 100 --haploid --out South_America
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 100 --haploid --out South_America

### Oceania                               
mkdir GISAID_Oceania && cd GISAID_Oceania
wget -O gisaid_Oceania_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ad62fc70fed0a55b/display?to_ext=fasta.gz && gunzip gisaid_Oceania_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 299999
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Oceania_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 100 --haploid --out Oceania
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 100 --haploid --out Oceania
```


## V) inStrain analysis of SRA sequencing cohorts (microdiversity)
To estimate nucleotide diversity (microdiversity within a sequencing sample), analysis of SNV linkage and coverage analysis, we will employ the inStrain package, written in python (https://instrain.readthedocs.io/en/latest/). To analyze Sequence Read Archive datasets from USA as done in the manuscript, do the following: 

```
mkdir inStrain_USA
cd inStrain_USA

# Download and convert to fastq.gz all USA accessions from this repository (inStrain_USA.tabular), using 40 threads

git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git
cp ./SARS-CoV-2-freebayes/SARS-CoV-2.gb ./
cp ./SARS-CoV-2-freebayes/covid19-refseq.fasta* ./
samtools faidx covid19-refseq.fasta
bash SARS-CoV-2-freebayes/SARS-CoV-2-freebayes.sh ./SARS-CoV-2-freebayes/inStrain_USA.tabular covid19-refseq.fasta 40

#############################
### inStrain_USA analysis ###
#############################   

# Execute inStrain, using 40 threads (take a while)

bam= ls -1 *.bam
for bam in *.bam; do inStrain profile ${bam} covid19-refseq.fasta --gene_file SARS-CoV-2.gb --processes 40 -o ./${bam}.IS; done

# Copy outputs
IS= ls -1 *.IS/output/*.bam.IS_genome_info.tsv
for IS in *.IS/output/*.bam.IS_genome_info.tsv; do cp ${IS} ./ ; done

# Fill empty spaces with asterisks
IS= ls -1 *.bam.IS_genome_info.tsv
for IS in *.bam.IS_genome_info.tsv; do sed -i 's/\t\t/\t*\t/' ${IS}; done

# Transpose outputs
IS= ls -1 *.bam.IS_genome_info.tsv
for IS in *.bam.IS_genome_info.tsv; do python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < ${IS} > ${IS}.transpose
done

# convert delimiters to tab
IS= ls -1 *.bam.IS_genome_info.tsv.transpose
for IS in *.bam.IS_genome_info.tsv.transpose; do sed -i 's/ /\t/'g ${IS}; done
rm *.bam.IS_genome_info.tsv

# Renaming transposed files
for filename in *.bam.IS_genome_info.tsv.transpose; do mv "./$filename" "./$(echo "$filename" | sed -e 's/.bam.IS_genome_info.tsv.transpose/.transpose/g')";  done

# add FILENAME inside files and remove extension
transpose= ls -1 *transpose
for transpose in *transpose; do sed -i "s/all_scaffolds/${transpose}/"g ${transpose}; done
transpose= ls -1 *transpose
for transpose in *transpose; do sed -i 's/.transpose//'g ${transpose}; done

# Join files with a function
bash ./SARS-CoV-2-freebayes/join_inStrain.sh *.transpose > joined
sed -i 's/ /\t/'g joined

# transpose merged.tab and generate final file (merged.tabular)
python -c "import sys; print('\n'.join(' '.join(c) for c in zip(*(l.split() for l in sys.stdin.readlines() if l.strip()))))" < joined > merged.tabular
rm joined
sed -i 's/ /\t/'g merged.tabular
sed -i 's/NULL/*/'g merged.tabular

### Freebayes variants (VF>5%, related to figure 1G). 
a= ls -1 *.bam 
for a in *.bam; do freebayes-parallel <(fasta_generate_regions.py covid19-refseq.fasta.fai 2000) 50 -f covid19-refseq.fasta -F 0.0499 -b ${a} > ${a}.freebayes.vcf
done 

{
vcf= ls -1 *.freebayes.vcf
for vcf in *.freebayes.vcf; do grep -P 'NC_045512.2\t' ${vcf} -c
done
#
} | tee logfile_variants_AF_5%_freebayes
#
grep ".bam.freebayes.vcf" logfile_variants_AF_5%_freebayes > vcf_files
grep -v ".bam.freebayes.vcf" logfile_variants_AF_5%_freebayes > variants_per_sample
paste vcf_files variants_per_sample > logfile_variants_AF_5%_freebayes.tabular
rm vcf_files variants_per_sample
sed -i s'/.bam.freebayes.vcf//'g logfile_variants_AF_5%_freebayes.tabular

```
merged.tabular contain all parameters calculated by inStrain for each genome, aggregated in a single file. logfile_variants_AF_5%_freebayes.tabular contain variants with viral frequency over 5%, per genome.  
