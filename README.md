# SARS-CoV-2-freebayes
Analysis of SARS-CoV-2 genome variants collected with freebayes variant caller.

## Installation:  

### Option 1: Via conda/pip (recommended)
- requires miniconda, python2.7 and python>=3. To install miniconda, see: https://docs.conda.io/en/latest/miniconda.html
```
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git  # clone repository
cd SARS-CoV-2-freebayes
conda config --add channels conda-forge                        # add conda-forge channel (if you haven't already done so)
conda config --add channels bioconda                           # add bioconda channel (if you haven't already done so)
conda env update --file environment.yml                        # install required programs 
conda activate SARS-CoV-2-freebayes                            # activate SARS-CoV-2-freebayes enviroment 
git clone https://github.com/cfarkas/vcftools.git              # install forked vcftools repository
cd vcftools
./autogen.sh
./configure
make
sudo make install
cd ..
```
After these steps, a conda enviroment called SARS-CoV-2-freebayes can be managed as follows:
```
# To activate this environment, use
#
#     $ conda activate SARS-CoV-2-freebayes
#
# To deactivate an active environment, use
#
#     $ conda deactivate
```
By activating the enviroment, all scripts in the SARS-CoV-2-freebayes repository can be executed, without further installations.

### Option 2: Without using conda, program by program:

- see detailed installation steps in our wiki here: https://github.com/cfarkas/SARS-CoV-2-freebayes/wiki#without-using-conda-program-by-program

### Note of performance and ulimit before using the scripts
For users working with small number of genomes (i.e.: < 60000) or have unlimited values of -s and -n parameters in their machines, use "nolimit" versions of scripts:
```
SARS-CoV-2-NGS-freebayes-nolimit.sh
SARS-CoV-2-GISAID-freebayes-nolimit.sh
SARS-CoV-2-processing-fasta-nolimit.sh
SARS-CoV-2-merge-variants-nolimit.sh
```
For a proper performance working with high number of genomes (i.e. > 60000), use these scripts: 
```
SARS-CoV-2-NGS-freebayes.sh
SARS-CoV-2-GISAID-freebayes.sh
SARS-CoV-2-processing-fasta.sh
SARS-CoV-2-merge-variants.sh
```
The latter scripts asume you can set these values: 
```
ulimit -n 1000000    # Check your machine with: ulimit -n 
ulimit -s 1000000      # Check your machine with: ulimit -s
```
Check README_ulimit for how to change these values in Ubuntu.

## Execution:

## I) Colecting Variants (Sequence Read Archive datasets)

In order to obtain SARS-CoV-2 variants (viral frequency >= 0.5) users need to provide:  

- SRA_list: Sequence read archive accessions of each datasets (SRR prefix list, in tabular format or txt format). As example, see SRA_Accessions_Aug_03_2020.tabular, containing curated SARS-CoV-2 accessions from SRA until August 03, 2020, provided in this repository.
- SARS-CoV-2 reference in fasta format (covid19-refseq.fasta, provided in this repository)
- VF: viral frequency threshold for variant calling (a number between 0-1). As example 0.5, would mean 50% viral frequency variants as threshold. 
- Threads: number of threads for calculations 

Execution (from scratch): 
```
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git
cd SARS-CoV-2-freebayes
samtools faidx covid19-refseq.fasta
chmod 755 SARS-CoV-2* covid19-refseq.fasta* *.sh
./SARS-CoV-2-NGS-freebayes.sh SRA_list covid19-refseq.fasta VF Threads
```
This execution will:

- Download SRA datasets from the provided list, convert to fastq, trim adaptors and gzip reads, for each line of the provided list
- Align trimmed reads against SARS-CoV-2 reference genome (NC_045512.2) by using minimap2
- Call variants at a certain viral frequency by using freebayes as frequency-based pooled caller
- Merge all variants in a single VCF file by using jacquard. This VCF file also contains viral frequencies in the AF field (AF=AO/AO+RO). See AO and RO fields for alternative and reference allele counts.

### Example 

We provided a curated list of 17560 SARS-CoV-2 worldwide datasets until July 28, 2020 (see SRA_Accessions_Jul_28_2020.tabular). We also provided curated lists in txt format by continent (see July_28_2020_*.txt files). As an example, we will collect variants from July_28_2020_North_America.txt datasets using 30 threads. In Ubuntu it is recommended to increase open file limit and stack size accordingly at the number of genomes to process, otherwise these steps may crush. see README_ulimit for details in Ubuntu. 

```
ulimit -n 1000000 && ulimit -s 1000000    # check if you can increase stack size and open file limit, see README_ulimit for details.
./SARS-CoV-2-NGS-freebayes.sh July_28_2020_North_America.txt covid19-refseq.fasta 0.4999 30
```
will collect variants (VF>=0.5) in each Sample.

### Number of variants per genome

Users can compute the number of variants per genome as follows. Inside SARS-CoV-2 folder containing VCF files, do:

```
ulimit -n 1000000 && ulimit -s 1000000    # check if you can set these values
{
vcf= ls -1 *.bam.freebayes.vcf
for vcf in *.bam.freebayes.vcf; do grep -P 'NC_045512.2\t' ${vcf} -c
done
#
} | tee logfile_variants_NGS_freebayes
#
grep ".vcf" logfile_variants_NGS_freebayes > vcf_files
grep -v ".vcf" logfile_variants_NGS_freebayes > variants_per_sample
paste vcf_files variants_per_sample > logfile_variants_NGS
rm vcf_files variants_per_sample logfile_variants_NGS_freebayes
```

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

As an example for merged.GISAID.fasta.gz (containing worldwide GISAID genomes until August 03, 2020) we can obtain aggregated variants from merged.GISAID.fasta.gz FASTA dataset in a folder called "GISAID_merge" in Ubuntu. From scratch:

```
# Clone repository
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git
samtools faidx ./SARS-CoV-2-freebayes/covid19-refseq.fasta && chmod 755 ./SARS-CoV-2-freebayes/*.sh ./SARS-CoV-2-freebayes/covid19-refseq.fasta*

# Create GISAID_merge folder, enter it, download merged.GISAID.fasta.gz and decompress.
mkdir GISAID_merge && cd GISAID_merge
wget -O merged.GISAID.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b50b3becb49899ed42/display?to_ext=fasta.gz && gunzip merged.GISAID.fasta.gz

# Execute the pipeline from SARS-CoV-2-freebayes folder. Using 10 threads:
ulimit -n 1000000 && ulimit -s 1000000  # check if you can increase stack size and open file limit, see README_ulimit for details.
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh merged.GISAID.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
```

-This operation will obtain joint calls in a single vcf containing all samples (merged.GISAID.AF.vcf). In this matrix, zeros indicate absence of variant and ones indicate the presence of the variant, per sample. Viral frequencies (AF field) were also added. An intermediate file will be also generated: combined_sites.vcf containing just merged variants
- Users can change the name of the folder (i.e.: GISAID_North_America for North America GISAID genomes). Execution takes 98000 cpu seconds (~27.2 hrs) with a peak of ~ 4GB of RAM.   

-NOTE: It is recommended to process larger FASTA collections by chunks (i.e. chunks of 100000 genomes). We provide up to date analysis of GISAID genomes here: https://github.com/cfarkas/SARS-CoV-2-freebayes/wiki (November 2020) using ~230000 genomes. Ubuntu users need to change ulimit -s and -n parameters, see README_ulimit file in this repository for details.

### Number of variants per genome

Users can compute the number of variants per genome as follows. Inside all_variants folder, do:

```
ulimit -n 1000000 && ulimit -s 1000000  # check if you can set these values
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


## III) Annotate and collect variants per protein on SnpEff-classified VCF variants
merged.GISAID.AF.vcf files can weight several gigabytes and therefore we compress it with gzip for storage. Compressed merged VCF files from GISAID genomes and SRA, including effect annotation are available here: https://usegalaxy.org/u/carlosfarkas/h/snpeffsars-cov-2. To speed-up things, prior to SnpEff annotation, it is preferable to lightweight these files by selecting one individual sample in the vcf (dropping all the others samples) or annotate combined_sites.vcf file. We will execute the first option.

```
mkdir SnpEff-SARS-CoV-2
cd SnpEff-SARS-CoV-2
wget -O merged.GISAID.Aug-03-2020.vcf.gz https://usegalaxy.org/datasets/bbd44e69cb8906b580d64fc18488c56a/display?to_ext=gff3.gz
wget -O merged.SRA.Jul-28-2020.vcf.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5f816689b9933202f/display?to_ext=gff3.gz
wget -O merged.GISAID.Nov-30-2020.vcf.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5509c81bce034b0a2/display?to_ext=gff3.gz
gunzip merged.GISAID.Aug-03-2020.vcf.gz # ~ 3.4 GB
gunzip merged.GISAID.Nov-30-2020.vcf.gz # ~ 17 GB
gunzip merged.SRA.Jul-28-2020.vcf.gz    # ~ 547 Mb
vcfkeepsamples merged.GISAID.Aug-03-2020.vcf EPI_ISL_402119 > merged.GISAID.Aug-03-2020.EPI_ISL_402119.vcf
gzip merged.GISAID.Aug-03-2020.vcf
vcfkeepsamples merged.SRA.Jul-28-2020.vcf ERR4082713 > merged.SRA.Jul-28-2020.ERR4082713.vcf
gzip merged.SRA.Jul-28-2020.vcf
vcfkeepsamples merged.GISAID.Nov-30-2020.vcf EPI_ISL_402119 > merged.GISAID.Nov-30-2020.EPI_ISL_402119.vcf
gzip merged.GISAID.Nov-30-2020.vcf
```
Then, merged.GISAID.Aug-03-2020.EPI_ISL_402119.vcf, merged.SRA.Jul-28-2020.ERR4082713.vcf and merged.GISAID.Nov-30-2020.EPI_ISL_402119.vcf files can be uploaded here: https://usegalaxy.org/ and annotated using the tool SnpEff eff: annotate variants for SARS-CoV-2 (default mode), outputting an annotated VCF file including an associated HTML file. As example, SnpEff-annotated VCF outputs are available here: https://usegalaxy.org/u/carlosfarkas/h/snpeffsars-cov-2 and can be processed as follows:
```
# SRA variants
mkdir SnpEff-Jul-28-2020.SRA && cd SnpEff-Jul-28-2020.SRA                                                               # 1) Create folder and enter it
wget -O SnpEff-Jul-28-2020.SRA.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b5b7e3ff0964c68fa2/display?to_ext=vcf   # 2) Download annotated vcf file
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git && chmod 755 ./SARS-CoV-2-freebayes/*.sh                  # 3) Download repo and change permissions
./SARS-CoV-2-freebayes/SnpEff_processing.sh SnpEff-Jul-28-2020.SRA.vcf                                                  # 4) execute SnpEff_processing.sh 
cd ..

# GISAID variants: August 03, 2020                        
mkdir SnpEff-Aug-03-2020.GISAID && cd SnpEff-Aug-03-2020.GISAID                                                          # 1) Create folder and enter it
wget -O SnpEff-Aug-03-2020.GISAID.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b5117bba070d8c6bca/display?to_ext=vcf # 2) Download annotated vcf file
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git && chmod 755 ./SARS-CoV-2-freebayes/*.sh                   # 3) Download repo and change permissions
./SARS-CoV-2-freebayes/SnpEff_processing.sh SnpEff-Aug-03-2020.GISAID.vcf                                                # 4) execute SnpEff_processing.sh 
cd ..

# GISAID variants: November 30, 2020
mkdir SnpEff-Nov-30-2020.GISAID && cd SnpEff-Nov-30-2020.GISAID                                                          # 1) Create folder and enter it
wget -O SnpEff-Nov-30-2020.GISAID.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b5743f34d0337d9459/display?to_ext=vcf # 2) Download annotated vcf file
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git && chmod 755 ./SARS-CoV-2-freebayes/*.sh                   # 3) Download repo and change permissions
./SARS-CoV-2-freebayes/SnpEff_processing.sh SnpEff-Nov-30-2020.GISAID.vcf                                                # 4) execute SnpEff_processing.sh 
```
In each folder, variants_per_protein subfolder contain variants per protein. Files ending in ".SnpEff" contains parsed variants per consequence and ".counts" contains associated counts, respectively. Also, the script computed aminoacid changes (see SnpEff.AAchanges files). We suggest user-provided vcf files should be processed in a likewise manner, placing the vcf file in a specific folder and executing steps 3) and 4).

## IV) Nucleotide diversity and Tajima's D test calculation per geographical region
To estimate nucleotide diversity (π) and Tajima's D test, we will employ vcftools program version from Julien Y. Dutheil (accepting --haploid flag) (https://github.com/jydu/vcftools). We will download with wget FASTA genomes from each continent submitted to GISAID until August 03, 2020 and we will execute variant calling and vcftools analysis, using a sliding window of 50 bp (can be changed). An excellent explanation of Tajima's D test can be found here: https://www.youtube.com/watch?v=wiyay4YMq2A .

From scratch, the whole analysis can be done in a folder (i.e. diversity), as presented below. If users already executed SARS-CoV-2-GISAID-freebayes.sh, can skip this step and proceed directly to execute the vcftools commands.  

```
mkdir diversity
cd diversity
git clone https://github.com/cfarkas/SARS-CoV-2-freebayes.git
samtools faidx ./SARS-CoV-2-freebayes/covid19-refseq.fasta && chmod 755 ./SARS-CoV-2-freebayes/*.sh ./SARS-CoV-2-freebayes/covid19-refseq.fasta*

### Full analysis, using 10 threads 

### Africa                                 
mkdir GISAID_Africa && cd GISAID_Africa
wget -O gisaid_Africa_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5df5a9de556b60745/display?to_ext=fasta.gz && gunzip gisaid_Africa_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                           # check if you can set these values
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Africa_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Africa.50    # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Africa.50      # 50 bp sliding window
cd ..

### Asia                                 
mkdir GISAID_Asia && cd GISAID_Asia
wget -O gisaid_Asia_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5c7bff6a669e318dc/display?to_ext=fasta.gz && gunzip gisaid_Asia_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                         # check if you can set these values
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Asia_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Asia.50    # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Asia.50      # 50 bp sliding window
cd ..

### Europe                                 
mkdir GISAID_Europe && cd GISAID_Europe
wget -O gisaid_Europe_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b507b027e055bf2df9/display?to_ext=fasta.gz && gunzip gisaid_Europe_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                           # check if you can set these values
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Europe_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Europe.50    # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Europe.50      # 50 bp sliding window
cd ..

### North_America                                
mkdir GISAID_North_America && cd GISAID_North_America  
wget -O gisaid_North_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ee919645a4a97d76/display?to_ext=fasta.gz && gunzip gisaid_North_America_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                                  # check if you can set these values
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_North_America_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out North_America.50    # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out North_America.50      # 50 bp sliding window
cd ..

### South_America                                
mkdir GISAID_South_America && cd GISAID_South_America
wget -O gisaid_South_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5134c7103a63c1db1/display?to_ext=fasta.gz && gunzip gisaid_South_America_08_03_2020.fasta.gz 
ulimit -n 1000000 && ulimit -s 1000000                                                  # check if you can set these values
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_South_America_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out South_America.50    # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out South_America.50      # 50 bp sliding window
cd ..

### Oceania                               
mkdir GISAID_Oceania && cd GISAID_Oceania
wget -O gisaid_Oceania_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ad62fc70fed0a55b/display?to_ext=fasta.gz && gunzip gisaid_Oceania_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                            # check if you can set these values
../SARS-CoV-2-freebayes/SARS-CoV-2-GISAID-freebayes.sh gisaid_Oceania_08_03_2020.fasta ../SARS-CoV-2-freebayes/covid19-refseq.fasta 10
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Oceania.50    # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Oceania.50      # 50 bp sliding window
```

## V) π versus Tajima's D values per bin
Inside every geographical folder, 50 bp bins containing π and Tajima's D values can be joined by bin for plotting purposes, as depicted here: https://doi.org/10.1016/j.tig.2006.06.005.  As example for π and Tajima's D values from every geographical region:
```  
# Africa
awk '{print $1"\t"$3"\t"$4"\t"$5}' Africa.50.windowed.pi > Africa.50.subset.windowed.pi
join -j 2 <(sort -k2 Africa.50.subset.windowed.pi) <(sort -k2 Africa.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Africa.50.pi.D && rm joined Africa.50.subset.windowed.pi

# Asia
awk '{print $1"\t"$3"\t"$4"\t"$5}' Asia.50.windowed.pi > Asia.50.subset.windowed.pi
join -j 2 <(sort -k2 Asia.50.subset.windowed.pi) <(sort -k2 Asia.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Asia.50.pi.D && rm joined Asia.50.subset.windowed.pi

# Europe
awk '{print $1"\t"$3"\t"$4"\t"$5}' Europe.50.windowed.pi > Europe.50.subset.windowed.pi
join -j 2 <(sort -k2 Europe.50.subset.windowed.pi) <(sort -k2 Europe.50000.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Europe.50.pi.D && rm joined Europe.50.subset.windowed.pi

# North_America
awk '{print $1"\t"$3"\t"$4"\t"$5}' North_America.50.windowed.pi > North_America.50.subset.windowed.pi
join -j 2 <(sort -k2 North_America.50.subset.windowed.pi) <(sort -k2 North_America.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > North_America.50.pi.D && rm joined North_America.50.subset.windowed.pi

# South_America
awk '{print $1"\t"$3"\t"$4"\t"$5}' South_America.50.windowed.pi > South_America.50.subset.windowed.pi
join -j 2 <(sort -k2 South_America.50.subset.windowed.pi) <(sort -k2 South_America.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > South_America.50.pi.D && rm joined South_America.50.subset.windowed.pi

# Oceania
awk '{print $1"\t"$3"\t"$4"\t"$5}' Oceania.50.windowed.pi > Oceania.50.subset.windowed.pi
join -j 2 <(sort -k2 Oceania.50.subset.windowed.pi) <(sort -k2 Oceania.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Oceania.50.pi.D && rm joined Oceania.50.subset.windowed.pi
```  

Files ended in "50.pi.D" contains π and Tajima's D values by bin and can be plotted in any sofware. We provided in this repository a script called pi-tajima.sh to process these files to obtain regions and variants falling outside 95% confidence interval of Tajima's D, useful for further study. This script requires vcfstats and some R libraries installed, as depicted here: https://github.com/cfarkas/SARS-CoV-2-freebayes#10-install-ggplot2-ggrepel-and-vcfr-r-libraries. 

For example to process Africa.50.pi.D and merged.GISAID.AF.vcf (from Africa), do the following:

```
/full/path/to/SARS-CoV-2-freebayes/pi-tajima.sh Africa.50.pi.D merged.GISAID.AF.vcf  
```
where /full/path/to/ is the full path to this repository. merged.GISAID.AF.vcf can be replaced by the annotated SnpEff version of this file as well.

after executing this script, check "postprocessing_pi_D_output_files"  folder, containing:
- bins outside 2.5 and 97.5% confidence intervals are outliers for further study (check bins_2.5%_confidencce.bed and bins_2.5%_confidencce.bed files). 
- pi_tajima.pdf plot, for a visual inspection of these regions.
- 2.5_CI_confidence.recode.vcf.gz and 97.5_CI_confidence.recode.vcf.gz, gzipped vcf files containing variants falling in these regions
- vcfstats plots and associated data.
