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
bash makefile.sh                                               # make & install
sudo cp ./bin/* /usr/local/bin/                                # Copy binaries into /usr/local/bin/ (require sudo privileges)
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

- Finally, we need to install a vcftools version that includes the --haploid flag as follows:
```
git clone https://github.com/cfarkas/vcftools.git       # Forked Julien Y. Dutheil version: https://github.com/jydu/vcftools
cd vcftools
./autogen.sh
./configure
make                                                    # make
sudo make install                                       # requires sudo privileges
```

### Option 2: Without using conda, program by program:

- see detailed installation steps in our wiki here: https://github.com/cfarkas/SARS-CoV-2-freebayes/wiki#without-using-conda-program-by-program

### Note of performance and ulimit before using the binaries
For users working with small number of genomes (i.e.: < 60000) or have unlimited values of -s and -n parameters in their machines, use "nolimit" versions of binaries:
```
SARS-CoV-2-NGS-freebayes-nolimit
SARS-CoV-2-GISAID-freebayes-nolimit
SARS-CoV-2-processing-fasta-nolimit
SARS-CoV-2-merge-variants-nolimit
```
For a proper performance working with high number of genomes (i.e. > 60000), use these binaries: 
```
SARS-CoV-2-NGS-freebayes
SARS-CoV-2-GISAID-freebayes
SARS-CoV-2-processing-fasta
SARS-CoV-2-merge-variants
```
The latter binaries asume you can set these values: 
```
ulimit -n 1000000    # Check your machine with: ulimit -n 
ulimit -s 1000000    # Check your machine with: ulimit -s
```
Check README_ulimit for how to change these values in Ubuntu.

## Execution:

## I) Colecting Variants (Sequence Read Archive datasets)

In order to obtain SARS-CoV-2 variants (viral frequency >= 0.5) users need to provide:  

- l: SRA_list: Sequence read archive accessions of each datasets (SRR prefix list, in tabular format or txt format). As example, see SRA_Accessions_Aug_03_2020.tabular, containing curated SARS-CoV-2 accessions from SRA until August 03, 2020, provided in this repository.
- g: SARS-CoV-2 reference in fasta format (covid19-refseq.fasta, provided in this repository)
- a: viral frequency threshold for variant calling (a number between 0-1). As example 0.5, would mean 50% viral frequency variants as threshold. 
- t: number of threads for calculations 

Execution, assuming that binaries are in ```/usr/local/bin```: 
```
samtools faidx covid19-refseq.fasta
SARS-CoV-2-NGS-freebayes -l <SRA_list> -g <covid19-refseq.fasta> -a <VF> -t <Threads>
```
This execution will:

- Download SRA datasets from the provided list, convert to fastq, trim adaptors and gzip reads, for each line of the provided list
- Align trimmed reads against SARS-CoV-2 reference genome (NC_045512.2) by using minimap2
- Call variants at a certain viral frequency by using freebayes as frequency-based pooled caller
- Merge all variants in a single VCF file by using jacquard. This VCF file also contains viral frequencies in the AF field (AF=AO/AO+RO). See AO and RO fields for alternative and reference allele counts.

### Example 

We provided a curated list of 17560 SARS-CoV-2 worldwide datasets until July 28, 2020 (see SRA_Accessions_Jul_28_2020.tabular). We also provided curated lists in txt format by continent (see July_28_2020_*.txt files). As an example, we will collect variants from July_28_2020_North_America.txt datasets using 30 threads. In Ubuntu it is recommended to increase open file limit and stack size accordingly at the number of genomes to process, otherwise these steps may crush. see README_ulimit for details in Ubuntu. In a given folder: 

```
ulimit -n 1000000 && ulimit -s 1000000    # check if you can increase stack size and open file limit, see README_ulimit for details.
# Execute the pipeline, providing full path to July_28_2020_North_America.txt and covid19-refseq.fasta files:
SARS-CoV-2-NGS-freebayes -l July_28_2020_North_America.txt -g covid19-refseq.fasta -a 0.4999 -t 30
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

As an example for merged.GISAID.fasta.gz (containing worldwide GISAID genomes until August 03, 2020) we can obtain aggregated variants from merged.GISAID.fasta.gz FASTA dataset. Assuming binaries are in ```/usr/local/bin``` and users previously did ```samtools faidx covid19-refseq.fasta```:

```
# In a given folder, download merged.GISAID.fasta.gz and decompress.
wget -O merged.GISAID.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b50b3becb49899ed42/display?to_ext=fasta.gz && gunzip merged.GISAID.fasta.gz

# Execute the pipeline, providing full path to merged.GISAID.fasta and covid19-refseq.fasta sequences:
ulimit -n 1000000 && ulimit -s 1000000  # check if you can increase stack size and open file limit, see README_ulimit for details.
SARS-CoV-2-GISAID-freebayes -f merged.GISAID.fasta -g covid19-refseq.fasta -t 10 
```
Execution takes 98000 cpu seconds (~27.2 hrs) in a regular ubuntu workstation with a peak of ~ 4GB of RAM.

-This operation will obtain joint calls in a single vcf containing all samples (merged.GISAID.AF.vcf). In this matrix, zeros indicate absence of variant and ones indicate the presence of the variant, per sample. Viral frequencies (AF field) were also added. An intermediate file will be also generated: combined_sites.vcf containing just merged variants

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
- Then, merged.GISAID.Aug-03-2020.EPI_ISL_402119.vcf, merged.SRA.Jul-28-2020.ERR4082713.vcf and merged.GISAID.Nov-30-2020.EPI_ISL_402119.vcf files can be uploaded here: https://usegalaxy.org/ and annotated using the tool SnpEff eff: annotate variants for SARS-CoV-2 (default mode), outputting an annotated VCF file including an associated HTML file. As example, SnpEff-annotated VCF outputs are available here: https://usegalaxy.org/u/carlosfarkas/h/snpeffsars-cov-2 and can be processed as follows. 
- Assuming binaries are in ```/usr/local/bin``` and users previously did ```samtools faidx covid19-refseq.fasta```:
```
# SRA variants
mkdir SnpEff-Jul-28-2020.SRA && cd SnpEff-Jul-28-2020.SRA                                                               # 1) Create folder and enter it
wget -O SnpEff-Jul-28-2020.SRA.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b5b7e3ff0964c68fa2/display?to_ext=vcf   # 2) Download annotated vcf file
SnpEff_processing -v SnpEff-Jul-28-2020.SRA.vcf                                                                         # 3) execute SnpEff_processing 
cd ..

# GISAID variants: August 03, 2020                        
mkdir SnpEff-Aug-03-2020.GISAID && cd SnpEff-Aug-03-2020.GISAID                                                          # 1) Create folder and enter it
wget -O SnpEff-Aug-03-2020.GISAID.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b5117bba070d8c6bca/display?to_ext=vcf # 2) Download annotated vcf file
SnpEff_processing -v SnpEff-Aug-03-2020.GISAID.vcf                                                                       # 3) execute SnpEff_processing 
cd ..

# GISAID variants: November 30, 2020
mkdir SnpEff-Nov-30-2020.GISAID && cd SnpEff-Nov-30-2020.GISAID                                                          # 1) Create folder and enter it
wget -O SnpEff-Nov-30-2020.GISAID.vcf https://usegalaxy.org/datasets/bbd44e69cb8906b5743f34d0337d9459/display?to_ext=vcf # 2) Download annotated vcf file
SnpEff_processing -v SnpEff-Nov-30-2020.GISAID.vcf                                                                       # 3) execute SnpEff_processing 
```
In each folder, variants_per_protein subfolder contain variants per protein. Files ending in ".SnpEff" contains parsed variants per consequence and ".counts" contains associated counts, respectively. Also, the script computed aminoacid changes (see SnpEff.AAchanges files). We suggest user-provided vcf files should be processed in a likewise manner, placing the vcf file in a specific folder and executing steps 3) and 4).

## IV) Nucleotide diversity and Tajima's D test calculation per geographical region
- To estimate nucleotide diversity (π) and Tajima's D test, we will employ vcftools program version from Julien Y. Dutheil (accepting --haploid flag) (https://github.com/jydu/vcftools). 
- We will download with wget FASTA genomes from each continent submitted to GISAID until August 03, 2020 and we will execute variant calling and vcftools analysis, using a sliding window of 50 bp across SARS-CoV-2 genome (can be changed). Then, we will we will use the pi-tajima.sh script to create a plot of pi and Tajima's D values per bin, as depicted here: https://doi.org/10.1016/j.tig.2006.06.005
- An excellent explanation of Tajima's D test can be found here: https://www.youtube.com/watch?v=wiyay4YMq2A. 
- pi-tajima.sh script requires vcfstats and some R libraries installed, as depicted here: https://github.com/cfarkas/SARS-CoV-2-freebayes#10-install-ggplot2-ggrepel-and-vcfr-r-libraries. 

From scratch, the whole analysis can be done in a folder (i.e. diversity), as presented below. Assuming binaries are in ```/usr/local/bin``` and users previously did ```samtools faidx covid19-refseq.fasta```: 

```
mkdir diversity
cd diversity

### Full analysis, using 10 threads 

### Africa                                 
mkdir GISAID_Africa && cd GISAID_Africa
wget -O gisaid_Africa_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5df5a9de556b60745/display?to_ext=fasta.gz && gunzip gisaid_Africa_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                                       # check if you can set these values
SARS-CoV-2-GISAID-freebayes -f gisaid_Africa_08_03_2020.fasta -g covid19-refseq.fasta -t 10  # call variants, provide full path to files
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Africa.50                 # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Africa.50                   # 50 bp sliding window

# joining pi with Tajima's D values and calculate Tajima's D confidence intervals
awk '{print $1"\t"$3"\t"$4"\t"$5}' Africa.50.windowed.pi > Africa.50.subset.windowed.pi
join -j 2 <(sort -k2 Africa.50.subset.windowed.pi) <(sort -k2 Africa.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Africa.50.pi.D && rm joined Africa.50.subset.windowed.pi
pi-tajima -f Africa.50.pi.D -v merged.GISAID.AF.vcf                                          # Run pi-tajima               
cd ..

### Asia                                 
mkdir GISAID_Asia && cd GISAID_Asia
wget -O gisaid_Asia_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5c7bff6a669e318dc/display?to_ext=fasta.gz && gunzip gisaid_Asia_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                                       # check if you can set these values
SARS-CoV-2-GISAID-freebayes -f gisaid_Asia_08_03_2020.fasta -g covid19-refseq.fasta -t 10    # call variants, provide full path to files
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Asia.50                   # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Asia.50                     # 50 bp sliding window

# joining pi with Tajima's D values and calculate Tajima's D confidence intervals
awk '{print $1"\t"$3"\t"$4"\t"$5}' Asia.50.windowed.pi > Asia.50.subset.windowed.pi
join -j 2 <(sort -k2 Asia.50.subset.windowed.pi) <(sort -k2 Asia.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Asia.50.pi.D && rm joined Asia.50.subset.windowed.pi
pi-tajima -f Asia.50.pi.D -v merged.GISAID.AF.vcf                                            # Run pi-tajima 
cd ..

### Europe                                 
mkdir GISAID_Europe && cd GISAID_Europe
wget -O gisaid_Europe_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b507b027e055bf2df9/display?to_ext=fasta.gz && gunzip gisaid_Europe_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                                       # check if you can set these values
SARS-CoV-2-GISAID-freebayes -f gisaid_Europe_08_03_2020.fasta -g covid19-refseq.fasta -t 10  # call variants, provide full path to files
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Europe.50                 # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Europe.50                   # 50 bp sliding window

# joining pi with Tajima's D values and calculate Tajima's D confidence intervals
awk '{print $1"\t"$3"\t"$4"\t"$5}' Europe.50.windowed.pi > Europe.50.subset.windowed.pi
join -j 2 <(sort -k2 Europe.50.subset.windowed.pi) <(sort -k2 Europe.50000.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Europe.50.pi.D && rm joined Europe.50.subset.windowed.pi
pi-tajima.sh -f Europe.50.pi.D -v merged.GISAID.AF.vcf                                       # Run pi-tajima 
cd ..

### North_America                                
mkdir GISAID_North_America && cd GISAID_North_America  
wget -O gisaid_North_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ee919645a4a97d76/display?to_ext=fasta.gz && gunzip gisaid_North_America_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                                       # check if you can set these values
SARS-CoV-2-GISAID-freebayes -f gisaid_North_America_08_03_2020.fasta -g covid19-refseq.fasta -t 10   # call variants, provide full path to files
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out North_America.50          # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out North_America.50            # 50 bp sliding window

# joining pi with Tajima's D values and calculate Tajima's D confidence intervals
awk '{print $1"\t"$3"\t"$4"\t"$5}' North_America.50.windowed.pi > North_America.50.subset.windowed.pi
join -j 2 <(sort -k2 North_America.50.subset.windowed.pi) <(sort -k2 North_America.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > North_America.50.pi.D && rm joined North_America.50.subset.windowed.pi
pi-tajima -f North_America.50.pi.D -v merged.GISAID.AF.vcf                                   # Run pi-tajima 
cd ..

### South_America                                
mkdir GISAID_South_America && cd GISAID_South_America
wget -O gisaid_South_America_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5134c7103a63c1db1/display?to_ext=fasta.gz && gunzip gisaid_South_America_08_03_2020.fasta.gz 
ulimit -n 1000000 && ulimit -s 1000000                                                       # check if you can set these values
SARS-CoV-2-GISAID-freebayes -f gisaid_South_America_08_03_2020.fasta -g covid19-refseq.fasta -t 10   # call variants, provide full path to files
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out South_America.50          # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out South_America.50            # 50 bp sliding window

# joining pi with Tajima's D values and calculate Tajima's D confidence intervals
awk '{print $1"\t"$3"\t"$4"\t"$5}' South_America.50.windowed.pi > South_America.50.subset.windowed.pi
join -j 2 <(sort -k2 South_America.50.subset.windowed.pi) <(sort -k2 South_America.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > South_America.50.pi.D && rm joined South_America.50.subset.windowed.pi
pi-tajima -f South_America.50.pi.D -v merged.GISAID.AF.vcf                                   # Run pi-tajima 
cd ..

### Oceania                               
mkdir GISAID_Oceania && cd GISAID_Oceania
wget -O gisaid_Oceania_08_03_2020.fasta.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5ad62fc70fed0a55b/display?to_ext=fasta.gz && gunzip gisaid_Oceania_08_03_2020.fasta.gz
ulimit -n 1000000 && ulimit -s 1000000                                                       # check if you can set these values
SARS-CoV-2-GISAID-freebayes -f gisaid_Oceania_08_03_2020.fasta -g covid19-refseq.fasta -t 10    # call variants, provide full path to files
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out Oceania.50                # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out Oceania.50                  # 50 bp sliding window

# joining pi with Tajima's D values and calculate Tajima's D confidence intervals
awk '{print $1"\t"$3"\t"$4"\t"$5}' Oceania.50.windowed.pi > Oceania.50.subset.windowed.pi
join -j 2 <(sort -k2 Oceania.50.subset.windowed.pi) <(sort -k2 Oceania.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > Oceania.50.pi.D && rm joined Oceania.50.subset.windowed.pi
pi-tajima.sh -f Oceania.50.pi.D -v merged.GISAID.AF.vcf                                      # Run pi-tajima 
```
Files ended in "50.pi.D" contains π and Tajima's D values by bin and can be plotted in any sofware. 

check "postprocessing_pi_D_output_files"  subdirectory, containing:
- bins outside 2.5 and 97.5% confidence intervals are outliers for further study (check bins_2.5%_confidencce.bed and bins_2.5%_confidencce.bed files). 
- pi_tajima.pdf plot, for a visual inspection of these regions.
- 2.5_CI_confidence.recode.vcf.gz and 97.5_CI_confidence.recode.vcf.gz, gzipped vcf files containing variants falling in these regions
