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

## Task Scenarios:

## I) Colecting Variants (Sequence Read Archive datasets) with SARS-CoV-2-NGS-freebayes

In order to obtain SARS-CoV-2 variants (viral frequency >= 0.5) users need to use SARS-CoV-2-NGS-freebayes and provide:  

- l: path to SRA_list, Sequence read archive accessions of each datasets (SRR prefix list, in tabular format or txt format). As example, see SRA_Accessions_Aug_03_2020.tabular, containing curated SARS-CoV-2 accessions from SRA until August 03, 2020, provided in this repository.
- g: path to SARS-CoV-2 reference in fasta format (covid19-refseq.fasta, provided in this repository)
- a: viral frequency threshold for variant calling (a number between 0-1). As example 0.5, would mean 50% viral frequency variants as threshold. 
- t: number of threads for calculations 

### Execution

Assuming that binaries are in ```/usr/local/bin```, users can do the following:
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

## II) Collecting Variants (GISAID FASTA genomes) with SARS-CoV-2-GISAID-freebayes
- By using freebayes, we can call variants from GISAID SARS-CoV-2 FASTA genomes (contained in a single file) and aggregate these variants into a single VCF for downstream analysis. 

### Execution

In order to obtain SARS-CoV-2 variants from FASTA GISAID genomes, users need to use SARS-CoV-2-GISAID-freebayes and provide:  

- f: GISAID genomes in FASTA format
- g: path to SARS-CoV-2 reference in fasta format (covid19-refseq.fasta, provided in this repository)
- t: number of threads for calculations 

Assuming that binaries are in ```/usr/local/bin```, do the following: 

```
samtools faidx /full/path/to/covid19-refseq.fasta
# Execute the pipeline, providing full path to merged.GISAID.fasta and covid19-refseq.fasta sequences:
ulimit -n 1000000 && ulimit -s 1000000  # check if you can increase stack size and open file limit, see README_ulimit for details.
SARS-CoV-2-GISAID-freebayes -f merged.GISAID.fasta -g covid19-refseq.fasta -t 10 
```
- This operation will obtain joint calls in a single vcf containing all samples (merged.GISAID.AF.vcf). In this matrix, zeros indicate absence of variant and ones indicate the presence of the variant, per sample. Viral frequencies (AF field) were also added. An intermediate file will be also generated: combined_sites.vcf containing just merged variants.  


## III) Nucleotide diversity and Tajima's D test calculation with pi-tajima
- To estimate nucleotide diversity (π) and Tajima's D test, we will employ vcftools program version from Julien Y. Dutheil (accepting --haploid flag) (https://github.com/jydu/vcftools). 
- As an example, we will execute variant calling and vcftools analysis on GISAID FASTA sequences (i.e. merged.GISAID.fasta), using a sliding window of 50 bp across SARS-CoV-2 genome (can be changed). Then, we will we will use the pi-tajima binary to create a plot of pi and Tajima's D values per bin, as depicted here: https://doi.org/10.1016/j.tig.2006.06.005
- An excellent explanation of Tajima's D test can be found here: https://www.youtube.com/watch?v=wiyay4YMq2A. 
- pi-tajima.sh script requires some R libraries installed, as depicted here: https://github.com/cfarkas/SARS-CoV-2-freebayes#10-install-ggplot2-ggrepel-and-vcfr-r-libraries. 

Assuming users previously runned  ```samtools faidx /full/path/to/covid19-refseq.fasta``` and binaries are in ```/usr/local/bin```, do: 

```
ulimit -n 1000000 && ulimit -s 1000000                                                       # check if you can set these values
SARS-CoV-2-GISAID-freebayes -f merged.GISAID.fasta -g covid19-refseq.fasta -t 10             # call variants, provide full path to files
vcftools --vcf merged.GISAID.AF.vcf --window-pi 50 --haploid --out merged.50                 # 50 bp sliding window
vcftools --vcf merged.GISAID.AF.vcf --TajimaD 50 --haploid --out merged.50                   # 50 bp sliding window

# joining pi with Tajima's D values and calculate Tajima's D confidence intervals
awk '{print $1"\t"$3"\t"$4"\t"$5}' merged.50.windowed.pi > merged.50.subset.windowed.pi
join -j 2 <(sort -k2 merged.50.subset.windowed.pi) <(sort -k2 merged.50.Tajima.D) -o 1.2,1.4,2.4> joined && sed -i 's/ /\t/'g joined
sort -k1 -n joined > merged.50.pi.D && rm joined merged.50.subset.windowed.pi
pi-tajima -f merged.50.pi.D -v merged.GISAID.AF.vcf                                          # Run pi-tajima               
```
The file "merged.50.pi.D" contains binned π and Tajima's D values and can be plotted in any sofware. 

check "postprocessing_pi_D_output_files"  subdirectory, containing:
- bins outside 2.5 and 97.5% confidence intervals are outliers for further study (check bins_2.5%_confidencce.bed and bins_2.5%_confidencce.bed files). 
- pi_tajima.pdf plot, for a visual inspection of these regions.
- 2.5_CI_confidence.recode.vcf.gz and 97.5_CI_confidence.recode.vcf.gz, gzipped vcf files containing variants falling in these regions


## IV) Annotate and collect variants per protein on SnpEff-classified VCF variants

- merged.GISAID.AF.vcf files, containing aggregated variants, can be annotated by SnpEff here in Galaxy: https://usegalaxy.org/root?tool_id=toolshed.g2.bx.psu.edu/repos/iuc/snpeff_sars_cov_2/snpeff_sars_cov_2/4.5covid19. These files can be further processed to obtain a summary of variants per viral gene. As an example, we will analyze worldwide GISAID variants until November 30, 2020 as follows: 
- Assuming binaries are in ```/usr/local/bin``` and users previously runned ```samtools faidx covid19-refseq.fasta```:

```
mkdir SnpEff-SARS-CoV-2
cd SnpEff-SARS-CoV-2
wget -O merged.GISAID.Nov-30-2020.vcf.gz https://usegalaxy.org/datasets/bbd44e69cb8906b5509c81bce034b0a2/display?to_ext=gff3.gz
gunzip merged.GISAID.Nov-30-2020.vcf.gz # ~ 17 GB
SnpEff_processing -v SnpEff-Nov-30-2020.GISAID.vcf    
```
In each folder, variants_per_protein subfolder contain variants per protein. Files ending in ".SnpEff" contains parsed variants per consequence and ".counts" contains associated counts, respectively. Also, the script computed aminoacid changes (see SnpEff.AAchanges files). We suggest user-provided vcf files should be processed in a likewise manner, placing the vcf file in a specific folder and executing steps 3) and 4).

## V) More than 100000 GISAID sequences to analyze?
please visit our wiki page here: https://github.com/cfarkas/SARS-CoV-2-freebayes/wiki#i-two-step-full-pipeline
