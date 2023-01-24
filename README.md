# Rare Variant Burden Testing

*This is a work in progress pipeline to showcase the scripts I have been writing and working on*

Scripts for running QC on exome data and harmonizing with publically available gnomAD controls and then counting rare variants 

A pipeline expanding on TRAPD (Testing Rare vAriants using Public Data) pipeline described in Guo MH, Plummer L, Chan Y-M, Hirschhorn JN, Lippincott MF. Burden testing of rare variants identified through exome sequencing using publicly available control data. American Journal of Human Genetics. 2018. 103(4):522-534.

The original TRAPD scripts can be found here: https://github.com/mhguo1/TRAPD

Harmonization of data when using publically available controls is crucial as even very small systematic differences between datasets can have large effects on results. My lab has been working to expand on the TRAPD methods by incorporating more rigorous QC steps and implementing a new counting algorithm.

For each QC step we generate a list of SNPs to exclude from our cases and the gnomAD controls. Then we combine the lists of SNPs to exclude and generate a new VCF by using the Exclude_snps_and_make_new_vcf script. 

QC1 involves generating a list of positions where 90% of samples have a depth>=10 at that position. We do this using GATK-depthofcoverage and bedtools intersect to generate a VCF of only positions meeting our criteria. 

# QC2

Parses the VEP annotations and excludes SNPs that do not contain a canonical annotation that is HIGH impact and protein coding 

# QC3 

Splits the VCF by X Chromosome 

# QC4 

QC4 involves 2 scripts. The first script creates an exclusion list of SNPs that fail the "PASS" filter in the VCF. The second script finds SNPs in our cases that have no valid alleles after filtering individual genotypes by depth, genotype quality, and allele balance. Individual genotypes need to have a depth>=10, genotype quality>=20, and 0.2< allele balance <0.8 to be considered valid. This follows gnomADs methods. 

# QC5 

Script to generate a list of SNPs where more than 2% of individuals have a bad allele balances (0.2>Allele balance>0.8)

# QC6 

Script to generate a list of SNPs from our cases with an allele frequency greater than an input threshold. Allele frequency is calculated after filtering by individual level quality filters (depth>=10, genotype quality>=20, and 0.2< allele balance <0.8). We also use a list of SNPs with gnomAD POPMAX allele frequencies greater than our input threshold. 

# Make SNP file 

Script to generate a list of genes and the qualifying SNPs used for counting variants. Finds SNPs with a VEP annotation that matches our criteria (Ex. Canonical annotation that is HIGH impact and protein coding)

# Counting Scripts

Counting Cases Script generates a variant count per gene using genotypes that pass the individual level quality filters (depth>=10, genotype quality>=20, and 0.2< allele balance <0.8). The script calculates a heterzygote and homozygote frequency. It then sums the frequencies per gene and calculates a count by mutliplying frequencies by the total number of individuals of the SNP with the smallest number for that gene. This is done to try to account for the 10% of individuals that might not be covered well (Depth>=10) at each SNP. We are currently trying to tackle the issue of Linkage Disequilibrium (LD). The current method we are testing is to count an individual once per gene (even if they have multiple variants of a gene) to avoid overcounting caused by LD. This causes issues when using controls that only offer summary data (like gnomAD).

Counting Controls Script parses the gnomAD VCF and calculates a heterzygote and homozygote frequency. It then sums the frequencies per gene and calculates a count by mutliplying frequencies by the total number of individuals of the SNP with the smallest number for that gene. This is done to try to account for the 10% of individuals that might not be covered well (Depth>=10) at each SNP.     
