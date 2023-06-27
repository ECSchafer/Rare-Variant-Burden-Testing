# Rare Variant Burden Testing

*This is a work in progress pipeline to showcase the scripts I have been writing and working on*

Scripts for running QC on exome data and harmonizing with publicly available gnomAD controls and then counting rare variants 

A pipeline expanding on TRAPD (Testing Rare vAriants using Public Data) pipeline described in Guo MH, Plummer L, Chan Y-M, Hirschhorn JN, Lippincott MF. Burden testing of rare variants identified through exome sequencing using publicly available control data. American Journal of Human Genetics. 2018. 103(4):522-534.

The original TRAPD scripts can be found here: https://github.com/mhguo1/TRAPD

Harmonization of data when using publically available controls is crucial as even very small systematic differences between datasets can have large effects on results. My lab has been working to expand on the TRAPD methods by incorporating more rigorous QC steps and implementing a new counting algorithm. I have also optimized our pipeline by incorporating multiprocesing in python.

For each QC step we generate a list of SNPs to exclude from our cases and the gnomAD controls. Then we combine the lists of SNPs to exclude and generate a new VCF by using the Exclude_snps_and_make_new_vcf script. 

Once we have our variant counts for each gene, we perform a fischer's exact test to find enrichment of genes in our cases.

We use synoynmous variants to calibrate our QC conditions and see how well we can harmonize our data with gnomAD. Synoynmous variants are thought to be benign and there should be no enrichment of genes in either dataset. We use a QQ plot to test our conditions and we are aiming for a lamda of 1. 

I will show our QQ plots at each QC step as we clean up the data and harmonize our dataset with gnomAD.
 
# QC1

Parses the VEP annotations and excludes SNPs that for example don't have a canonical, protein coding, synoynmous variant annotation (If running calibration step). We also use this to exclude SNPs that don't have a canonical, protein coding, HIGH impact or missense annotation when running the actual analysis 

# QC2

Script to generate a list of SNPs from our cases with an allele frequency greater than an input threshold. Allele frequency is calculated after filtering by individual level quality filters (depth>=10, genotype quality>=20, and 0.2< allele balance <0.8). We also use a list of SNPs with gnomAD POPMAX allele frequencies greater than our input threshold. 

We chose to exclude SNPs with a case allele frequency greater than 0.1. SNPs with high allele frequency in our cases are likely due to sequencing artifacts. 

# QC3 

Splits the VCF by X/Y/MT Chromosome 

# QC4 

QC4 involves generating a list of positions where 90% of samples have a depth>=10 at that position. We do this using GATK-depthofcoverage and bedtools intersect to generate a VCF of only positions meeting our criteria. 

![image](https://github.com/ECSchafer/Rare-Variant-Burden-Testing/assets/123387175/a2490315-8a8f-4e4d-b85e-cd2787b8a8b6)

With just the coverage filter applied we see genes with a significance of as high as -e89. We should be seeing a lamda of 1 and the genes falling close to the line of harmony. We need to apply more filters.

# QC5 

QC5 involves 2 scripts. The first script creates an exclusion list of SNPs that fail the "PASS" filter in the VCF. The second script finds SNPs in our cases that have no valid alleles after filtering individual genotypes by depth, genotype quality, and allele balance. Individual genotypes need to have a depth>=10, genotype quality>=20, and 0.2< allele balance <0.8 to be considered valid. This follows gnomADs methods.

![image](https://github.com/ECSchafer/Rare-Variant-Burden-Testing/assets/123387175/71b0528a-d641-46c5-aba9-a38ec688f9a0)

This filter cleans up things up a bit, but our lamda is still over 1 and we have genes with a significance greater than -e20

# QC6 

Script to generate a list of SNPs where more than 2% of individuals have a bad allele balance (0.2>Allele balance>0.8)

![image](https://github.com/ECSchafer/Rare-Variant-Burden-Testing/assets/123387175/3881ab68-4138-4038-b759-662d1242e035)

With this step we have lowered the significance of the top genes to -e10. We have also filtered out OR4A16 as the top gene.

# QC7

A script that removes variants with a QD greater than 2. This follows one of gnomAD's filter criteria.

![image](https://github.com/ECSchafer/Rare-Variant-Burden-Testing/assets/123387175/3eb6c95c-92be-4f31-bc8a-234b3e126e6f)

This cleans up the graph really well and we see genes falling close to the line of harmony.

# QC8

With this QC step we removed related individuals.

![image](https://github.com/ECSchafer/Rare-Variant-Burden-Testing/assets/123387175/ce8ee2c0-32aa-4e5b-841c-7b5f61e80efe)

Our lamda is 0.98 and there is only one outlier gene.

# Next Step: Run The Rare-Variant Burden Test On Protein Coding HIGH Impact Variants

We have harmonized our data with gnomAD and are very satisfied with our QC filters. We can now apply these QC filters to protein coding HIGH impact and missense variants to find enriched genes in our cases.

We currently have very promising preliminary results.

# Make SNP file 

Script to generate a list of genes and the qualifying SNPs used for counting variants. Finds SNPs with a VEP annotation that matches our criteria (Ex. Canonical annotation that is HIGH impact and protein coding)

# Counting Scripts

Counting Cases Script generates a variant count per gene using genotypes that pass the individual level quality filters (depth>=10, genotype quality>=20, and 0.2< allele balance <0.8). The script calculates a heterzygote and homozygote frequency. It then sums the frequencies per gene and calculates a count by mutliplying frequencies by the total number of individuals of the SNP with the smallest number for that gene. This is done to try to account for the 10% of individuals that might not be covered well (Depth>=10) at each SNP. We are currently trying to tackle the issue of Linkage Disequilibrium (LD). The current method we are testing is to count an individual once per gene (even if they have multiple variants of a gene) to avoid overcounting caused by LD. This is probably not our final method and currently causes issues when using controls that only offer summary data (like gnomAD).

Counting Controls Script parses the gnomAD VCF and calculates a heterzygote and homozygote frequency. It then sums the frequencies per gene and calculates a count by mutliplying frequencies by the total number of individuals of the SNP with the smallest number for that gene. This is done to try to account for the 10% of individuals that might not be covered well (Depth>=10) at each SNP.     
