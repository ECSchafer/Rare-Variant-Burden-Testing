#!/usr/bin/python
import optparse
import operator
import re
import sys
import gzip 
import time
import math

#Parses case VCF and generates variant counts per gene using a list of SNPs

#Command Line Arguments
parser = optparse.OptionParser()
parser.add_option("-s", "--snpfile", action="store",dest="snpfilename") #File matching SNPs to genes
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename") #Path to vcf file
parser.add_option("-e", "--exclusionsnpfile", action="store",dest="exlusionsnpfilename")
parser.add_option("--snpoutfile","--snpoutfile", action="store",dest="snpoutfilename", default="case_snp_counts.txt") #Output file name
parser.add_option("-o","--outfile", action="store",dest="outfilename", default="case_counts.txt") #Output file name
parser.add_option("--samplefile", action="store",dest="samplefilename", default="ALL")
options, args = parser.parse_args()

#Creating an output file that has variant counts per SNP 
snpoutfile=open(options.snpoutfilename, "w")
snpoutfile.write("#SNP\tHET_CARRIERS_BEFORE\tHOM_CARRIERS_BEFORE\tNON_CARRIERS_NUMBER\tHET_CARRIERS_COUNTED\tHOM_CARRIERS_COUNTED\tHET_IND_FREQ\tHOM_IND_FREQ\tTOTAL_IND_NUMBER\tGENE_INDEX\tGENE\n")
snpoutfile.close()

#Reads in a file with SNPs to exclude 
excluded_snps=[]
if options.exlusionsnpfilename is not None:
	excluded_snpfile=open(options.exlusionsnpfilename,'r')
	for snp_line in excluded_snpfile:
		snp=snp_line.rstrip()
		excluded_snps.append(str(snp))
	excluded_snpfile.close()

#Find annotation header
vcffile=open(options.vcffilename, "r")
csq_found=0
for line_vcf1 in vcffile:
	if line_vcf1[0]=="#":
		if ("ID=CSQ" in line_vcf1) or ("ID=vep" in line_vcf1):
			csq_anno=line_vcf1.rstrip('\n').replace('"', '').strip('>').split("Format: ")[1].split("|")
			break
vcffile.close()       

#Find info header
vcffile=open(options.vcffilename, "r")
chrformat="number"
for line_vcf1 in vcffile:
	line_vcf=line_vcf1.split("\t")
	if "##contig" in line_vcf[0]:
		if "ID=chr" in line_vcf[0]:
			chrformat="chr"
	elif line_vcf[0]=="#CHROM":
		#This takes the vcf header line and finds the indices corresponding to the individuals present in the sample file
		samplenames=line_vcf[9:]

		#If User doesn't provide sample list, assume all samples in vcf
		if options.samplefilename=="ALL":
			sampleindices=range(0, len(samplenames),1)

		#else, find the indices corresponding to the samples in the user-provided list
		else:
			#Generate sample list
			sample_list=[]
			sample_file=open(options.samplefilename, "r")
			for line_s1 in sample_file:
					sample_list.append(line_s1.rstrip())
			sample_file.close()
			sampleindices=[i for i,val in enumerate(samplenames) if str(val) in sample_list]
		break
vcffile.close()

print("sampleindices is",sampleindices)


def find_vep_gene(annot, csq_anno):
	csq_index=csq_anno.index("SYMBOL")
	genename=annot.split("|")[csq_index]
	return genename


def findcarriers(vcfline,samplelist):

	#Find the column in the genotype field corresponding to the genotypes, genotype quality, depth, allele depth
	gtcol=vcfline[8].split(":").index('GT')
	gqcol=vcfline[8].split(":").index('GQ')
	dpcol=vcfline[8].split(":").index('DP')
	adcol=vcfline[8].split(":").index("AD")

	snpid=str(vcfline[0]).lstrip("chr")+":"+str(vcfline[1])+":"+str(vcfline[3])+":"+str(vcfline[4])
	#Calculate allele balance
	allele_balance=[]
	for i in vcfline[9:]:
		alt_rc=(i.split(':')[adcol].split(",")[1])
		ref_rc=(i.split(':')[adcol].split(",")[0])
		sum_ad=float(alt_rc)+float(ref_rc)
		gt=i.split(':')[gtcol]
		if sum_ad==0:
			ab=0
			allele_balance.append(ab)
		else:
			ab=float(alt_rc)/sum_ad
			allele_balance.append(ab)
	#Grab individuals' genotypes if they have a depth >= 10 and genotype quality >= 20
	gt_list=[]
	for i in vcfline[9:]:
		gt=i.split(':')[gtcol]
		dp=i.split(':')[dpcol]
		gq=i.split(':')[gqcol]
		alt_rc=(i.split(':')[adcol].split(",")[1])
		ref_rc=(i.split(':')[adcol].split(",")[0])
		sum_ad=float(alt_rc)+float(ref_rc)
		if (str(gq)!="."):
			gq_final=float(gq)
		else:
			gq_final=0
		if (str(dp)!="."):
			dp_final=float(dp)
		else:
			dp_final=0
		if (dp_final>=10) and (float(gq_final)>=20):
			gt_list.append(gt)
		else:
			gt=5
			gt_list.append(gt)

	print("length of gt is", len(gt_list))
	print("length of ab is", len(allele_balance))
	#Generate a list of het individuals if they have an allele balance between 0.2 and .8
	hets=[i for i,val in enumerate(gt_list) if ((str(val) in ["0/1", "1/0", "0|1", "1|0"]) and (0.2 < allele_balance[i] < 0.8))]
	print(hets)
	hetcarriers=list(set(hets) & set(samplelist))
	print(hetcarriers)
	af=0
	
	#Generate a list of hom individuals 
	homs=[i for i,val in enumerate(gt_list) if str(val) in ["1/1", "1|1"]]
	homcarriers=list(set(homs) & set(samplelist))
	
	#Generate a list of noncarriers 
	nons=[i for i,val in enumerate(gt_list) if str(val) in ["0/0", "0/0", "0|0", "0|0"]] #, "./0", "0/.", ".|0", "0|.", "./.", ".|."]
	noncarriers=list(set(nons) & set(samplelist))
	
	#Calculate total number of individuals, the number of het individuals, and the number of hom individuals for this SNP
	total_ind_number=(float(len(noncarriers)+len(homcarriers)+len(hetcarriers)))
	het_ind_n = len(hetcarriers)
	hom_ind_n = (float(len(homcarriers)))
	ac_out=len(hets)+(2*len(homs))
	an=(float(2*total_ind_number))
	if an>0:
	#Calculate allele freq, freq of hom individuals, and freq of het individuals (previously used when generating counts)
		af=(float(ac_out/an))
		if total_ind_number>0:
			hom_ind_freq=float(hom_ind_n/total_ind_number)
			het_ind_freq=float(het_ind_n/total_ind_number)
		else:
			hom_ind_freq=0
			het_ind_freq=0
			total_ind_number=0
	else:
		hom_ind_freq=0
		het_ind_freq=0
		total_ind_number=0

	return [het_ind_freq, hom_ind_freq, total_ind_number, af, an, hetcarriers, homcarriers, noncarriers]

	


def makesnplist(snpfile):
	#Makes a list of SNPs present in the snpfile
	snplist=[]
	#Read in snpfile
	snp_file=open(snpfile, "r")
	
	for line_snp1 in snp_file:
		line_snp=line_snp1.rstrip('\n').split('\t')

		#Find column corresponding to desired snps
		if line_snp[0][0]!="#":
			snplist=snplist+line_snp[1].split(",")
	return set(snplist)
	snp_file.close()

def calculatecount(genesnps, snptable):
	#Calculate a freq of het and hom individuals per gene, but only consider an individual once (Avoid double counting individuals with multiple variants per gene)
	total_het_ind_freq=0
	total_hom_ind_freq=0
	current_smallest_ind_n=None
	gene_index=[]
	no_snp_dups=[]
	for s in range(0, len(genesnps), 1):
		if genesnps[s] not in no_snp_dups:
			no_snp_dups.append(genesnps[s])
	for s in no_snp_dups:
		if s in snptable:
			tempsnp=s
			hetcarriers=[]
			homcarriers=[]
			noncarriers=snptable[tempsnp][3]
			for het in snptable[tempsnp][1]:
				if het not in gene_index:
					gene_index.append(het)
					hetcarriers.append(het)
			for hom in snptable[tempsnp][2]:
				if hom not in gene_index:
					gene_index.append(hom)
					homcarriers.append(hom)			
			total_ind_number=(float(len(noncarriers)+len(homcarriers)+len(hetcarriers)))
			het_ind_n = len(hetcarriers)
			hom_ind_n = (float(len(homcarriers)))
			ac_out=len(hetcarriers)+(2*len(homcarriers))
			an=(float(2*total_ind_number))
			if an>0:
				af=(float(ac_out/an))
				if total_ind_number>0:
					hom_ind_freq=float(hom_ind_n/total_ind_number)
					het_ind_freq=float(het_ind_n/total_ind_number)
				else:
					hom_ind_freq=0
					het_ind_freq=0
					total_ind_number=0
			else:
				hom_ind_freq=0
				het_ind_freq=0
				total_ind_number=0
			
			total_het_ind_freq=total_het_ind_freq+het_ind_freq
			total_hom_ind_freq=total_hom_ind_freq+hom_ind_freq

			#Find the smallest total individual number seen per gene 
			if hom_ind_freq>0 or het_ind_freq>0:
				if current_smallest_ind_n == None and (total_ind_number!=0):	
					current_smallest_ind_n=total_ind_number
				elif total_ind_number <= current_smallest_ind_n:
					current_smallest_ind_n=total_ind_number

			with open(options.snpoutfilename, "a") as snpoutfile:
				snpoutfile.write(str(s)+'\t'+str(snptable[tempsnp][1])+'\t'+str(snptable[tempsnp][2])+'\t'+str(snptable[tempsnp][3])+'\t'+str(hetcarriers)+'\t'+str(homcarriers)+'\t'+str(het_ind_freq)+'\t'+str(hom_ind_freq)+'\t'+str(total_ind_number)+'\t'+str(gene_index)+'\t'+str(snptable[tempsnp][4])+'\n')
	
	return [total_het_ind_freq, total_hom_ind_freq, current_smallest_ind_n]

#Make list of all SNPs across all genes present in snpfile
allsnplist=makesnplist(options.snpfilename)

count_table={} 

ts=time.time()
#Read the vcf and generate dictionary with SNPID and list of het, hom, and total individuals seen for that variant
vcffile=open(options.vcffilename, "r")
for line_vcf1 in vcffile:
	line_vcf=line_vcf1.rstrip().split('\t')
	if line_vcf[0][0]!="#" and ("," not in line_vcf[4]):
		snpid=str(line_vcf[0]).lower().replace("chr", "")+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		vcfline=line_vcf[7].replace("vep=", "CSQ=")
		if "CSQ=" in vcfline:
			annots=(";"+vcfline).split(";CSQ=")[1].split(";")[0].split(",")
			annot_gene_list=None
			for i in range(0,len(annots),1):
				current_gene_name=find_vep_gene(annots[i],csq_anno)
				annot_gene_list=str(str(annot_gene_list)+","+str(current_gene_name))
		if snpid in allsnplist and (snpid not in excluded_snps):
			counts=findcarriers(line_vcf,sampleindices)
			if (counts[0]!=0 or counts[1]!=0) and (counts[2]!=0):
				count_table[snpid]=[snpid, counts[5], counts[6], counts[7], annot_gene_list]

vcffile.close()
# snpoutfile.close()
print(time.time()-ts)

# Generate output counts by multiplying the freq of hets and freq of homs by the smallest total number per gene
outfile=open(options.outfilename, "w")
outfile.write("#GENE\tCASE_HET_IND_FREQ\tCASE_HOM_IND_FREQ\tCASE_SMALLEST_IND_NUMBER\tCASE_HET_IND_COUNT\tCASE_HOM_IND_COUNT\n")
snpfile=open(options.snpfilename, "r")
for line_s1 in snpfile:
	line_s=line_s1.rstrip('\n').split('\t')
	if line_s[0][0]!="#":
		genesnplist=list(set(line_s[1].split(',')))
		sumcounts=calculatecount(genesnplist, count_table)

		if sumcounts[2] is not None:
			het_count=float(math.floor(float(sumcounts[0])*float(sumcounts[2])))
			hom_count=float(math.floor(float(sumcounts[1])*float(sumcounts[2])))			
			outfile.write(str(line_s[0])+"\t"+str(sumcounts[0])+"\t"+str(sumcounts[1])+"\t"+str(sumcounts[2])+"\t"+str(het_count)+"\t"+str(hom_count)+'\n')
print(time.time()-ts)
outfile.close()
