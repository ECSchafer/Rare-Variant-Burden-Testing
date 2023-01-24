#!/usr/bin/python
import optparse
import operator
import re
import sys
import gzip 
import time
import math

#Parses gnomAD VCF and generates variant counts per gene based on list of SNPs 

#Command Line Arguments
parser = optparse.OptionParser()
parser.add_option("-s", "--snpfile", action="store",dest="snpfilename")
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("--pop", action="store",dest="pop", default="ALL")
parser.add_option("-e", "--exclusionsnpfile", action="store",dest="exlusionsnpfilename")
parser.add_option("--snpoutfile","--snpoutfile", action="store",dest="snpoutfilename", default="control_snp_counts.txt") #Output file name
parser.add_option("-o","--outfile", action="store",dest="outfilename", default="control_counts.txt") #Output file name
options, args = parser.parse_args()

#Specify the population to use when parsing gnomAD, default is ALL
if options.pop is not None:
	pops=str(options.pop).split(',')
else:
	pops=["ALL"]

#Reads in a file with SNPs to exclude 
excluded_snps=[]
if options.exlusionsnpfilename is not None:
	excluded_snpfile=open(options.exlusionsnpfilename,'r')
	for snp_line in excluded_snpfile:
		snp=snp_line.rstrip()
		excluded_snps.append(str(snp))
	excluded_snpfile.close()

#Find annotation header
if str(options.vcffilename).endswith(".gz") is True:
	vcffile=gzip.open(options.vcffilename, "rb")
else:
	vcffile=open(options.vcffilename, "r")
csq_found=0
for line_vcf1 in vcffile:
	if line_vcf1[0]=="#":
		if ("ID=CSQ" in line_vcf1) or ("ID=vep" in line_vcf1):
			csq_anno=line_vcf1.rstrip('\n').replace('"', '').strip('>').split("Format: ")[1].split("|")
			break
vcffile.close()     


def find_vep_gene(annot, csq_anno):
	csq_index=csq_anno.index("SYMBOL")
	genename=annot.split("|")[csq_index]
	return genename

#Find out chromosome format
if str(options.vcffilename).endswith(".gz") is True:
	vcffile=gzip.open(options.vcffilename, "rb")
else:
	vcffile=open(options.vcffilename, "r")
chrformat="number"
for line_vcf1 in vcffile:
	line_vcf=line_vcf1.split("\t")
	if "##contig" in line_vcf[0]:
		if "ID=chr" in line_vcf[0]:
			chrformat="chr"
vcffile.close()

def makesnplist(snpfile):
	#Makes a list of SNPs present in the snpfile
	snplist=[]
	#Read in snpfile
	snp_file=open(snpfile, "r")
	
	for line_snp1 in snp_file:
		line_snp=line_snp1.rstrip('\n').split('\t')

		#Find column corresponding to desired snps
		if line_snp[0]!="GENE":
			snplist=snplist+line_snp[1].split(",")
	return set(snplist)
	snp_file.close()

def num_convert(val_in, val_def):
	#Checks if something is int or float; if not, returns default
	try:
		float(val_in) or int(val_in)
		val_out=float(val_in)
	except ValueError:
		val_out=float(val_def)
	return val_out

def sumcount(genesnps, snptable):
	#Calculates a freq of het and hom individuals used for the final gene count
	het_ind_freq_sum=0
	hom_ind_freq_sum=0
	current_smallest_ind_n=None
	for s in range(0, len(genesnps), 1):
		if genesnps[s] in snptable:
			tempsnp=genesnps[s]
			het_ind_freq_sum=het_ind_freq_sum+float(snptable[tempsnp][1])
			hom_ind_freq_sum=hom_ind_freq_sum+float(snptable[tempsnp][2])
			
			#Finds the smallest total number seen per gene
			if current_smallest_ind_n == None and (snptable[tempsnp][3]!=0):	
				current_smallest_ind_n=snptable[tempsnp][3]
			elif snptable[tempsnp][3] <= current_smallest_ind_n:
				current_smallest_ind_n=snptable[tempsnp][3]

	
	return [het_ind_freq_sum, hom_ind_freq_sum, current_smallest_ind_n]


def extractcounts(pops, vcfline):
	#Parses the gnomAD VCF and grabs the counts 
	vcfline=vcfline.lower()
	ac_out=0
	if "ALL" in pops:
		ac_out=num_convert((";"+vcfline).split((";controls_ac="))[1].split(";")[0].split(",")[0],0)
		an=num_convert((";"+vcfline).split((";controls_an="))[1].split(";")[0].split(",")[0],0)
		if ";nhomalt=" in (";"+vcfline):
			hom_out=(";"+vcfline).split((";controls_nhomalt="))[1].split(";")[0].split(",")[0]
		else:
			hom_out=0
		hom_out=num_convert(hom_out, 0)
	elif "ALL" not in pops:
		ac_out=0
		hom_out=0
		for p in range(0, len(pops), 1):
			temp_pop=pops[p].lower()
			an=num_convert((";"+vcfline).split((";controls_an_"+temp_pop+"="))[1].split(";")[0].split(",")[0],0)
			if (";controls_ac_"+temp_pop+"=") in (";"+vcfline):
				ac_out=ac_out+num_convert((";"+vcfline).split((";controls_ac_"+temp_pop+"="))[1].split(";")[0].split(",")[0],0)
			if ";nhomalt=" in (";"+vcfline):
				hom_out=hom_out+num_convert((";"+vcfline).split((";controls_nhomalt_"+temp_pop+"="))[1].split(";")[0].split(",")[0], 0)
	
	#Calculates a frequency for het and hom individuals 
	if float(an) > 0:
		hom_ind_out=float(hom_out)
		het_ind_out=float(float(ac_out)-(hom_ind_out*2))
		total_ind_n=float(float(an)/2)
		af_out=(ac_out/an)
		if total_ind_n>0:
			hom_ind_freq=float(hom_ind_out/total_ind_n)
			het_ind_freq=float(het_ind_out/total_ind_n)

			return [het_ind_freq, hom_ind_freq, total_ind_n, ac_out, an, af_out]

#Make list of all SNPs across all genes present in snpfile
allsnplist=makesnplist(options.snpfilename)

#Make a hashtable with keys as each SNP, and stores a list of indices of carriers for that SNP
count_table={} 


ts=time.time()
snpoutfile=open(options.snpoutfilename, "w")
snpoutfile.write("#SNP\tCONTROL_HET_IND_FREQ\tCONTROL_HOM_IND_FREQ\tCONTROL_SMALLEST_IND_NUMBER\tCONTROL_HET_IND_COUNT\tCONTROL_HOM_IND_COUNT\tCONTROL_AC\tCONTROL_AN\tCONTROL_AF\tGENES\n")
if str(options.vcffilename).endswith(".gz") is True:
	vcffile=gzip.open(options.vcffilename, "rb")
else:
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
			if 'x' in snpid:
				pass
			else:
				counts=extractcounts(pops, line_vcf[7])
				if counts!=None:
					count_table[snpid]=[snpid, counts[0], counts[1], counts[2]]
					if (counts[0]!=0 or counts[1]!=0):
						het_count=float(math.ceil(float(counts[0])*float(counts[2])))
						hom_count=float(math.ceil(float(counts[1])*float(counts[2])))

						snpoutfile.write(str(snpid)+"\t"+str(counts[0])+"\t"+str(counts[1])+"\t"+str(counts[2])+"\t"+str(het_count)+"\t"+str(hom_count)+"\t"+str(counts[3])+"\t"+str(counts[4])+"\t"+str(counts[5])+"\t"+str(annot_gene_list)+'\n')
						print(count_table[snpid])
vcffile.close()
snpoutfile.close()
print(time.time()-ts)
# Generate output counts by multiplying the freq of hets and freq of homs by the smallest total number per gene
outfile=open(options.outfilename, "w")
outfile.write("#GENE\tCONTROL_HET_IND_FREQ\tCONTROL_HOM_IND_FREQ\tCONTROL_SMALLEST_IND_NUMBER\tCONTROL_HET_IND_COUNT\tCONTROL_HOM_IND_COUNT\n")
snpfile=open(options.snpfilename, "r")
for line_s1 in snpfile:
	line_s=line_s1.rstrip('\n').split('\t')
	if line_s[0][0]!="#":
		genesnplist=list(set(line_s[1].split(',')))
		sumcounts=sumcount(genesnplist, count_table)
		if sumcounts[2] is not None:
			het_count=float(math.ceil(float(sumcounts[0])*float(sumcounts[2])))
			hom_count=float(math.ceil(float(sumcounts[1])*float(sumcounts[2])))
			outfile.write(line_s[0]+"\t"+str(sumcounts[0])+"\t"+str(sumcounts[1])+"\t"+str(sumcounts[2])+"\t"+str(het_count)+"\t"+str(hom_count)+'\n')
outfile.close()
print(time.time()-ts)

