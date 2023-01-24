#!/usr/bin/python
import optparse
import operator
import re
import sys
import gzip 
import bisect
import time

#Script to generate a list of genes and the qualifying SNPs used for counting variants 

parser = optparse.OptionParser()
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="snpfile.txt")
parser.add_option("--genecolname", action="store", dest="genecolname")
parser.add_option("--impact", action="store", dest="impact", default="HIGH")
parser.add_option("--csq", action="store", dest="csq", default="frameshift_variant")
parser.add_option("--biotype", action="store", dest="biotype", default="protein_coding")
parser.add_option("--genenull", action="store", dest="genenull", default=".,NA")
options, args = parser.parse_args()

impact_option=str(options.impact)
biotype_option=str(options.biotype)
csq_option=str(options.csq)

ts=time.time()

#Read in vcf header and extract all INFO fields
info_fields=[]
chrformat="number"
vcffile=open(options.vcffilename, "r")
for line_vcf1 in vcffile:
	if line_vcf1[0]=="#":
		if "##INFO=<ID=" in line_vcf1:
			temp_field=line_vcf1.split("##INFO=<ID=")[1].split(",")[0]
			info_fields.append(temp_field)
		elif "##contig" in line_vcf1:
					if "ID=chr" in line_vcf1:
							chrformat="chr"
	else:
		break
vcffile.close()

#Read in vcf header to get VEP CSQ fields
vcffile=open(options.vcffilename, "r")
csq_found=0
for line_vcf1 in vcffile:
	if line_vcf1[0]=="#":
		if ("ID=CSQ" in line_vcf1):
			csq_anno=line_vcf1.rstrip('\n').replace('"', '').strip('>').split("Format: ")[1].split("|")
			break
vcffile.close()       

#Looking for VEP annotations that match the criteria we are looking for
def find_biotype(annots):
	csq_field_biotype_option="BIOTYPE"
	biotype_index=csq_anno.index(csq_field_biotype_option)
	field_value=annots.split("|")[biotype_index]
	out=2
	if str(field_value)==biotype_option:
		out=1
	return out

def find_impact(annots):
	impact_field_option="IMPACT"
	impact_index=csq_anno.index(impact_field_option)
	field_value=annots.split("|")[impact_index]
	out=2
	if str(field_value)==impact_option:
		out=1
	return out

def find_csq(annots):
	csq_field_option="Consequence"
	csq_index=csq_anno.index(csq_field_option)
	field_value=annots.split("|")[csq_index]
	out=2
	if str(field_value)!=csq_option:
		out=1
	return out

def canonical_vep(annots):
	canonical_index=csq_anno.index("CANONICAL")
	field_value=annots.split("|")[canonical_index]
	out=2
	if field_value.upper()=="YES":
		out=1
	return out

def find_vep_gene(genecolname, annot, csq_anno):
	csq_index=csq_anno.index(genecolname)
	genename=annot.split("|")[csq_index]
	return genename

#Create empty snptable
snptable={}

#Iterates over annotations per SNP to try to find an annotation that is HIGH impact, protein_coding, and canonical
#Creates a dictionary with gene and the qualifying SNPs 
vcffile=open(options.vcffilename, "r")
for line_vcf1 in vcffile:
	line_vcf=line_vcf1.rstrip().split('\t')
	if line_vcf[0][0]!="#" and ("," not in line_vcf[4]):
		snpid=str(line_vcf[0]).lower().replace("chr", "")+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		vcfline=line_vcf[7]
		if "CSQ=" in vcfline:
			annots=(";"+vcfline).split(";CSQ=")[1].split(";")[0].split(",")
			gene=[]
			for i in range(0, len(annots), 1):
				anno_out=[]
				impact=find_impact(annots[i])
				anno_out.append(impact)
				consequence=find_csq(annots[i])
				anno_out.append(consequence)
				biotype=find_biotype(annots[i])
				anno_out.append(biotype)
				canonical=canonical_vep(annots[i])
				anno_out.append(canonical)
				print(consequence)
				print(biotype)
				print(canonical)
				print(anno_out)
				if 2 not in anno_out:
					gene.append(find_vep_gene(options.genecolname, annots[i], csq_anno))
				else:
					pass
			
			if len(gene)>0:
				for i in range(0, len(gene), 1):
					if gene[i] not in options.genenull.split(","):
						if gene[i] not in snptable:
							snptable[gene[i]]=[gene[i], [snpid]]
						else:
							snptable[gene[i]][1].append(snpid)
			


#Write Output
outfile=open(options.outfilename, "w")
outfile.write("#GENE\tSNPS\n")
for x in snptable:
	if len(x)>0:
			snp_out=','.join(snptable[x][1])
			outfile.write(str(x)+"\t"+snp_out+"\n")
outfile.close()
print(time.time()-ts)
