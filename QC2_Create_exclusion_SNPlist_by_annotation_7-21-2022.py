#!/usr/bin/python
import optparse
import operator
import re
import sys
import gzip 
import bisect
import multiprocessing
import time

#Script that reads in VEP annotations and creates a list of SNPs to exclude for SNPs that do not include an annotation with canonical, protien_coding, HIGH impact variant

#Command Line Arguments
parser = optparse.OptionParser()
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("--excludeoutfile", action="store",dest="excludeoutfilename", default="snps_to_exclude.txt")
parser.add_option("--csq", action="store", dest="csq", default="HIGH")
parser.add_option("--biotype", action="store", dest="biotype", default="protein_coding")
parser.add_option("-l", "--lines", action="store",dest="lines", default=5000)
options, args = parser.parse_args()

##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|MANE|TSL|APPRIS|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">
csq_option=str(options.csq)
biotype_option=str(options.biotype)
ts=time.time()

#Read in vcf header and extract all INFO fields
info_fields=[]
chrformat="number"
if str(options.vcffilename).endswith(".gz") is True:
	vcffile=gzip.open(options.vcffilename, "rb")
else:
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

#Looking for VEP annotations that match the criteria we are looking for
def find_biotype(annots):
	csq_field_biotype_option="BIOTYPE"
	biotype_index=csq_anno.index(csq_field_biotype_option)
	field_value=annots.split("|")[biotype_index]
	out=2
	if str(field_value)==biotype_option:
		out=1
	return out

def find_consequence(annots):
	csq_field_csq_option="IMPACT"
	csq_index=csq_anno.index(csq_field_csq_option)
	field_value=annots.split("|")[csq_index]
	out=2
	if str(field_value)==csq_option:
		out=1
	return out

def canonical_vep(annots):
	canonical_index=csq_anno.index("CANONICAL")
	field_value=annots.split("|")[canonical_index]
	out=2
	if field_value.upper()=="YES":
		out=1
	return out

exclude_outfile=open(options.excludeoutfilename,'w')
exclude_outfile.close()

#Iterates over annotations per SNP to try to find an annotation that is HIGH impact, protein_coding, and canonical
#Writes to exclude file if all annotations do not match the criteria we are looking for 
def write_output(line):
	line_vcf=line.rstrip().split('\t')
	if line_vcf[0][0]!="#" and ("," not in line_vcf[4]):
		snpid=str(line_vcf[0]).lower().replace("chr", "")+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		vcfline=line_vcf[7].replace("vep=", "CSQ=")
		if "CSQ=" in vcfline:
			annots=(";"+vcfline).split(";CSQ=")[1].split(";")[0].split(",")
			keep_var=[]
			for i in range(0, len(annots), 1):
				anno_out=[]
				biotype=find_biotype(annots[i])
				anno_out.append(biotype)
				canonical=canonical_vep(annots[i])
				anno_out.append(canonical)
				consequence_out=find_consequence(annots[i])
				anno_out.append(consequence_out)
				print(biotype)
				print(canonical)
				print(anno_out)
				if 2 in anno_out:
					keep_var.append(2)
				else:
					keep_var.append(1)
			if 1 not in keep_var:
				with open(options.excludeoutfilename, "a") as excludedout:
					excludedout.write(str(snpid)+'\n')

#Using multiprocessing to speed up script
if __name__ == '__main__':
	# Variables

	pool = multiprocessing.Pool()
	LINES_PER_PROCESS = int(options.lines)
	if str(options.vcffilename).endswith(".gz") is True:
		with gzip.open(options.vcffilename, "rb") as infile:
			next(iter((pool.imap(write_output, infile, LINES_PER_PROCESS))))
			pool.close()
			pool.join()
	else:
		with open(options.vcffilename, "r") as infile:
			next(iter((pool.imap(write_output, infile, LINES_PER_PROCESS))))
			pool.close()
			pool.join()
print(time.time()-ts)

