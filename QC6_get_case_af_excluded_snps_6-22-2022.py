#!/usr/bin/python
import optparse
import operator
import re
import time
import sys
import gzip
import multiprocessing

#Generate list of SNPs with allele frequency greater than an input threshold 
#Generates allele frequency using individual level quality filters (Depth>=10, Genotype Quality>=20, 0.2<Allele Balance<0.8)

ts=time.time()
#Command Line Arguments
parser = optparse.OptionParser()
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o","--outfile", action="store",dest="outfilename", default="cases_af_excluded_snps.txt")
parser.add_option("-m","--moreinfooutfile", action="store",dest="moreinfooutfile", default="cases_af_excluded_snps_more_info.txt")
parser.add_option("--samplefile", action="store",dest="samplefilename", default="ALL")
parser.add_option("-f","--filterval", action="store", dest="filterval")
parser.add_option("-l", "--lines", action="store",dest="lines", default=5000)
options, args = parser.parse_args()


outfile=open(options.outfilename,'w')
outfile.close()
moreinfooutfile=open(options.moreinfooutfile,'w')
moreinfooutfile.close()

filter_value=float(options.filterval)

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



def get_af(vcfline,samplelist):
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

	hets=[i for i,val in enumerate(gt_list) if ((str(val) in ["0/1", "1/0", "0|1", "1|0"]) and (0.2 < allele_balance[i] < 0.8))]
	hetcarriers=list(set(hets) & set(samplelist))
	homs=[i for i,val in enumerate(gt_list) if str(val) in ["1/1", "1|1"]]
	homcarriers=list(set(homs) & set(samplelist))
	nons=[i for i,val in enumerate(gt_list) if str(val) in ["0/0", "0/0", "0|0", "0|0"]] #, "./0", "0/.", ".|0", "0|.", "./.", ".|."]
	noncarriers=list(set(nons) & set(samplelist))
	total_ind_number=(float(len(noncarriers)+len(homcarriers)+len(hetcarriers)))
	
	ac_out=len(hetcarriers)+(2*len(homcarriers))
	an=(float(2*total_ind_number))
	if an>0:
		af=(float(ac_out/an))
	else:
		af=2
	
	return snpid,af

#Generates list of SNPs with an allele freq greater than an input threshold 
def get_case_af(line):
	line_vcf=line.rstrip().split('\t')
	if line_vcf[0][0]!="#":
		snpid,af=get_af(line_vcf,sampleindices)
		if af>=filter_value:
			print(str(snpid)+""+"AF greater than"+""+str(filter_value)+":"+str(af))
			with open(options.outfilename, "a") as afoutfile:
				afoutfile.write(str(snpid)+'\n')
                # afoutfile.write(str(snpid)+'\t'+str(af)+'\n')
			with open(options.moreinfooutfile, "a") as afmoreinfooutfile:
				afmoreinfooutfile.write(str(snpid)+'\t'+str(af)+'\t'+str(line))

#Using multiprocessing to speed up script
if __name__ == '__main__':
	# Variables
	
	pool = multiprocessing.Pool()
	LINES_PER_PROCESS = int(options.lines)
	if str(options.vcffilename).endswith(".gz") is True:
		with gzip.open(options.vcffilename, "rb") as infile:
			next(iter((pool.imap(get_case_af, infile, LINES_PER_PROCESS))))
			pool.close()
			pool.join()
	else:
		with open(options.vcffilename, "r") as infile:
			next(iter((pool.imap(get_case_af, infile, LINES_PER_PROCESS))))
			pool.close()
			pool.join()
print(time.time()-ts)