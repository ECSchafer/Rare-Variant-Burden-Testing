#!/usr/bin/python
import optparse
import operator
import re
import sys
import time
import gzip
import multiprocessing

#Reads in VCF and a list of SNPs to remove and outputs a new VCF

#Options
parser = optparse.OptionParser()
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="outvcf.vcf")
parser.add_option("-m","--removedoutfile", action="store",dest="removedoutfilename", default="removed_lines.txt")
parser.add_option("-s", "--snpfile", action="store",dest="snpfilename")
parser.add_option("-l", "--lines", action="store",dest="lines", default=5000)
options, args = parser.parse_args()

ts=time.time()

#Reads in a file with SNPs to exclude 
excluded_snps=[]
snpfile=open(options.snpfilename,'r')
for snp_line in snpfile:
	snp=snp_line.rstrip()
	excluded_snps.append(str(snp))
snpfile.close()

#Creates an output file and keeps the header
outfile=open(options.outfilename,'w')
if str(options.vcffilename).endswith(".gz") is True:
	vcffile=gzip.open(options.vcffilename, "rb")
else:
	vcffile=open(options.vcffilename, "r")
for line in vcffile:
	line_vcf=line.rstrip().split('\t')
	if line_vcf[0][0]=="#":
		outfile.write(str(line))
	else:
		break
vcffile.close()
outfile.close()

removedoutfile=open(options.removedoutfilename,'w')
removedoutfile.close()


def write_output(line):
	line_vcf=line.rstrip().split('\t')
	if line_vcf[0][0]!="#" and ("," not in line_vcf[4]):
		snpid=str(line_vcf[0]).lower().replace("chr", "")+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		if snpid in excluded_snps:
			with open(options.removedoutfilename, "a") as removed_lines:
				removed_lines.write(str(line))	
		else:
			with open(options.outfilename, "a") as outfile:
				outfile.write(str(line))


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