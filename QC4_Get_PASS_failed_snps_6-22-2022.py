#!/usr/bin/python
import optparse
import operator
import re
import time
import sys
import gzip
import multiprocessing

#Script to create an exclusion list of all SNPs that fail the PASS filter 

ts=time.time()

#Command Line Arguments
parser = optparse.OptionParser()
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o","--outfile", action="store",dest="outfilename", default="pass_excluded_snps.txt")
parser.add_option("-m","--moreinfoexcludeoutfile", action="store",dest="excludedmoreinfooutfilename", default="snps_to_exclude.more_info.txt")
parser.add_option("-l", "--lines", action="store",dest="lines", default=5000)

options, args = parser.parse_args()


outfile=open(options.outfilename,'w')
outfile.close()
more_info_exclude_outfile=open(options.excludedmoreinfooutfilename,'w')
more_info_exclude_outfile.close()


def get_pass(line):
	line_vcf=line.rstrip().split('\t')
	if line_vcf[0][0]!="#":
		snpid=str(line_vcf[0]).lower().replace("chr", "")+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		if line_vcf[6]!="PASS" and line_vcf[6]!="ExcessHet":
			with open(options.outfilename, "a") as excludedoutfile:
				excludedoutfile.write(str(snpid)+'\n')
			with open(options.excludedmoreinfooutfilename, "a") as excludedmoreinfooutfilename:
				excludedmoreinfooutfilename.write(str(snpid)+'\t'+str(line_vcf[6])+'\n')

#Using multiprocessing to speed up script
if __name__ == '__main__':
	# Variables
	pool = multiprocessing.Pool()
	LINES_PER_PROCESS = int(options.lines)
	if str(options.vcffilename).endswith(".gz") is True:
		with gzip.open(options.vcffilename, "rb") as infile:
			next(iter((pool.imap(get_pass, infile, LINES_PER_PROCESS))))
			pool.close()
			pool.join()
	else:
		with open(options.vcffilename, "r") as infile:
			next(iter((pool.imap(get_pass, infile, LINES_PER_PROCESS))))
			pool.close()
			pool.join()
print(time.time()-ts)



