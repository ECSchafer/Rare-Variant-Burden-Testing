#!/usr/bin/python
import optparse
import operator
import re
import sys
import gzip 
import bisect
import multiprocessing
import time

#Script that splits the VCF to remove X Chromosome 

#Command Line Arguments
parser = optparse.OptionParser()
parser.add_option("-v", "--vcffile", action="store",dest="vcffilename")
parser.add_option("-o", "--outfile", action="store",dest="outfilename", default="outvcf.vcf")
parser.add_option("--xchroutfile", action="store",dest="xchroutfilename", default="xchr_outvcf.vcf")
parser.add_option("--pass", action="store_true", dest="passfilter")
parser.add_option("-l", "--lines", action="store",dest="lines", default=5000)
options, args = parser.parse_args()


ts=time.time()
outfile=open(options.outfilename,'w')
xchroutfile=open(options.xchroutfilename,'w')
vcffile=open(options.vcffilename,'r')
for line in vcffile:
	line_vcf=line.rstrip().split('\t')
	if line_vcf[0][0]=="#":
		outfile.write(str(line))
		xchroutfile.write(str(line))
	else:
		break
vcffile.close()
outfile.close()    
xchroutfile.close()

def write_output(line):
	line_vcf=line.rstrip().split('\t')
	if line_vcf[0][0]!="#" and ("," not in line_vcf[4]):
		snpid=str(line_vcf[0]).lower().replace("chr", "")+":"+str(line_vcf[1])+":"+str(line_vcf[3])+":"+str(line_vcf[4])
		if 'x' in snpid:
			with open(options.xchroutfilename, "a") as xchrvcf:
					xchrvcf.write(str(line))	
		else:
			with open(options.outfilename, "a") as outfile:
					outfile.write(str(line))

#Using multiprocessing to speed up script
if __name__ == '__main__':
	# Variables

	pool = multiprocessing.Pool()
	LINES_PER_PROCESS = int(options.lines)
	with open(options.vcffilename,"r") as infile:
		next(iter((pool.imap(write_output, infile, LINES_PER_PROCESS))))
		pool.close()
		pool.join()
print(time.time()-ts)

