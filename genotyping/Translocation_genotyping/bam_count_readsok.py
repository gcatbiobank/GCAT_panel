#!/usr/bin/python
# -*- coding: utf-8 -*-
import re
import pysam
import subprocess

def count_reads(chrm,chr2,pos,pos_fin,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length):
	content=subprocess.call(["./samtools-1.14/samtools", 'view', str(samfile), str(chrm)+":"+str(pos)+"-"+str(pos_fin)], stdout=open(str(path_content_sam)+"/sam_content_"+str(sample), "w"))
	all_content=open(str(path_content_sam)+"/sam_content_"+str(sample),"r")
	subprocess.call(["rm", str(path_content_sam)+"/sam_content_"+str(sample)])
	for read in all_content:
		read1=str(read)
		ids = read1.split()[0]
		cigar= read1.split()[5]
		if re.findall ('H',cigar):
			continue
		else:
			if ids in dict_reads:
				continue
			else:
				dict_reads[ids]={}
				chr_pair = read1.split()[6] 
				qual= read1.split()[4]
				if int(qual) <= 20:
					continue
				else:
					if str(chr_pair) =="=" or str(chr_pair) == str(chr2):
						#print chrm
						if str(cigar) == read_length:
							if str(chr_pair) == str(chr2):	
								cov_total=cov_total+1
								cov_mut=cov_mut+1
							else:
								cov_total=cov_total+1
								#print cov_total
						else:
							cov_total=cov_total+1
							cov_mut=cov_mut+1
					#else:
						#print read1
	return(cov_total,cov_mut,dict_reads)
