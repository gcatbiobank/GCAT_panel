#! /usr/bin/python
# -*- coding: utf-8 -*-
from decimal import *
import re
import sys
import pysam

filesok=open(sys.argv[1], "r")
filesok1=filesok.readlines()
saver=open(sys.argv[2], "w")
samfile = pysam.AlignmentFile(sys.argv[3], "rb")

for i in filesok1:
	if re.findall("#",i):
		saver.write("%s" % i)
	else:
		cov_total=0
		cov_mut=0
		line=i.split("\n")[0]
		line1=i.split()
		gt=line.split()[9]
		ad=gt.split(":")[1]
		chrm= i.split()[0]
		pos=i.split()[1]
		pos_fin= int(pos)+10
		for read in samfile.fetch(str(chrm), int(pos), int(pos_fin)):
			read1=str(read)
			cigar= read1.split()[5]
			qual= read1.split()[4]
			if int(qual) <= 20:
				continue
			else:
				if str(cigar) == sys.argv[4]:
					cov_total=cov_total+1
				else:
					cov_total=cov_total+1
					cov_mut=cov_mut+1
		heterdown= int(cov_total)*0.2
		heterdown=int(heterdown)
		heterup=int(cov_total)*0.8
		heterup=int(heterup)
		if int(cov_mut) > int(heterdown) and int(cov_mut)< int(heterup):
			linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) > int(heterdown) and int(cov_mut) == int(heterup):
			linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) == int(heterdown) and int(cov_mut)< int(heterup):
			linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut)> int(heterup):
			linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"1/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) < int(heterdown):
			linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/0:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) == int(heterdown) and int(cov_mut) == int(heterup):
			linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/0"+":"+str(cov_mut)
			saver.write("%s\n" % linefinal)
