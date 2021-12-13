#! /usr/bin/python
# -*- coding: utf-8 -*-
from decimal import *
import re
import sys
import pysam
import gzip
saver=open(sys.argv[4], "w")
samfile = pysam.AlignmentFile(sys.argv[1], "rb")

cov_median=open(sys.argv[3] ,"r")
covmedian_ok = cov_median.readlines()

merge_callers_file = open(sys.argv[2],"r")
merge_callers_fileok = merge_callers_file.readlines()
read_len=sys.argv[5]

dict_gt = {}

for i in covmedian_ok:
	line=i.split("\n")[0]
	line1=i.split()
	chrm= line1[0]
	pos= line1[1]
	cov_final =line1[7]
	if chrm in dict_gt:
		if pos in dict_gt[chrm]:
			dict_gt[chrm][pos].append(cov_final)
		else:
			dict_gt[chrm][pos] = []
			dict_gt[chrm][pos].append(cov_final)
	else:
		dict_gt[chrm] = {}
		dict_gt[chrm][pos] = []
		dict_gt[chrm][pos].append(cov_final)
w=0
for i in merge_callers_fileok:
	if re.findall(r"start_",i):
		j = i.split("start_")[1]
		j1 =j.split()[0]
		i = i.split("\n")[0]
		lineok = str(i) + " "+str(j1)+"_GT_inhouse"+" "+str(j1)+"_coverage_total"+" "+str(j1)+"_coverage_splitreads"	
		saver.write("%s\n" % lineok)
	else:
		cov_mut=0
		cov_total=0
		line=i.split("\n")[0]
		line1=i.split()
		chrm= i.split()[0]
		if str(chrm) == "23":
			chrm="X"
		elif str(chrm) == "24":
			chrm = "Y"
		pos=i.split()[1]
		if re.findall(r"\.", pos):
			m = pos.split(".")[0]
			pos= int(m)+1
		if int(pos) <= 10:
			pos_ini = int(pos)
		else:
			pos_ini= int(pos)-10
		pos_fin= int(pos)+10
		for read in samfile.fetch(str(chrm), int(pos_ini), int(pos_fin)):
			read1=str(read)
			cigar= read1.split()[5]
			qual= read1.split()[4]
			if str(cigar) == read_len:
				cov_total=cov_total+1
			elif re.findall("H",cigar):
				cov_total=cov_total+1
			elif re.findall("D",cigar):
				cov_total=cov_total+1
			elif re.findall("I",cigar):
				cov_total=cov_total+1
			else:
				cov_total=cov_total+1
				cov_mut=cov_mut+1
		if re.findall(r'\.',line1[3]):
			d = line1[3].split(".")[0]
			leng = int(d)+1
			pos2=int(pos)+int(leng)
		else:
			if re.findall(r'\+',line1[3]):
				leng=int(float(str(line1[3])))
				print leng
				print pos
				pos2=int(pos)+int(leng)
				print pos2
			else:
				pos2=int(pos)+int(line1[3])
		pos2_ini=int(pos2)-10
		pos2_fin=int(pos2)+10
		for read2 in samfile.fetch(str(chrm), int(pos2_ini), int(pos2_fin)):
			read3=str(read2)
			cigar= read3.split()[5]
			qual= read3.split()[4]
			if str(cigar) == read_len:
				cov_total=cov_total+1
			elif re.findall("H",cigar):
				cov_total=cov_total+1
			elif re.findall("D",cigar):
				cov_total=cov_total+1
			elif re.findall("I",cigar):
				cov_total=cov_total+1
			else:
				cov_total=cov_total+1
				cov_mut=cov_mut+1
		pos=str(pos)
		if str(chrm) == "X":
			chrm="23"
		elif str(chrm) == "Y":
			chrm="24"
		if pos in dict_gt[chrm]:
			cov_final = dict_gt[chrm][pos][0]
		else:
			print pos
		heterdown= int(cov_final)*0.2 #cov_final
		heterdown=int(heterdown)
		heterup=int(cov_final)*0.8 #cov_final
		heterup=int(heterup)
		win = line1[5]
		if int(cov_mut) > int(heterdown) and int(cov_mut)< int(heterup):
			linefinal = str(line)+" "+"0/1"+" "+str(cov_final)+" "+str(cov_mut)
	#		linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) > int(heterdown) and int(cov_mut) == int(heterup):
			linefinal = str(line)+" "+"0/1"+" "+str(cov_final)+" "+str(cov_mut)
		#	linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) == int(heterdown) and int(cov_mut)< int(heterup):
			linefinal = str(line)+" "+"0/1"+" "+str(cov_final)+" "+str(cov_mut)
	#		linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut)> int(heterup):
			linefinal = str(line)+" "+"1/1"+" "+str(cov_final)+" "+str(cov_mut)
		#	linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"1/1:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) < int(heterdown):
			linefinal = str(line)+" "+"0/0"+" "+str(cov_final)+" "+str(cov_mut)
	#		linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/0:"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		elif int(cov_mut) == int(heterdown) and int(cov_mut) == int(heterup):
			linefinal = str(line)+" "+"0/0"+" "+str(cov_final)+" "+str(cov_mut)
		#	linefinal=str(line1[0])+"\t"+str(line1[1])+"\t"+str(line1[2])+"\t"+str(line1[3])+"\t"+str(line1[4])+"\t"+str(line1[5])+"\t"+str(line1[6])+"\t"+str(line1[7])+"\t"+str(line1[8])+"\t"+"0/0"+":"+str(cov_mut)
			saver.write("%s\n" % linefinal)
		else:
			print line
