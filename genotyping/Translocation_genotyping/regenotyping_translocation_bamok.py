#! /usr/bin/python
# -*- coding: utf-8 -*-
from decimal import *
import re
import sys
import pysam
import gzip
from bam_count_readsok import *
import subprocess
filesok=open(sys.argv[1], "r")
filesok1=filesok.readlines()
saver=open(sys.argv[2], "w")
saver1 = open(sys.argv[3], "w")
samfile = sys.argv[4]
sample = sys.argv[5]
read_length = sys.argv[6]
print sample
path_content_sam = sys.argv[7]

n=0
for i in filesok1:
	if re.findall("start_",i):
		saver1.write("%s" % i)
		i=i.split("\n")[0]
		i= str(i)+" "+"our_genotyper"
		saver.write("%s\n" % i)
	else:
		n=n+1
		print n
		cov_total=0
		cov_mut=0
		line=i.split("\n")[0]
		line1=i.split()
		chr1= line1[0]
		if str(chr1) == "23":
			chr1="X"
		if str(chr1) == "24":
			chr1="Y"
		pos1_ini = line1[1]
		if re.findall('.',pos1_ini):
			pos1_ini = pos1_ini.split(".")[0]
		pos1_fin = line1[2]
		if re.findall('.',pos1_fin):
			pos1_fin =pos1_fin.split(".")[0]
		chr2= line1[3]
		if str(chr2) == "23":
			chr2="X"
		if str(chr2) == "24":
			chr2="Y"
		pos2_ini =line1[4]
		if re.findall('.',pos2_ini):
			pos2_ini = pos2_ini.split(".")[0]
		pos2_fin = line1[5]
		if re.findall('.',pos2_fin):
			pos2_fin =pos2_fin.split(".")[0]
		gt_golden = line1[6]
		gt_merge = line1[7]
		leng=line1[8].split("\n")[0]
		if str(leng) == "NA":
			leng=0
		dif_pos1 = int(pos1_fin)-int(pos1_ini)
		dif_pos2 = int(pos2_fin)-int(pos2_ini)
		pos1_bam_ini_1 = int(pos1_ini)-3
		pos1_bam_fin_1 = int(pos1_ini)+3
		pos1_bam_ini_2 = int(pos1_fin)-3
		pos1_bam_fin_2 = int(pos1_fin)+3
		pos2_bam_ini_1 = int(pos2_ini)-3
		pos2_bam_fin_1 = int(pos2_ini)+3
		pos2_bam_ini_2 = int(pos2_fin)-3
		pos2_bam_fin_2 = int(pos2_fin)+3
		dict_reads={}
		if int(leng) <= 150:
			cov_total,cov_mut,dict_reads= count_reads(chr1,chr2,pos1_bam_ini_1,pos1_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
			cov_total,cov_mut,dict_reads= count_reads(chr2,chr1,pos2_bam_ini_1,pos2_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
		else:
			if int(dif_pos1) <=150 and int(dif_pos2)<=150:
				#DOS POSICIONS
				cov_total,cov_mut,dict_reads= count_reads(chr1,chr2,pos1_bam_ini_1,pos1_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
				cov_total,cov_mut,dict_reads= count_reads(chr2,chr1,pos2_bam_ini_1,pos2_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
			elif int(dif_pos1) > 150 and int(dif_pos2)>150:
				#QUATRE POSICIONS
				cov_total,cov_mut,dict_reads= count_reads(chr1,chr2,pos1_bam_ini_1,pos1_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
				cov_total,cov_mut,dict_reads= count_reads(chr2,chr1,pos2_bam_ini_1,pos2_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
				cov_total,cov_mut,dict_reads= count_reads(chr2,chr1,pos2_bam_ini_2,pos2_bam_fin_2,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
				cov_total,cov_mut,dict_reads= count_reads(chr1,chr2,pos1_bam_ini_2,pos1_bam_fin_2,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
			else:
				if int(dif_pos1) < int(dif_pos2):
					#TRES POSICIONS
					cov_total,cov_mut,dict_reads= count_reads(chr1,chr2,pos1_bam_ini_1,pos1_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
					cov_total,cov_mut,dict_reads= count_reads(chr2,chr1,pos2_bam_ini_1,pos2_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
					cov_total,cov_mut,dict_reads= count_reads(chr2,chr1,pos2_bam_ini_2,pos2_bam_fin_2,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
				else:
					cov_total,cov_mut,dict_reads= count_reads(chr1,chr2,pos1_bam_ini_1,pos1_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
					cov_total,cov_mut,dict_reads= count_reads(chr1,chr2,pos1_bam_ini_2,pos1_bam_fin_2,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
					cov_total,cov_mut,dict_reads= count_reads(chr2,chr1,pos2_bam_ini_1,pos2_bam_fin_1,cov_total,cov_mut,samfile,dict_reads,path_content_sam,sample,read_length)
		heterdown= int(cov_total)*0.2
		heterdown=int(heterdown)
		heterup=int(cov_total)*0.8
		heterup=int(heterup)
		if int(cov_mut) > int(heterdown) and int(cov_mut)< int(heterup):
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+str(line1[6])+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])+" "+"0/1:"+str(cov_mut)+":"+str(cov_total)
			saver.write("%s\n" % linefinal)
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+"0/1"+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])
			saver1.write("%s\n" % linefinal)
		elif int(cov_mut) > int(heterdown) and int(cov_mut) == int(heterup):
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+str(line1[6])+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])+" "+"0/1:"+str(cov_mut)+":"+str(cov_total)
			saver.write("%s\n" % linefinal)
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+"0/1"+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])
			saver1.write("%s\n" % linefinal)
		elif int(cov_mut) == int(heterdown) and int(cov_mut)< int(heterup):
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+str(line1[6])+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])+" "+"0/1:"+str(cov_mut)+":"+str(cov_total)
			saver.write("%s\n" % linefinal)
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+"0/1"+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])
			saver1.write("%s\n" % linefinal)
		elif int(cov_mut)> int(heterup):
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+str(line1[6])+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])+" "+"1/1:"+str(cov_mut)+":"+str(cov_total)
			saver.write("%s\n" % linefinal)
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+"1/1"+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])
			saver1.write("%s\n" % linefinal)
		elif int(cov_mut) < int(heterdown):
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+str(line1[6])+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])+" "+"0/0:"+str(cov_mut)+":"+str(cov_total)
			saver.write("%s\n" % linefinal)
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+"0/0"+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])
			saver1.write("%s\n" % linefinal)
		elif int(cov_mut) == int(heterdown) and int(cov_mut) == int(heterup):
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+str(line1[6])+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])+" "+"0/0"+":"+str(cov_mut)+":"+str(cov_total)
			saver.write("%s\n" % linefinal)
			linefinal=str(line1[0])+" "+str(line1[1])+" "+str(line1[2])+" "+str(line1[3])+" "+str(line1[4])+" "+str(line1[5])+" "+"0/0"+" "+str(line1[7])+" "+str(line1[8])+" "+str(line1[9])+" "+str(line1[10])+" "+str(line1[11])+" "+str(line1[12])+" "+str(line1[13])+" "+str(line1[14])
			saver1.write("%s\n" % linefinal)
