#! /usr/bin/python
# -*- coding: utf-8 -*-
import re
import sys
import decimal



saver=open(sys.argv[4]+str(sys.argv[2])+"_chr"+str(sys.argv[3])+".txt","w")

variants = open(sys.argv[1],"r")
variants = variants.readlines()
var= len(variants)
for k in xrange(1,var):
	chrm= variants[k].split()[0]
	ini = variants[k].split()[1]
	if re.findall(r'\.',ini):
		f = ini.split(".")[0]
		ini=int(f)+1
	leng = variants[k].split()[3]
	if re.findall(r'\.',leng):
		d = leng.split(".")[0]
		leng = int(d)+1
	elif re.findall(r'e',leng):
		x = int(decimal.Decimal(leng))
		leng = x
	win = variants[k].split()[9]
	gt = variants[k].split()[2]
	end = int(ini)+int(leng)
	if int(k) == 1:
		lineok = str(chrm)+" "+str(ini)+" "+str(end) + " "+str(leng)+" " +"0"+" "+str(win)+" "+str(gt)
		saver.write("%s\n" % lineok)
		pos_antic = ini
	else:
		diff=int(ini)-int(pos_antic)
		lineok = str(chrm)+" "+str(ini)+" "+str(end) + " "+str(leng)+" " +str(diff)+" "+str(win)+" "+str(gt)
		pos_antic = ini
		saver.write("%s\n" % lineok)
