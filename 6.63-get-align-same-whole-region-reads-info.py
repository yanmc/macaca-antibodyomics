#!/usr/bin/env python
# encoding: utf-8
"""
6.80-get-align-same-whole-region-reads-info.py -i infile -org orgnism

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  the input file contain rearrengment info
orgism:  human,mouse or rabbit

"""

import re,csv
import sys,os
from mytools import *

def count(ID_pair):
	unique = []
	for n in ID_pair:
		flag=1
		#print n,type(n)
		#n = str(n)
		#print n,type(n)
		for i in range(0,len(unique)):
			if n[1] == unique[i][0]:
				unique[i][1] += 1
				#unique[i][2] = "&".join(unique[i][2],n[0])
				unique[i][2].append(n[0])
				flag=0
				break
		if flag==1:
			unique.append([n[1],1,[n[0]]])
	return unique

def main():

	outfile = open("%s_get-align-same-whole-region-reads-info.txt"%(orgnism),"wa+")
	writer = csv.writer(outfile,delimiter = "\t")
	writer.writerow(["Cover VDJ","V_gene","D_gene","J_gene","reads_IDs","# reads","Percentage"])

	print "Processing %s " %infile
	reader = open(infile,"rU")
	lines = reader.readlines()
	VDJ = []
	for index,line in enumerate(lines):
		line = line.replace('\n','')
		line = line.split("\t")
		#delete the allele info: *01
		#line[1] = line[1].split('-')[0]
		#line[2] = line[2].split('-')[0]
		line[1] = line[1].split(',')[0]
		line[2] = line[2].split(',')[0]
		line[3] = line[3].split(',')[0]
		vd = '&'.join([line[1],line[2],line[3]])
		VDJ.append([line[0],vd])


	unique_VDJ = count(VDJ)

	for n in unique_VDJ:
		writer.writerow([str(n[0]),str(n[0]).split('&')[0],str(n[0]).split('&')[1],str(n[0]).split('&')[2],str(n[2]),str(n[1]),str(float(n[1])/float(len(VDJ))*100)])

if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, org="orgnism",i="infile")
	orgnism, infile = getParas(dict_args, "orgnism", "infile")

	main()
	print "finished"




