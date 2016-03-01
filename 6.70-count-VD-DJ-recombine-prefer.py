#!/usr/bin/env python
# encoding: utf-8
"""
6.6-count-VD-DJ-recombine-prefer.py -i infile -org orgnism

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
			if n == unique[i][0]:
				unique[i][1] += 1
				flag=0
				break
		if flag==1:
			unique.append([n,1])
	return unique

def main():

	outfile = open("%s_count-VD-recombine-prefer.txt"%(orgnism),"wa+")
	writer = csv.writer(outfile,delimiter = "\t")
	writer.writerow(["Cover VD","V_gene","D_gene","# reads","Percentage"])
	
	outfile2 = open("%s_count-DJ-recombine-prefer.txt"%(orgnism),"wa+")
	writer2 = csv.writer(outfile2,delimiter = "\t")
	writer2.writerow(["Cover DJ","D_gene","J_gene","# reads","Percentage"])
	
	
	print "Processing %s " %infile
	reader = open(infile,"rU")
	lines = reader.readlines()
	VD, DJ = [],[]
	for index,line in enumerate(lines):
		line = line.replace('\n','')
		line = line.split("\t")
		#delete the allele info: *01
		#line[1] = line[1].split('-')[0]
		#line[2] = line[2].split('-')[0]
		vd = '&'.join([line[1],line[3]])
		VD.append(vd)

	unique_VD = count(VD)
	for n in unique_VD:
		writer.writerow([str(n[0]),str(n[0]).split('&')[0],str(n[0]).split('&')[1],str(n[1]),str(float(n[1])/float(len(VD))*100)])


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

	'''
		for index,line in enumerate(lines):
			line = line.replace('\n','')
			line = line.split("\t")
			#delete the allele info: *01
			#line[3] = line[3].split('*')[0]
			#line[3] = line[3].split('-')[0]
			#line[2] = line[2].split('-')[0]

			dj = '&'.join([line[2],line[3]])
			DJ.append(dj)
		unique_VD = count(VD)
		unique_DJ = count(DJ)
		for n in unique_VD:
			writer.writerow([str(n[0]),str(n[0]).split('&')[0],str(n[0]).split('&')[1],str(n[1]),str(float(n[1])/float(len(VD))*100)])
	'''


