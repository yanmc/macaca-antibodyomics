#!/usr/bin/env python
# encoding: utf-8
"""
6.61-get-rearr-summary-produc-cover-VJ-info.py -i infile -org orgnism

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

get the productive and cover vj info from the rearrengment info

infile:  the input file contain rearrengment info
orgism:  human,mouse or rabbit

HV
"""

import re,csv
import sys,os
from mytools import *

def choose_best_id(line):
	for i in (1,2,3):
		line[i] = line[i].split(',')
		line[i] = line[i][0]
	result = [line[0],line[1],line[2],line[3]]
	return result
def main():
	#infile = "human_get_rearrangement_summary.txt"
	#orgnism = 'human'
	outfile = open("%s_get_produ_cover_align_containD_info.txt"%(orgnism),"wa+")
	
	writer = csv.writer(outfile,delimiter = "\t")
	
	outfile2 = open("%s_get_produ_cover_align_containD_ID_info.txt"%(orgnism),"wa+")
	
	writer2 = csv.writer(outfile2,delimiter = "\t")

	print "Processing %s " %infile
	reader = open(infile,"rU")
	lines = reader.readlines()
	for index,line in enumerate(lines):
		line = line.replace('\n','')
		line = line.split("\t")
		if len(line) == 8 and   line[-2] == "Yes" and line[-3] == "In-frame":
		#if len(line) == 7 and   line[-2] == "Yes" and line[-3] == "In-frame" and 'IGKV' in line[0]:
		#if len(line) == 7 and   line[-2] == "Yes" and line[-3] == "In-frame" and 'IGLV' in line[0]:
			ID = lines[index-1].replace('# Query: ','')
			ID = ID.replace('\n','')
			result = []
			if line[-1] == "+":
				line.insert(0,ID)
				writer.writerow(line)
				result = choose_best_id(line)
				writer2.writerow(result)
			elif line[-1] == "-":
				ID = 'reversed|'+ID
				line.insert(0,ID)
				writer.writerow(line)
				result = choose_best_id(line)
				writer2.writerow(result)
			else:
				print 'warning'

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




