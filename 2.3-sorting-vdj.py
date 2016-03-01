#!/usr/bin/env python
# encoding: utf-8
"""
2.3-sorting-vdj.py -i infile -org orgnism

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  the input file, orgism_get_vdj_add.txt
orgism:  human,mouse or rabbit...

"""

import os, sys, glob, re, csv
from Bio import SeqIO
from mytools_ymc import *

def sorting_vdj(model,infile):
	result = []
	for index, line in enumerate(open(infile, "rU")):
		match = re.findall(model,line)
		if match != []:
			result.append(line)
	return result

def processer(infile):
	for model in ['HV|VH','KV|VK','LV|VL','HD|DH','HJ|JH','KJ|JK','LJ|JL']:
		
		print "Processing %s ,grep %s..." %(infile,model)
		sort_result = sorting_vdj(model,infile)
		fname, suffix = os.path.splitext(infile)
		write_file = "%s_%s.txt" %(fname,model)
		writer = csv.writer(open(write_file, "wt"), delimiter = "\t")
		for line in sort_result:
			line = line.split()
			writer.writerow(line)
		print "Got %d reads align to %s" %(len(sort_result),model)

def main():
	infiles = glob.glob(infile_model)
	for infile in infiles:
		processer(infile)
		

if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)

			# get parameters from input
	dict_args = processParas(sys.argv, i="infile_model", org="orgnism")
	infile_model, orgnism = getParas(dict_args, "infile_model", "orgnism")

	main()

