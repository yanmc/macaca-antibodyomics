#!/usr/bin/env python
# encoding: utf-8
"""
count-len-ratio.py -i infile_model

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile_model:  the input file from 2.3 output

"""

import os, sys, glob, re, csv
from mytools_ymc import *
	
def main():
	infiles = glob.glob(infile_model)
	for infile in infiles:
		header, tail = os.path.splitext(infile)
		write_file = "%s_align_ratio.txt" %header
		writer = csv.writer(open(write_file, "wt"), delimiter = "\t")
		reader = csv.reader(open(infile,"rU"),delimiter = "\t")
		print "Processing %s " %infile
		for line in reader:
			alignment_length = int(line[4])-int(line[7])
			ratio = float(alignment_length)/float(line[15])*100
			line.append(ratio)
			writer.writerow(line)
		print "****************** Done with file %s *****************" %infile

if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)
	
	# get parameters from input
	dict_args = processParas(sys.argv, i="infile_model")
	infile_model = getParas(dict_args, "infile_model")
	
	main()

