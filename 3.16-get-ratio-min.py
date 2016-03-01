#!/usr/bin/env python
# encoding: utf-8
"""
get-ratio-min.py -i infile_model -min min_value

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile_model:  the input file from 3.0 output, eg: \*ratio.txt
min_value: the cutoff value, eg: 80

"""

import csv,glob,os
from mytools import *

def main():
	infiles = glob.glob(infile_model)
	for infile in infiles:
		print "Processing %s " %infile
		header, tail = os.path.splitext(infile)
		write_file = "%s_min%d.txt" %(header,min_value)
		writer = csv.writer(open(write_file, "wt"), delimiter = "\t")
		total = 0
		f = open(infile, "rU")
		for line in f:
			line = line.replace('\n','')
			line = line.split('\t')
			if float(line[-1]) >= float(min_value):
				total += 1
				writer.writerow(line)
		print "*****************Done with file %s, %d reads ratio larger than %d ******************" %(infile,total,min_value)
	print "Finished"

if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, i="infile_model", min="min_value")
	infile_model, min_value = getParas(dict_args, "infile_model", "min_value")

	main()
	