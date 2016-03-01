#!/usr/bin/env python
# encoding: utf-8
"""
6.4-get-trimmed-reads-fasta.py -i infile_model

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  the input file contain trimmed region

"""

import sys,os,glob,csv
from mytools import *

	
def main():
	infiles = glob.glob(infile_model)
	for infile in infiles:
		print "processing %s"%infile
		header,tail = os.path.splitext(infile)
	
		output_handle = open("%s-sequence.fasta"%header, "w")
		reader = csv.reader(open(infile,"rU"),delimiter = "\t")
		for line in reader:
			output_handle.write(">"+line[0]+'\n'+line[3]+'\n')
		output_handle.close()


if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)
	
	# get parameters from input
	dict_args = processParas(sys.argv, i="infile_model",)
	infile_model = getParas(dict_args, "infile_model")
	
	main()
	print "finished"
			
		
