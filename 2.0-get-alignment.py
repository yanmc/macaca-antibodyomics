#!/usr/bin/env python
# encoding: utf-8
"""
2.0-get-alignment.py -i infile_model -org orgnism

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  the input file name' model from IGBLAST output, eg IGBLAST\*
orgism:  human,mouse or rabbit

"""

import re,csv
import sys,os,glob
from mytools_ymc import *

def main():
	infiles = glob.glob(infile_model)
	
	

	for infile in infiles:
		fname, suffix = os.path.splitext(infile)
		count_v,count_d,count_j = 0,0,0 
		c_v,c_d,c_j='','',''
		outfile = open("%s_get_vdj.txt"%(fname),"wa+")
		writer = csv.writer(outfile,delimiter = "\t")
		
		print "Processing %s " %infile
		reader = csv.reader(open(infile,"rU"),delimiter = "\t")

		for line in reader:

			con = str(line)
			con=con.replace('\'','')
			con=con[1:-1]
		
			hit_v=re.match(r'^V.+',con)
			hit_d=re.match(r'^D.+',con)
			hit_j=re.match(r'^J.+',con)
			if hit_v and len(line) == 30:
				a_v=hit_v.group()
				b_v = a_v.split(',')
				if not(c_v == b_v[1]):
					c_v=b_v[1]
					count_v += 1
					writer.writerow(b_v)
			if hit_d and len(line) == 30:
				a_d=hit_d.group()
				b_d = a_d.split(',')
				if not(c_d == b_d[1]):
					c_d=b_d[1]
					count_d += 1
					writer.writerow(b_d)
			if hit_j and len(line) == 30:
				a_j=hit_j.group()
				b_j = a_j.split(',')
				if not(c_j == b_j[1]):
					c_j=b_j[1]
					count_j += 1
					writer.writerow(b_j)
		print "DONE! There are %d align to V, %d align to D, %d align to J." %(count_v,count_d,count_j)


if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)
	
	# get parameters from input
	dict_args = processParas(sys.argv, org="orgnism",i="infile_model")
	orgnism, infile_model = getParas(dict_args, "orgnism", "infile_model")
	
	main()
	print "finished"