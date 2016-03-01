#!/usr/bin/env python
# encoding: utf-8
"""
get-ratio-min.py -i infile_model 

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile_model:  the input file from 3.1 output, eg: \*min80.txt

"""

import os,glob,csv
from mytools import *

def main():
	infiles = glob.glob(infile_model)
	for infile in infiles:
		print "processing %s"%infile
		header,tail = os.path.splitext(infile)
		percent_file = "%s_maturation.txt" %header
		writer = csv.writer(open(percent_file,"wt"),delimiter = "\t")
		writer.writerow(['maturation','total_number','percentage(%)'])

		mismatch=[]
		unique=[]
		f = open(infile,'rU')
		for line in f:
			line = line.split('\t')
			mismatch.append((float(100)-float(line[3])))
		f.close

		for n in mismatch:
			flag = 1
			for i in range(0,len(unique)):
				if n == unique[i][0]:
					unique[i][1] += 1
					flag = 0
					break
			if flag == 1:
				unique.append([n,1])

		for n in unique:
			writer.writerow([str(n[0]),str(n[1]),str(float(n[1])/float(len(mismatch))*100)])
	
		print "***************** Done with file %s, %d occurances ******************" %(infile,len(unique))
	print "Finished"
	
	
if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, i="infile_model")
	infile_model = getParas(dict_args, "infile_model")

	main()
