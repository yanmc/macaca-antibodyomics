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
		percent_file = "%s_subject_percent.txt" %header
		writer = csv.writer(open(percent_file,"wt"),delimiter = "\t")
		writer.writerow(['subject_id','total_number','percentage(%)'])

		germ_ID=[]
		unique_ID=[]
		f = open(infile,'rU')
		for line in f:
			line = line.split('\t')
			#line[2] = line[2]
			germ_ID.append(line[2])
		f.close

		for n in germ_ID:
			flag=1
			for i in range(0,len(unique_ID)):
				if n == unique_ID[i][0]:
					unique_ID[i][1] += 1
					flag=0
					break
			if flag==1:
				unique_ID.append([n,1])

		for n in unique_ID:
			writer.writerow([str(n[0]),str(n[1]),str(float(n[1])/float(len(germ_ID))*100)])

		print "******************* Done with file %s, %d occurances ****************" %(infile,len(unique_ID))
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

