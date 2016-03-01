#!/usr/bin/env python
# encoding: utf-8
"""
6.3-get-min-ID-region.py -i infile -idfile id_file

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  the input file contain alignment region info in 6.1-get-trimmed-trans-region
			for D-region, the file in 6.20-get-D-region
orgnism:  human,mouse,rabbit
id_file:  the file contain query ID which cover V and J and productive from 6.62-get-coverVJ-productively-min80-ID
			for D-region, the file in 6.63-get-NA-D-queryID

note: Also be used to get min D-region  
"""
import csv
from mytools import *
def main():	
	fhead, ftail = os.path.splitext(infile)
	outputfile = open('%s-min.txt'%fhead,"wa+")
	writer = csv.writer(outputfile,delimiter = "\t")
	IDs = open(id_file,"rU")
	f = open(infile,"rU")
	reader = csv.reader(f,delimiter = "\t")
	result = []
	infile_dict = dict()
	print "Reading the inputfile..."
	for line in reader:
		infile_dict[line[0]] = line
	print "Reading the ID file..."
	for ID in IDs:
		ID = ID.replace('\n','')
		result.append(infile_dict[ID])
	print "Writing the result..."
	for line in result:
		writer.writerow(line)
	print"There are %d (%.2f%%)reads are min."%(len(result),float(len(result))/len(infile_dict)*100)
	outputfile.close()
	IDs.close()
	f.close()
if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)
	
	# get parameters from input
	dict_args = processParas(sys.argv, i="infile", idfile="id_file")
	infile, id_file = getParas(dict_args, "infile", "id_file")
	
	main()
	print "finished"