#!/usr/bin/env python
# encoding: utf-8
"""
6.01-get-rearrangement-summary-info-and-ID.py -i infile_model -org orgnism

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  the input file name' model from IGBLAST output, eg IGBLAST\*
orgism:  human,mouse or rabbit

"""

import re,csv
import sys,os,glob
from mytools import *

def get_info(reader):
	result = []
	for index, line in enumerate(reader):
		con = str(line)
		con=con.replace('\'','')
		con=con[1:-1] #con=con[1:-1]
		hit = re.findall('# Query',con)
		if hit:
			result.append(line)
		if len(line) >= 7: 
			if line[-2] == 'Yes' or line[-2] == 'No':
				result.append(line)
	return result
	
	
def main():
	infiles = glob.glob(infile_model)
	outfile = open("%s_get_rearrangement_summary.txt"%(orgnism),"wa+")
	writer = csv.writer(outfile,delimiter = "\t")
	
	outfile2 = open("%s_get_productively_ID.txt"%(orgnism),"wa+")
	writer2 = csv.writer(outfile2,delimiter = "\t")

	for infile in infiles:
		print "Processing %s " %infile
		reader = csv.reader(open(infile,"rU"),delimiter = "\t")

		hit = get_info(reader)
		#print hit,type(hit)
		for index,line in enumerate(hit):
			#print index,line,hit[index+1]
			writer.writerow(line)
			if len(line) > 1 and line[-2] == "Yes":
				#print line,hit[index],hit[index-1]
				ID = hit[index-1][0].replace('# Query: ','')
				ID_writer = []
				if line[-1] == "+":
					ID_writer.append(ID)
					writer2.writerow(ID_writer)
				elif line[-1] == "-":
					ID = 'reversed|'+ID
					ID_writer.append(ID)
					writer2.writerow(ID_writer)
		
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
	
	
