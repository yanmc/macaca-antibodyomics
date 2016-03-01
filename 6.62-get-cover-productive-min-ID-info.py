#!/usr/bin/env python
# encoding: utf-8
"""
6.62-get-cover-productive-min-ID-info.py -i infile -org orgnism -ref referencefile -c chain

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

get the productive and cover vj info from the rearrengment info

infile:  the input file contain rearrengment info in 6.61
orgism:  human,mouse or rabbit
referencefile: the file contain cover VJ ID in 5.0
chain: heavy, kappa,lambda

"""
from mytools import *

def convert_ID_to_set(infile):
	query_ID = []
	for line in infile:
		#line = line.replace(' ','')
		#line = line.replace('\n','')
		line = line.rstrip()
		query_ID.append(line)
	ID_set = set(query_ID)
	return ID_set
def main():
	#infile = "human_get_produ_cover_align_containD_ID_info.txt"
	#orgnism = 'human'
	#chain ='heavy'
	#referencefile = "human_kappa_cover_VJ_ID.txt"
	
	outfile = open("%s_%s_cover_productive_min_ID.txt"%(orgnism,chain),"w")
	#writer = csv.writer(outfile,delimiter = "\t")
	outfile2 = open("%s_%s_cover_productive_min_info.txt"%(orgnism,chain),"w")
	writer2 = csv.writer(outfile2,delimiter = "\t")
	
	refer = open(referencefile,"rU")
	refer_id_set = convert_ID_to_set(refer)
	#print refer_id_set,type(refer_id_set)
	#if 'M03098:1:000000000-ADRR6:1:2108:27237:14766' in refer_id_set:
	#	print 'yes'
	query_ID = []
	count = 0
	print "Processing %s " %infile
	reader = open(infile,"rU")
	lines = reader.readlines()
	for index,line in enumerate(lines):
		line = line.replace('\n','')
		line = line.replace(' ','')
		line = line.split("\t")
		#print line,type(line)
		#print index
		#line[0] = line[0].replace('\n','')
		line[0] = line[0].rstrip()
		#line[0] = line[0].replace(' ','')
		#print line[0],type(line[0])
		if line[0] in refer_id_set:
			#count += 1
			#print index,count, line
			outfile.write(line[0]+'\n')
			writer2.writerow(line)

if __name__ == '__main__':
#	# check parameters
	if len(sys.argv) < 9 :
		print __doc__
		sys.exit(0)
	
#	# get parameters from input
	dict_args = processParas(sys.argv, i="infile", org="orgnism", ref="referencefile", c="chain")
	infile,orgnism,referencefile,chain = getParas(dict_args, "infile", "orgnism", "referencefile", "chain")
	
	main()