#!/usr/bin/env python
# encoding: utf-8
"""
get_both_VJ_ID.py -v infile_a -j infile_b -o outfile 

Created by Mingchen on 2014-10-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile_a:  the input IG*V file
infile_b:  the input IG*J file
outfile:  the outfile

"""

import os, sys
from mytools import *

def convert_queryID_to_set(infile):
	query_ID = []
	for line in infile:
		line = line.split('\t')
		query_ID.append(line[1])
	ID_set = set(query_ID)
	return ID_set

def main():
	v = open(infile_a,"rU")
	j = open(infile_b,"rU")
	writer = open(outfile,"wa+")

	ID_set_a = convert_queryID_to_set(v)
	ID_set_b = convert_queryID_to_set(j)
	vj_ID_set = ID_set_a.intersection(ID_set_b)
	for ID in vj_ID_set:
		writer.write(ID + '\n')
	print("There are %d reads have both V and J"%len(vj_ID_set))


if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 7 :
		print __doc__
		sys.exit(0)
	
	# get parameters from input
	dict_args = processParas(sys.argv, v="infile_a", j="infile_b", o="outfile")
	infile_a,infile_b,outfile = getParas(dict_args, "infile_a", "infile_b", "outfile")
	
	main()
