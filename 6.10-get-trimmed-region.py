#!/usr/bin/env python
# encoding: utf-8
"""
6.1-get-trimmed-region.py -i infile -ref reference_file -org orgnism

Created by Mingchen on 2014-12-11.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  original fasta file
orgnism:  human,mouse or rabbit
referecce_file:  the file contain variable region start and end info eg: human_get_vdj.txt 

note: already trans the reversed reads 
"""
import glob, csv, Bio, re
from mytools import *

def main():		
	outfile = open("%s_get_trimmed_Variable_region.txt"%(orgnism),"wa+")
	writer = csv.writer(outfile,delimiter = "\t")
	
	record_dict = SeqIO.index(infile,  "fasta")
	
	ref_infile = open(reference_file,"rU")
	reader = csv.reader(ref_infile,delimiter = "\t")
	result = ['Query_ID','query_seq','query_length','trimmed_reads','tr_length']
	ID = ''
	count = 0
	count1 = 0
	start, end = int(), int()
	for line in reader:
		query_ID = line[1].replace(' ','')
		nr_query_ID = query_ID.replace('reversed|','')
		if query_ID != ID:
			count1 += 1
			if (count1) % 10000 == 0:
				print "%d reads processed" %(count1)
			query_seq = record_dict.get(nr_query_ID).seq
			if 'reversed' in query_ID:
				count += 1
				if (count) % 1000 == 0:
					print "%d reversed reads processed" %(count)
				a = record_dict.get(nr_query_ID)
				rc = a.reverse_complement(id = query_ID )
				query_seq = rc.seq
			writer.writerow(result)
			ID = query_ID
			start = int(line[8])
			trans_start = int(line[10])
			end = int(line[9])
		else:
			end = int(line[9])
		result[0] = query_ID
		result[1] = query_seq
		result[2] = len(query_seq)
		if start-trans_start < 0:
			
			if (trans_start-1) % 3 == 0:
				result[3] = query_seq[:end]
			else:
				result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
		else:
			
			print "Maybe a wrong reads which not start at first nucle.%s"%query_ID
			#result[3] = query_seq[start-(trans_start-1)-1:end]
			#result[3] = query_seq[(start - 1 + 3 - (trans_start-1)%3 -1):end]
		result[4] = len(result[3])
	writer.writerow(result)
	print 'Find %d seq and %d reversed seq...'%(count1,count)

if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 7 :
		print __doc__
		sys.exit(0)
	
	# get parameters from input
	dict_args = processParas(sys.argv, org="orgnism", i="infile", ref="reference_file")
	orgnism, infile, reference_file = getParas(dict_args, "orgnism", "infile", "reference_file")
	
	main()
	print "finished"