#!/usr/bin/env python
# encoding: utf-8
"""
unique-fastq.py

Created by Mingchen on 2015-05-27.
Copyright (c) 2015 __MyCompanyName__. All rights reserved

"""
import glob, os, sys, subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
#from mytools_ymc import *

def unique_fasta(prj_folder):
	handle = "H23RSS.fasta"
	#handle = "%s/1.2-merged-fastq-file/test0513.extendedFrags.fasta"%prj_folder
	reader = SeqIO.parse(handle, "fasta")
	fname, suffix = os.path.splitext(handle)
	writer = open("%s_unique.fasta"%fname,"w")
	handle_dict, handle_dict_unique, dict_unique = {}, {}, {}
	for index, record in enumerate(reader):
		handle_dict[record.id] = record.seq
		handle_dict_unique.setdefault(record.seq, []).append(record.id)
	for seq, ID in handle_dict_unique.items():
		if len(ID) >= 2:
			print len(ID)
		dict_unique["%s_%d"%(ID[0], len(ID))] = seq
	for ID, seq in dict_unique.items():
		seqrecord = SeqRecord(seq, id =ID)
		SeqIO.write(seqrecord, writer, "fasta")
	print "The number of unique reads in fasta file is %d"%len(dict_unique)

def main():
	print "Begin!"
	prj_folder = os.getcwd()
	
	unique_fasta(prj_folder)

if __name__ == '__main__':
	'''
	# get parameters from input
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, p = "project")
	project = getParas(dict_args, "project")
	
	# create 1st and 2nd subfolders
	'''
	#create_folders(os.getcwd())
	
	main()
	print "Finished"
