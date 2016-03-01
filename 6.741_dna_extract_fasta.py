#!/usr/bin/env python
# encoding: utf-8
"""
misc_edit_dist.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc_edit_dist.py -i infile -r referencefile -o outfile


"""
import sys, os, csv, glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#from mytools import *
#def write_fasta_region():
def extracter(infile):
	#infile = "human-get-trimmed-trans-region-min-sequence_heavy_dna.txt"
	fhead, ftail = os.path.splitext(infile)
	bad_outfile = csv.writer(open("bad_%s.txt"%fhead,"w"),delimiter = "\t")
	ref_infile = open(infile,"rU")
	reader = csv.reader(ref_infile,delimiter = "\t")
	region,FR1_list,CDR1_list,FR2_list,CDR2_list,FR3_list,CDR3_list,FR4_list = [],[],[],[],[],[],[],[]
	total, bad = 0, 0
	reader.next()
	
	for index,line in enumerate(reader):
		if '*' in line[1]:
		#print line
			bad_outfile.writerow(line)
		else:
			#print line[29]
			#print type(line[2]), line[2]
			FR1 = SeqRecord(Seq(line[2]),id=line[0],description="%s-%s | %s"%(line[3],line[4],line[5]),name='FR1-'+line[0])
			#print FR1
			FR2 = SeqRecord(Seq(line[10]),id=line[0],description="%s-%s | %s"%(line[11],line[12],line[13]),name='FR2-'+line[0])
			#print FR2
			FR3 = SeqRecord(Seq(line[18]),id=line[0],description="%s-%s | %s"%(line[19],line[20],line[21]),name='FR3-'+line[0])
			#print FR3
			FR4 = SeqRecord(Seq(line[26]),id=line[0],description="%s-%s | %s"%(line[27],line[28],line[29]),name='FR4-'+line[0])
			#print FR4
			CDR1 = SeqRecord(Seq(line[6]),id=line[0],description="%s-%s | %s"%(line[7],line[8],line[9]),name='CDR1-'+line[0])
			
			CDR2 = SeqRecord(Seq(line[14]),id=line[0],description="%s-%s | %s"%(line[15],line[16],line[17]),name='CDR2-'+line[0])
			
			CDR3 = SeqRecord(Seq(line[22]),id=line[0],description="%s-%s | %s"%(line[23],line[24],line[25]),name='CDR3-'+line[0])
			
			FR1_list.append(FR1),CDR1_list.append(CDR1),FR2_list.append(FR2),CDR2_list.append(CDR2),FR3_list.append(FR3),CDR3_list.append(CDR3),FR4_list.append(FR4)
	#print FR1_list
	Zs = zip(['FR1','CDR1','FR2','CDR2','FR3','CDR3','FR4'],[FR1_list,CDR1_list,FR2_list,CDR2_list,FR3_list,CDR3_list,FR4_list])
	for index,(z1,z2) in enumerate(Zs):
		#print index,(z1,z2)
		#for records in [FR1_list,CDR1_list,FR2_list,CDR2_list,FR3_list,CDR3_list,FR4_list]:
	
		outfile = open("%s_%s.fasta"%(fhead,z1),"w")
		SeqIO.write(z2,outfile,"fasta")

def main():
	#the_file = 'human_get_align_region_heavy_min_trimmed.fasta'
	#infile = 'test.fasta'
	the_files = glob.glob("*heavy_CDR_dna.txt")
	for the_file in the_files:
		print "processing %s"%the_file
		extracter(the_file)
main()
