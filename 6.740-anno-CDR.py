#!/usr/bin/env python
# encoding: utf-8
"""
6.5-anno-CDR.py -i infile -c chain

Created by Mingchen on 2014-09-01.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

infile:  the input file contain trans region info
chain: heavy or light 

"""
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
import csv, os, sys,glob
from mytools import *
def extracts(row):
	return row[3], row[4], row[7], row[8], row[13], row[14], row[17], row[18], row[23], row[24], row[27], row[28], row[33], row[34]
def getcolumn(ref_file):
	#print type(ref_file)
	the_values = map(extracts, ref_file)
	#the_values = extracts(ref_file)
	return the_values
def code_site(solo_pos):
	dna_solo_pos = 3 * int(solo_pos)
	#print solo_pos, dna_solo_pos
	#print type(solo_pos), type(solo_pos)
	return dna_solo_pos
def code_pos_site(line):
	record = list(line)
	#print type(line),type(record)
	#print line,record
	dna_pos = map(code_site, record)
	return dna_pos

def annotation_dna(the_file, record_prot):
	handle = open(the_file, 'rU')
	the_handle = SeqIO.parse(handle, 'fasta')
	print "Annotating dna"
	record_dna = []
	for seq_index, seq_record in enumerate(the_handle):
		if int(seq_index+1) % 100 == 0:
			print "Annotated %d read..."%int(seq_index+1)
		#print seq_record
		dna_str = seq_record.seq.tostring()
		dna_id = seq_record.id
		ref_file = record_prot
		pro_position = getcolumn(ref_file)
		#print pro_position
		for line_index, line in enumerate(pro_position):
			dna_pos = code_pos_site(line)
			#print dna_pos
			if line_index == seq_index :
				fr1 = dna_str[dna_pos[0]-3 : dna_pos[1]]
				#print "fr1:",fr1, dna_pos[0]-2, dna_pos[1], len(fr1)
				fr_1 = [fr1, dna_pos[0]-2, dna_pos[1], len(fr1)]

				cdr1 = dna_str[dna_pos[2]-3 : dna_pos[3]] 
				#print "cdr1:",cdr1, dna_pos[2]-2, dna_pos[3], len(cdr1)
				cdr_1 = [cdr1, dna_pos[2]-2, dna_pos[3], len(cdr1)]

				fr2 = dna_str[dna_pos[4]-3:dna_pos[5]]
				#print "fr2:",fr2, dna_pos[4]-2, dna_pos[5], len(fr2)
				fr_2 = [fr2, dna_pos[4]-2, dna_pos[5], len(fr2)]

				cdr2 = dna_str[dna_pos[6]-3:dna_pos[7]]
				#print "cdr2:",cdr2, dna_pos[0], dna_pos[1], len(cdr2)
				cdr_2 = [cdr2, dna_pos[6]-2, dna_pos[7], len(cdr2)]

				fr3 = dna_str[dna_pos[8]-3:dna_pos[9]] 
				#print "fr3 :",fr3, dna_pos[2], dna_pos[3], len(fr3)
				fr_3 = [fr3, dna_pos[8]-2, dna_pos[9], len(fr3)]

				cdr3 = dna_str[dna_pos[10]-3:dna_pos[11]]
				#print "cdr3:",cdr3, dna_pos[4], dna_pos[5],len(cdr3)
				cdr_3 = [cdr3, dna_pos[10]-2, dna_pos[11], len(cdr3)]
				
				fr4 = dna_str[dna_pos[12]-3:dna_pos[13]] 
				#print "fr4 :",fr4, dna_pos[2], dna_pos[3], len(fr4)
				fr_4 = [fr4, dna_pos[12]-2, dna_pos[13], len(fr4)]

			else:
				pass

		record = [dna_id, dna_str]
		fr_1.extend(cdr_1)
		fr_1.extend(fr_2)
		fr_1.extend(cdr_2)
		fr_1.extend(fr_3)
		fr_1.extend(cdr_3)
		fr_1.extend(fr_4)
		record.extend(fr_1)
		record_dna.append(record)
	return record_dna

def annotation_heavy_chain(the_file):
	'''
	CDRH1, CDRH2, CDRH3 = 31 to 35, 50 to 65, and 95 to 102
	'''

	#CDR notation of IGHV  protein  
	the_handle = open(the_file, 'rU')
	record_prot = []
	for seq_index, seq_record in enumerate(SeqIO.parse(the_handle, 'fasta')):
		coding_dna = seq_record.seq
		pro_str = str(coding_dna.translate())#priductive and no stop codon
		pro_id = seq_record.id
	
		#annotation cdrh1
		seq_start_H1 = pro_str.find('C', 15, 30) #CDR1 start from 30  #Expected position : 31-35B
		match_H1 = re.search('W(V|I|A)', pro_str[30:50])#CDR2 start from 50 residue
		H1_end = str(match_H1)
		if H1_end != 'None' and seq_start_H1 != -1:
			match_end_H1 = match_H1.start() + 30
			CDR_H1 = pro_str[seq_start_H1+5 : match_end_H1]
			#print "CDRH1:Have C and W(V|I|A) in %s"%pro_id, CDR_H1, seq_start_H1+5+1, match_end_H1, len(CDR_H1)
			cdr_h1 = [CDR_H1, seq_start_H1+5+1, match_end_H1, len(CDR_H1),'very good','have C and W(V|I|A)']
			seq_start_H1 = seq_start_H1+5
			seq_end_H1 = match_end_H1-1
		elif seq_start_H1 != -1 and H1_end == 'None':
			seq_end_H1 = pro_str.find('W',30,50)
			if seq_end_H1 != -1:
				CDR_H1 = pro_str[seq_start_H1+5 : seq_end_H1] #[)
				#print "CDRH1:Have C and W in %s"%pro_id, CDR_H1, seq_start_H1+6, seq_end_H1, len(CDR_H1)
				cdr_h1 = [CDR_H1, seq_start_H1+5+1, seq_end_H1, len(CDR_H1),'good', 'have C and W']
				seq_start_H1 = seq_start_H1+5
				seq_end_H1 = seq_end_H1-1
			elif  seq_end_H1 == -1:
				CDR_H1 = pro_str[seq_start_H1+5 :37]
				#print "CDRH1:Have C in left %s"%pro_id, CDR_H1, seq_start_H1+5, 37, len(CDR_H1)
				cdr_h1 = [CDR_H1,seq_start_H1+5+1, 37, len(CDR_H1),'right bad','have C in left and no W in right']
				seq_start_H1 = seq_start_H1+5
				seq_end_H1 = 36
		elif seq_start_H1 == -1 and H1_end != 'None':
			match_end_H1 = match_H1.start() + 30
			CDR_H1 = pro_str[30 : match_end_H1]
			#print "CDRH1:Have W(V|I|A) in right %s"%pro_id, CDR_H1, 31, match_end_H1, len(CDR_H1)
			cdr_h1 = [CDR_H1, 31, match_end_H1, len(CDR_H1),'left bad','have W(V|I|A) in right and no C in left']
			seq_start_H1 = 30
			seq_end_H1 = match_end_H1-1
		elif seq_start_H1 == -1 and H1_end == 'None':
			seq_end_H1 = pro_str.find('W',30,50)
			if seq_end_H1 != -1:
				CDR_H1 = pro_str[30 : seq_end_H1] #[)
				#print "CDRH1:Have W in right %s"%pro_id, CDR_H1, 31, seq_end_H1, len(CDR_H1)
				cdr_h1 = [CDR_H1, 31, seq_end_H1, len(CDR_H1),'left bad','have W in right and no C in left']
				seq_start_H1 = 30
				seq_end_H1 = seq_end_H1-1
			elif seq_end_H1 == -1:
				CDR_H1 = pro_str[30:37]
				#print "CDRH1:No C and W in %s"%pro_id, CDR_H1,31, 37, len(CDR_H1)#31-37(35B)
				cdr_h1 = [CDR_H1,31, 37, len(CDR_H1),'bad','no C and W']
				seq_start_H1 = 30
				seq_end_H1 = 36

		#annotation cdrh3
		seq_start_H3 = pro_str.find('C',71,110)#Expected position : 95-102 #H2_end  25+12+15+19 =71 
		match_start_H3 = re.search('CAR', pro_str[71:110])
		H3_start = str(match_start_H3)
		if H3_start != 'None' and seq_start_H3 != -1:
			match_start_H3 = match_start_H3.start() + 71
			match_end_H3 = re.search('WG[A-Z]G',pro_str[71:120])
			H3_end = str(match_end_H3)
			if H3_end !=  'None':
				match_end_H3 = match_end_H3.start()+71
				CDR_H3 = pro_str[match_start_H3+3 : match_end_H3]#CDRH3 contain VH ,DH and JH germline gene
				#print "CDRH3:Have CAR in left ,WG[A-Z]G in right %s"%pro_id, CDR_H3, match_start_H3+3+1, match_end_H3, len(CDR_H3)
				cdr_h3 = [CDR_H3, match_start_H3+3+1, match_end_H3, len(CDR_H3),'very good','have CAR in left ,WG[A-Z]G in right']
				seq_start_H3 = match_start_H3+3
				seq_end_H3 = match_end_H3-1
			elif H3_end == 'None':
				endH3 = pro_str.find('W', 71, 120)
				if endH3 != -1:
					CDR_H3 = pro_str[match_start_H3+3 : endH3]#Expected position : 95-102
					#print "CDRH3:Have CAR in left %s"%pro_id, CDR_H3, match_start_H3+3+1, endH3, len(CDR_H3)
					cdr_h3 = [CDR_H3, match_start_H3+3+1, endH3, len(CDR_H3),'good','have CAR in left ,have W in right']
					seq_start_H3 = match_start_H3+3
					seq_end_H3 = endH3-1
				else:
					CDR_H3 = pro_str[match_start_H3+3 : 102]#Expected position : 95-102
					#print "CDRH3:Have CAR in left %s"%pro_id, CDR_H3, match_start_H3+3+1, 102, len(CDR_H3)
					cdr_h3 = [CDR_H3, match_start_H3+3+1, 102, len(CDR_H3),'right bad','have CAR in left ,no WG[A-Z]G and W in right']
					seq_start_H3 = match_start_H3+3
					seq_end_H3 = 101
				
		elif H3_start == 'None' and seq_start_H3 != -1:
			match_end_H3 = re.search('WG[A-Z]G',pro_str[71:120])
			H3_end = str(match_end_H3)
			if H3_end != 'None':
				match_end_H3 = match_end_H3.start()+71
				CDR_H3 = pro_str[seq_start_H3+3 :match_end_H3]#CDRH3 contain VH ,DH and JH germline gene
				#print "CDRH3:Have CXX in left ,WG[A-Z]G in right %s"%pro_id, CDR_H3, seq_start_H3+3+1, match_end_H3, len(CDR_H3)
				cdr_h3 = [CDR_H3, seq_start_H3+3+1, match_end_H3, len(CDR_H3),'good','have CXX in left ,WG[A-Z]G in right']
				seq_start_H3 = seq_start_H3+3
				seq_end_H3 = match_end_H3-1
			elif H3_end == 'None':
				endH3 = pro_str.find('W', 71, 120)
				if endH3 != -1:
					CDR_H3 = pro_str[seq_start_H3+3 : endH3]#Expected position : 95-102
					#print "CDRH3:Have CXX in left %s"%pro_id, CDR_H3, seq_start_H3+3+1, endH3, len(CDR_H3)
					cdr_h3 = [CDR_H3, seq_start_H3+3+1, endH3, len(CDR_H3),'good','have CXX in left ,have W in right']
					seq_start_H3 = seq_start_H3+3
					seq_end_H3 =  endH3-1
				else:
					CDR_H3 = pro_str[seq_start_H3+3 : 102]#Expected position : 95-102
					#print "CDRH3:Have CXX in left %s"%pro_id, CDR_H3, seq_start_H3+3+1, 102, len(CDR_H3)
					cdr_h3 = [CDR_H3, seq_start_H3+3+1, 102, len(CDR_H3),'right bad','have CXX in left ,no WG[A-Z]G in right']
					seq_start_H3 = seq_start_H3+3
					seq_end_H3 = 101
		elif H3_start == 'None' and seq_start_H3 == -1:
			match_H2 = re.search('(K|R)(L|I|V|F|T|A)(T|S|I|A)', pro_str[30:50])
			H2_end = str(match_H2)
			if H2_end !=  'None':
				match_end_H2 = match_H2.start()+60
				seq_start_H3 = match_end_H2-1+33
				match_end_H3 = re.search('WG[A-Z]G',pro_str[71:120])
				H3_end = str(match_end_H3)
				if H3_end != 'None':
					match_end_H3 = match_end_H3.start()+71
					CDR_H3 = pro_str[seq_start_H3 : match_end_H3]#CDRH3 contain VH ,DH and JH germline gene
					#print "Have WG[A-Z]G in right%s"%pro_id, CDR_H3, seq_start_H3+1, match_end_H3, len(CDR_H3)
					cdr_h3 = [CDR_H3, seq_start_H3+1, match_end_H3, len(pro_str), 'good','have (K|R)(L|I|V|F|T|A)(T|S|I|A)in CDRH2 end,have WG[A-Z]G in right']
					seq_start_H3 = seq_start_H3
					seq_end_H3 = match_end_H3-1
				elif H3_end == 'None':
					endH3 = pro_str.find('W', 71, 120)
					if endH3 != -1:
						CDR_H3 = pro_str[seq_start_H3 : endH3]#Expected position : 95-102
						#print "CDRH3:No WG[A-Z]G in right%s"%pro_id, CDR_H3, seq_start_H3+1, 102, len(CDR_H3)
						cdr_h3 = [CDR_H3, seq_start_H3+1, endH3, len(CDR_H3),'right bad','have (K|R)(L|I|V|F|T|A)(T|S|I|A)in CDRH2 end and no WG[A-Z]G in right']
						seq_start_H3 = seq_start_H3
						seq_end_H3 = endH3-1
					else:
						CDR_H3 = pro_str[seq_start_H3 : 102]#Expected position : 95-102
						#print "CDRH3:No WG[A-Z]G in right%s"%pro_id, CDR_H3, seq_start_H3+1, 102, len(CDR_H3)
						cdr_h3 = [CDR_H3, seq_start_H3+1, 102, len(CDR_H3),'right bad','have (K|R)(L|I|V|F|T|A)(T|S|I|A)in CDRH2 end and no WG[A-Z]G in right']
						seq_start_H3 = seq_start_H3
						seq_end_H3 = 101
			elif H2_end ==  'None':
				match_end_H3 = re.search('WG[A-Z]G',pro_str[71:120])
				H3_end = str(match_end_H3)
				if H3_end != 'None':
					match_end_H3 = match_end_H3.start()+71
					CDR_H3 = pro_str[94 : match_end_H3]#CDRH3 contain VH ,DH and JH germline gene
					#print "Have WG[A-Z]G in right%s"%pro_id, CDR_H3, 95, match_end_H3, len(CDR_H3)
					cdr_h3 = [CDR_H3, 95, match_end_H3, len(pro_str), 'left bad','have WG[A-Z]G in right']
					seq_start_H3 = 94
					seq_end_H3 = match_end_H3-1
				elif H3_end == 'None':
					endH3 = pro_str.find('W', 71, 120)
					if endH3 != -1:
						CDR_H3 = pro_str[94 : endH3]#Expected position : 95-102
						#print "CDRH3:No WG[A-Z]G in right%s"%pro_id, CDR_H3, 95, endH3, len(CDR_H3)
						cdr_h3 = [CDR_H3, 95, endH3, len(CDR_H3),'bad','95-102(kabat)']
						seq_start_H3 = 94
						seq_end_H3 = endH3-1
					else:
						CDR_H3 = pro_str[94 : 102]#Expected position : 95-102
						#print "CDRH3:No WG[A-Z]G in right%s"%pro_id, CDR_H3, 95, 102, len(CDR_H3)
						cdr_h3 = [CDR_H3, 95, 102, len(CDR_H3),'bad','95-102(kabat)']
						seq_start_H3 = 94
						seq_end_H3 = 101
	
		#annotation cdrh2
		match_start_H2 = re.search('LEWIG', pro_str[30:50])  #Expected position : 50-65
		H2_start = str(match_start_H2)
		seq_start_H2 = seq_end_H1+15 #start = H1_end+15
		if  H2_start !=  'None':
			match_start_H2 = match_start_H2.start()+30
			if seq_start_H2 == match_start_H2+5:
				match_H2 = re.search('(K|R)(L|I|V|F|T|A)(T|S|I|A)', pro_str[60:80])#25+10+15+16 = 66,25+12+15+19 = 71
				H2_end = str(match_H2)
				seq_end_H2 = seq_start_H3 -33
				if H2_end !=  'None':
					match_end_H2 = match_H2.start()+60
					if match_end_H2-1 == seq_end_H2:
						CDR_h2 = pro_str[match_start_H2+5 : match_end_H2]
						#print "CDRH2:Have LEWIG in left,have K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_h2, match_start_H2+5+1, match_end_H2, len(CDR_h2)
						cdr_h2 = [CDR_h2, match_start_H2+5+1, match_end_H2, len(CDR_h2),'very good','Have LEWIG(15) in left,have K/R-L/I/V/F/T/A-T/S/I/A(-33) in right']
						seq_start_H2 = match_start_H2+5
						seq_end_H2 = match_end_H2-1
					elif match_end_H2-1 != seq_end_H2:
						CDR_H2 = pro_str[match_start_H2+5 : seq_end_H2+1]
						#print "CDRH2:Have LEWIG in left,no K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_H2,match_start_H2+5+1, seq_end_H2+1, len(CDR_H2)
						cdr_h2 = [CDR_H2, match_start_H2+5+1, seq_end_H2+1, len(CDR_H2),'good','Have LEWIG(15) in left,cdr3start-33 in right']
						seq_start_H2 = match_start_H2+5
						seq_end_H2 = seq_end_H2
				elif H2_end == 'None':
					CDR_H2 = pro_str[match_start_H2+5 : seq_end_H2+1]
					#print "CDRH2:Have LEWIG in left,no K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_H2,match_start_H2+5+1, seq_end_H2+1, len(CDR_H2)
					cdr_h2 = [CDR_H2, match_start_H2+5+1, seq_end_H2+1, len(CDR_H2),'good','Have LEWIG(15) in left,no K/R-L/I/V/F/T/A-T/S/I/A in right']
					seq_start_H2 = match_start_H2+5
					seq_end_H2 = seq_end_H2				
			elif  seq_start_H2 != match_start_H2+5:
				match_H2 = re.search('(K|R)(L|I|V|F|T|A)(T|S|I|A)', pro_str[60:80])#25+10+15+16 = 66,25+12+15+19 = 71
				H2_end = str(match_H2)
				seq_end_H2 = seq_start_H3 -33
				if H2_end !=  'None': 
					match_end_H2 = match_H2.start()+60
					if match_end_H2-1 == seq_end_H2:
						CDR_h2 = pro_str[seq_start_H2 : match_end_H2]
						#print "CDRH2:Have LEWIG in left,have K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_h2, seq_start_H2+1, match_end_H2, len(CDR_h2)
						cdr_h2 = [CDR_h2, seq_start_H2+1, match_end_H2, len(CDR_h2),'good','cdr1end+15 in left,have K/R-L/I/V/F/T/A-T/S/I/A(-33) in right']
						seq_start_H2 = seq_start_H2
						seq_end_H2 = match_end_H2-1
					elif match_end_H2-1 != seq_end_H2:
						CDR_H2 = pro_str[seq_start_H2 : seq_end_H2+1]
						#print "CDRH2:Have LEWIG in left,no K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_H2, seq_start_H2+1, seq_end_H2+1, len(CDR_H2)
						cdr_h2 = [CDR_H2, seq_start_H2+1, seq_end_H2+1, len(CDR_H2),'good','cdr1end+15 in left,cdr3start-33 in right']
						seq_start_H2 = seq_start_H2
						seq_end_H2 = seq_end_H2
				elif H2_end == 'None':
					CDR_H2 = pro_str[seq_start_H2 : seq_end_H2+1]
					#print "CDRH2:Have LEWIG in left,no K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_H2,seq_start_H2+1, seq_end_H2+1, len(CDR_H2)
					cdr_h2 = [CDR_H2, seq_start_H2+1, seq_end_H2+1, len(CDR_H2),'good','cdr1end+15 in left,no K/R-L/I/V/F/T/A-T/S/I/A in right']
					seq_start_H2 = seq_start_H2
					seq_end_H2 = seq_end_H2
			
		elif H2_start ==  'None' :
			match_H2 = re.search('(K|R)(L|I|V|F|T|A)(T|S|I|A)', pro_str[60:80])#25+10+15+16 = 66,25+12+15+19 = 71
			H2_end = str(match_H2)
			seq_end_H2 = seq_start_H3 -33
			if H2_end != 'None':
				match_end_H2 = match_H2.start()+60
				if match_end_H2-1 == seq_end_H2:
					CDR_h2 = pro_str[seq_start_H2 : match_end_H2]
					#print "CDRH2:Have LEWIG in left,have K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_h2, seq_start_H2+1, match_end_H2, len(CDR_h2)
					cdr_h2 = [CDR_h2, seq_start_H2+1, match_end_H2, len(CDR_h2),'good','no LEWIG in left,have K/R-L/I/V/F/T/A-T/S/I/A(-33) in right']
					seq_start_H2 = seq_start_H2
					seq_end_H2 = match_end_H2-1
				elif match_end_H2-1 != seq_end_H2:
					CDR_H2 = pro_str[seq_start_H2 : seq_end_H2+1]
					#print "CDRH2:Have LEWIG in left,no K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_H2,seq_start_H2+1, seq_end_H2+1, len(CDR_H2)
					cdr_h2 = [CDR_H2, seq_start_H2+1, seq_end_H2+1, len(CDR_H2),'good','no LEWIG in left,cdr3start-33 in right']
					seq_start_H2 = seq_start_H2
					seq_end_H2 = seq_end_H2
			elif H2_end == 'None':
				CDR_H2 = pro_str[seq_start_H2 : seq_end_H2+1]
				#print "CDRH2:Have LEWIG in left,no K/R-L/I/V/F/T/A-T/S/I/A in right %s"%pro_id, CDR_H2, seq_start_H2+1, seq_end_H2+1, len(CDR_H2)
				cdr_h2 = [CDR_H2, seq_start_H2+1, seq_end_H2+1, len(CDR_H2),'good','no LEWIG in left,no K/R-L/I/V/F/T/A-T/S/I/A in right']
				seq_start_H2 = seq_start_H2
				seq_end_H2 = seq_end_H2
		FR1 = pro_str[ :seq_start_H1]
		fr1 = [FR1, 1, seq_start_H1, len(FR1)]
		FR2 = pro_str[seq_end_H1+1 : seq_start_H2]
		fr2 = [FR2, seq_end_H1+2, seq_start_H2, len(FR2)]
		FR3 = pro_str[seq_end_H2+1 : seq_start_H3]
		fr3 = [FR3, seq_end_H2+2, seq_start_H3, len(FR3)]
		FR4 = pro_str[seq_end_H3+1: ]
		fr4 = [FR4, seq_end_H3+2, len(pro_str), len(FR4)]
		record = [pro_id, pro_str]
		fr1.extend(cdr_h1)
		fr1.extend(fr2)
		fr1.extend(cdr_h2)
		fr1.extend(fr3)
		fr1.extend(cdr_h3)
		fr1.extend(fr4)
		record.extend(fr1)
		record_prot.append(record)
		#writer.writerow(record_prot)
	print "There are %s reads be annotated..."% str(seq_index+1)
	return record_prot
	
	
def annotation_light_chain(the_file):
	the_handle = open(the_file, 'rU')
	record_prot = []
	for seq_index, seq_record in enumerate(SeqIO.parse(the_handle, 'fasta')):
		coding_dna = seq_record.seq
		pro_str = str(coding_dna.translate())#priductive and no stop codon
		pro_id = seq_record.id
		
		#annotation cdrl1
		seq_start_L1 = pro_str.find('C',0,23)#CDR1 start from 24 residue   #24+10 = 34 residue #Expected position : 24-34
		match_end_L1 = re.search('W(YQ|FQ|LQ|YL)',pro_str[33:50])
		L1_end = str(match_end_L1)
		if L1_end != 'None' and seq_start_L1 != -1:
			match_end_L1 = match_end_L1.start() + 33
			CDR_L1 = pro_str[seq_start_L1+1 : match_end_L1]
			#print "CDRL1:Have C and W(YQ|FQ|LQ|YL) in %s"%pro_id, CDR_L1, seq_start_L1+1+1, match_end_L1, len(CDR_L1)
			cdr_l1 = [CDR_L1, seq_start_L1+1+1, match_end_L1, len(CDR_L1),'very good','have C in left and W(YQ|FQ|LQ|YL)in right']
			seq_start_L1 = seq_start_L1+1
			seq_end_L1 = match_end_L1-1
		elif seq_start_L1 != -1 and L1_end == 'None':
			seq_end_L1 = pro_str.find('W',33,50)
			if seq_end_L1 != -1:
				CDR_L1 = pro_str[seq_start_L1+1 : seq_end_L1] #[)
				#print "CDRL1:Have C and W in %s"%pro_id, CDR_L1, seq_start_L1+1+1, seq_end_L1, len(CDR_L1)
				cdr_l1 = [CDR_L1, seq_start_L1+1+1, seq_end_L1, len(CDR_L1),'good','have C in left and W in right']
				seq_start_L1 = seq_start_L1+1
				seq_end_L1 = seq_end_L1-1
			elif  seq_end_H1 == -1:
				CDR_L1 = pro_str[seq_start_L1+1 :34]
				#print "CDRH1:Have C in left %s"%pro_id, CDR_L1, seq_start_L1+1, 34, len(CDR_L1)
				cdr_l1 = [CDR_L1, seq_start_L1+1+1, 34, len(CDR_L1),'right bad','only have C in left']
				seq_start_L1 = seq_start_L1+1
				seq_end_L1 = 33
		elif L1_end != 'None' and seq_start_L1 == -1 :
			match_end_L1 = match_end_L1.start() + 33
			CDR_L1 = pro_str[23 : match_end_L1]
			#print "CDRL1:Have W(YQ|FQ|LQ|YL) in right %s"%pro_id, CDR_L1, 24, match_end_L1, len(CDR_L1)
			cdr_l1 = [CDR_L1, 24, match_end_L1, len(CDR_L1),'left bad','only have W(YQ|FQ|LQ|YL)in right']
			seq_start_L1 = 23
			seq_end_L1 = match_end_L1-1
		elif seq_start_l1 == -1 and L1_end == 'None':
			seq_end_L1 = pro_str.find('W',33,49)
			if seq_end_L1 != -1:
				CDR_L1 = pro_str[23 : seq_end_L1] #[)
				#print "CDRL1:Have W in right %s"%pro_id, CDR_L1, 24, seq_end_L1, len(CDR_L1)
				cdr_l1 = [CDR_L1, 24, seq_end_L1, len(CDR_L1),'left bad', 'only have W in right']
				seq_start_L1 = 23
				seq_end_L1 = seq_end_L1-1
			elif  seq_end_H1 == -1:
				CDR_L1 = pro_str[23 :34]
				#print "CDRL1:No C and W in %s"%pro_id, CDR_L1, 24, 34, len(CDR_L1)
				cdr_l1 = [CDR_L1,24, 34, len(CDR_L1),'bad','24-34(kabat)']
				seq_start_L1 = 23
				seq_end_L1 = 33
	
		#annotation cdrl2
		match_start_L2 = re.search('IY|VY|IK|IF', pro_str[30:50]) #Expected position : 50-56
		L2_start = str(match_start_L2)
		seq_start_L2 = seq_end_L1+16 #start = L1_end+16
		if L2_start ==  'None':
			seq_end_L2 = seq_start_L2 + 7
			CDR_L2 = pro_str[seq_start_L2 : seq_end_L2]
			#print "CDRL2:No IY|VY|IK|IF in left %s"%pro_id, CDR_L2, seq_start_L2+1, seq_end_L2, len(CDR_L2)
			cdr_l2 = [CDR_L2, seq_start_L2+1, seq_end_L2, len(CDR_L2),'good','no IY|VY|IK|IF(generally) in left,length 16 and 7(always)']
			seq_start_L2 = seq_start_L2
			seq_end_L2 = seq_end_L2-1
		elif L2_start !=  'None':
			match_start_L2 = match_start_L2.start()+30
			if match_start_L2+2 == seq_start_L2:
				seq_end_L2 = match_start_L2+2+7
				CDR_L2 = pro_str[match_start_L2+2 : seq_end_L2]
				#print "CDRL2:Have IY|VY|IK|IF in left %s"%CDR_L2, match_start_L2+2+1, seq_end_L2, len(CDR_L2)
				cdr_l2 = [CDR_L2, match_start_L2+2+1, seq_end_L2, len(CDR_L2),'good','have IY|VY|IK|IF in left,length 7']
				seq_start_L2 = match_start_L2+2
				seq_end_L2 = seq_end_L2-1
			elif match_start_L2+2 != seq_start_L2:
				seq_end_L2 = seq_start_L2 + 7
				CDR_L2 = pro_str[seq_start_L2 : seq_end_L2]
				#print "CDRL2:Have IY|VY|IK|IF in left %s"%CDR_L2, seq_start_L2+1, seq_end_L2, len(CDR_L2)
				cdr_l2 = [CDR_L2, seq_start_L2+1, seq_end_L2, len(CDR_L2),'good','have IY|VY|IK|IF in left,length 16 and 7']
				seq_start_L2 = seq_start_L2
				seq_end_L2 = seq_end_L2-1

		#annotation cdrl3
		match_start_L3 = pro_str.find('C',56,90) #23+10+16+7+33=89, 23+17+16+7+33+11=107 #Expected position : 89-97
		start_3 = seq_end_L2 + 33
		if match_start_L3 != -1 and start_3 == match_start_L3+1:
			match_end_L3 = re.search('FG[A-Z]G',pro_str[90:110])
			L3_end = str(match_end_L3)
			if L3_end != 'None':
				match_end_L3 = match_end_3.start()+90
				CDR_L3 = pro_str[match_start_L3+1 : match_end_L3]
				#print "CDRL3:Have C in left and FG[A-Z]G in right %s"%CDR_L3, match_start_L3+1+1, match_end_L3, len(CDR_L3)
				cdr_l3 = [CDR_L3, match_start_L3+1+1, match_end_L3, len(CDR_L3),'very good','have C(33) in left and FG[A-Z]G in right']
				seq_start_L3 = match_start_L3+1
				seq_end_L3 = match_end_L3-1
			elif L3_end == 'None':
				CDR_L3 = pro_str[match_start_L3+1 : 97]
				#print "CDRL3:Have C in left ,no FG[A-Z]G in right 33%s"%CDR_L3, match_start_L3+1+1, 97, len(CDR_L3)
				cdr_l3 = [CDR_L3, match_start_L3+1+1, 97, len(CDR_L3),'right bad','have C(33) in left ,no FG[A-Z]G in right']
				seq_start_L3 = match_start_L3+1
				seq_end_L3 = 96
		elif match_start_L3 != -1 and start_3 != match_start_L3+1:
			match_end_L3 = re.search('FG[A-Z]G',pro_str[90:110])
			L3_end = str(match_end_L3)
			if L3_end != 'None':
				match_end_L3 = match_end_3.start()+90
				CDR_L3 = pro_str[match_start_L3+1 : match_end_L3]
				#print "CDRL3:Have C in left and FG[A-Z]G in right,!33%s"%CDR_L3, match_start_L3+1+1, match_end_L3, len(CDR_L3)
				cdr_l3 = [CDR_L3, match_start_L3+1+1, match_end_L3, len(CDR_L3),'good','have C in left and FG[A-Z]G in right']
				seq_start_L3 = match_start_L3+1
				seq_end_L3 = match_end_L3-1
			elif L3_end == 'None':
				CDR_L3 = pro_str[match_start_L3+1 : 97]
				#print "CDRL3:Have C in left ,no FG[A-Z]G in right %s"%CDR_L3, match_start_L3+1+1, 97, len(CDR_L3)
				cdr_l3 = [CDR_L3, match_start_L3+1+1, 97, len(CDR_L3),'right bad','have C in left ,no FG[A-Z]G in right']
				seq_start_L3 = match_start_L3+1
				seq_end_L3 = 96
		elif match_start_L3 == -1:
			match_end_L3 = re.search('FG[A-Z]G',pro_str[90:110])
			L3_end = str(match_end_L3)
			if L3_end != 'None':
				match_end_L3 = match_end_3.start()+90
				CDR_L3 = pro_str[start_3 : match_end_L3]
				#print "CDRL3:No C in left and have FG[A-Z]G in right %s"%CDR_L3,start_3+1, match_end_L3, len(CDR_L3)
				cdr_l3 = [CDR_L3, start_3+1, match_end_L3, len(CDR_L3),'good','no C in left and have FG[A-Z]G in right']
				seq_start_L3 = start_3
				seq_end_L3 = match_end_L3-1
			elif L3_end == 'None':
				CDR_L3 = pro_str[start_3 : 97]
				#print "CDRL3:no C in left ,no FG[A-Z]G in right =33%s"%CDR_L3, start_3+1, 97, len(CDR_L3)
				cdr_l3 = [CDR_L3, start_3+1, 97, len(CDR_L3),'right bad','no C in left ,no FG[A-Z]G in right']
				seq_start_L3 = start_3
				seq_end_L3 = 96	
		FR1 = pro_str[ :seq_start_L1]
		fr1 = [FR1, 1, seq_start_L1, len(FR1)]
		FR2 = pro_str[seq_end_L1+1 : seq_start_L2]
		fr2 = [FR2, seq_end_L1+2, seq_start_L2, len(FR2)]
		FR3 = pro_str[seq_end_L2+1 : seq_start_L3]
		fr3 = [FR3, seq_end_L2+2, seq_start_L3, len(FR3)]
		FR4 = pro_str[seq_end_L3+1: ]
		fr4 = [FR4, seq_end_L3+2, len(pro_str), len(FR4)]
		record = [pro_id, pro_str]
		fr1.extend(cdr_l1)
		fr1.extend(fr2)
		fr1.extend(cdr_l2)
		fr1.extend(fr3)
		fr1.extend(cdr_l3)
		fr1.extend(fr4)
		record.extend(fr1)
		record_prot.append(record)
		#writer.writerow(record)
	print "There are %s reads be annotated..."%str(seq_index+1)
	return record_prot

def main():
	#files = "426_411_diff.fasta"
	#chain = "heavy"
	#the_file = 'human_get_align_region_heavy_min_trimmed.fasta'
	#infile = 'test.fasta'
	the_files = glob.glob(files)
	for the_file in the_files:
		print "processing %s"%the_file
		header, tail = os.path.splitext(the_file)
		out_file = "%s_%s_CDR_pro.txt" %(header,chain)
		writer = csv.writer(open(out_file, "wt"), delimiter = "\t")
		writer.writerow(['Reads ID','Protein seq','FR1','F1.Start','F1.End','F1.Length','CDRH1','C1.Start','C1.End','C1.Length','C1.Quality','C1.Note','FR2','F2.Start','F2.End','F2.Length','CDRH2','C2.Start','C2.End','C2.Length','C2.Quality','C2.Note','FR3','F3.Start','F3.End','F3.Length','CDRH3','C3.Start','C3.End','C3.Length','C3.Quality','C3.Note','FR4','F4.Start','F4.End','F4.Length'])
	
		W_file = "%s_%s_CDR_dna.txt" %(header,chain)
		writer_dna = csv.writer(open(W_file, "wt"), delimiter = "\t")
		writer_dna.writerow(['Reads ID','DNA seq','FR1','F1.Start','F1.End','F1.Length','CDRH1','C1.Start','C1.End','C1.Length','FR2','F2.Start','F2.End','F2.Length','CDRH2','C2.Start','C2.End','C2.Length','FR3','F3.Start','F3.End','F3.Length','CDRH3','C3.Start','C3.End','C3.Length','FR4','F4.Start','F4.End','F4.Length'])
		if chain == "heavy":
			print "Doing anno heavy chain..."
			record_prot = annotation_heavy_chain(the_file)
			#print record_prot
			for item in record_prot:
				writer.writerow(item)

		if chain == "light":
			print "Doing anno light chain..."
			record_prot = annotation_light_chain(the_file)
			for item in record_prot:
				writer.writerow(item)

		record_dna = annotation_dna(the_file, record_prot)
		for item in record_dna:	
			writer_dna.writerow(item)

if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)
	
	# get parameters from input
	dict_args = processParas(sys.argv, i="infile", c="chain")
	files, chain = getParas(dict_args, "infile", "chain")
	
	main()
	print "finished"
