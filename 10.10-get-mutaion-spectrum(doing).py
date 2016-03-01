#!/usr/bin/env python
# encoding: utf-8
"""
8.1-Accepted-Replacement-Matrix.py -i infile_model -org orgnism

Created by Mingchen on 2014-01-28.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

orgnism:  human,mouse or rabbit
infile_model: trimmed_seq: the file contain clustalw output file in 7.0-clustal

"""
#from mytools import *
from  Bio import SeqIO
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
import os, glob, csv

def main():
	germ_file = "Human_germline_gene.txt"
	germ_dict = SeqIO.to_dict(SeqIO.parse(open(germ_file,"rU"), "fasta"))
	#print germ_dict
	ref_dict = dict()
	ref_file = "human_get_produ_cover_align_containD_ID_info.txt"
	refer = csv.reader(open(ref_file,'rU'),delimiter = "\t")
	for index,line in enumerate(refer):
		ref_dict[line[0]] = line[1:]
	#print ref_dict

	orgnism = "human"
	infile_model = "human_127_heavy_dna*.aln"
	infiles = glob.glob(infile_model)
	for infile in infiles:
		print "Processing %s " %infile
		name, suffix = os.path.splitext(infile)
		outfile1 = csv.writer(open("%s_info_result_pssm.txt"%(name),"w"),delimiter = "\t")
		# get an alignment object from a Clustalw alignment output
		c_align = AlignIO.read(infile, 'clustal')
		summary_align = AlignInfo.SummaryInfo(c_align)
		consensus = summary_align.dumb_consensus(consensus_alpha = Bio.Alphabet.IUPAC.IUPACAmbiguousDNA)
		my_pssm = summary_align.pos_specific_score_matrix()
		#print help(c_align)
		#print help(summary_align.pos_specific_score_matrix)
		print "The consensus sequence is %s"%consensus
		print my_pssm
		#refer_ID = ref_dict[str(c_align[0].id.replace('_',':'))][0]
		#refer_seq = germ_dict[c_align[0].id.split('&')[0]]
		#print refer_seq
		print len(c_align[0]),len(c_align)
		
		#outfile1.writerow(["Ref_nul","A","C","-","T","G"])
		zs = zip(c_align[0],my_pssm)
		for index,(i,j) in enumerate(zs):
			line = j.values()
			line.insert(0,i)

			if line[0] == 'A':
				line[1] = line[1]-len(c_align)
			elif line[0] == 'C':
				line[2] = line[2]-len(c_align)
			elif line[0] == '-':
				line[3] = line[3]-len(c_align)
			elif line[0] == 'T':
				line[4] = line[4]-len(c_align)
			elif line[0] == 'G':
				line[5] = line[5]-len(c_align)
			outfile1.writerow(line)

#			outfile1.writerow(line)
main()
'''	
if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, org="orgnism", i="infile_model")
	orgnism, infile_model = getParas(dict_args, "orgnism", "infile_model")

	main()
	print "finished"
'''