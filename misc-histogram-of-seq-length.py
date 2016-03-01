#!/usr/bin/env python
# encoding: utf-8
"""
8.0-histogram-of-seq-length.py -i infile -org orgnism

Created by Mingchen on 2014-12-16.
Copyright (c) 2014 __MyCompanyName__. All rights reserved.

orgnism:  human,mouse or rabbit
infile: trimmed_seq: the file contain trimmed sequence in 6.41-get-trimmed-reads-fasta

"""
from mytools import *
from Bio import SeqIO
import pylab

def main():
	name, suffix = os.path.splitext(infile)
	sizes = [len(rec) for rec in SeqIO.parse(infile, "fasta")]
	print len(sizes), min(sizes), max(sizes)
	#print sizes
	pylab.hist(sizes, bins=20)
	pylab.title("%i orchid sequences\nLengths %i to %i" % (len(sizes),min(sizes),max(sizes)))
	pylab.xlabel("Sequence length (bp)")
	pylab.ylabel("count")
	#pylab.show()
	pylab.savefig('%s_histogram_length.png'%name) #to save the figure to a file (e.g. as a PNG or PDF).
if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 5 :
		print __doc__
		sys.exit(0)

	# get parameters from input
	dict_args = processParas(sys.argv, org="orgnism", i="infile")
	orgnism, infile = getParas(dict_args, "orgnism", "infile")

	main()
	print "finished"