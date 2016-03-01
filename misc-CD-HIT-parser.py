#!/usr/bin/env python
# encoding: utf-8
"""
misc-CD-HIT-parser.py

Created by Zhenhai Zhang on 2013-02-27.
Copyright (c) 2013 __MyCompanyName__. All rights reserved.
"""

import sys
import os
from mytools import *
import dendropy
MIN_NUM = 20
def get_cluster_id(s):
	return int(s.split(" ")[-1])

def get_abid(s):
	abid = s[s.index(">") +1 : s.index(".")]
	if len(abid) != 8:
		print s, abid
	return abid

def sort_abids(ab_set):
	ab_list = sorted(list(ab_set))
	return ab_list

def get_rep_id(abids):
	abids = sorted(list(abids))
	heads = sorted(list(set([x[0] for x in abids])))
	head = heads[0] 
	if "O" in heads:
		head = "O"
	rep_abids = [x for x in abids if x.startswith(head)]
	return random.choice(rep_abids)
		

def main():
	infile = sys.argv[1]
	reader = open(infile, "rU")
	dict_cluster_abs = dict()
	cid = -1
	for line in reader:
		line = line.strip()  # remove blank chars
		if line.startswith(">Cluster"):
			# this is the head of the cluster
			cid = get_cluster_id(line)
			dict_cluster_abs[cid] = set()
		else:
			# this is a member
			abid = get_abid(line)
			#print "Cluster %d : %s" %(cid, abid)
			dict_cluster_abs[cid].add(abid)
	total = 0
	writer = csv.writer(open("cluster_id-99-95-min%d.txt" %MIN_NUM, "w"), delimiter = sep)
	writer.writerow(["Cluster", "Count", "Abs"])
	
	rep_ids = set()
	for cid, abids in dict_cluster_abs.items():
		if len(abids) >= MIN_NUM:
			writer.writerow([cid, len(abids), ",".join(sort_abids(abids))])
			total += 1
			rep_id = get_rep_id(abids)
			rep_ids.add(rep_id)
	print "%d clusters has at least %d members" %(total, MIN_NUM)
	print "%d representatives..." %len(rep_ids)
	handle = open("cluster_id-99-95-min%d-rep.fa" %MIN_NUM, "w")
	for ab in SeqIO.parse(open("G01.fa", "rU"), "fasta"):
		if ab.id in rep_ids:
			SeqIO.write(ab, handle, "fasta")
			
	

if __name__ == '__main__':
	
	main()

