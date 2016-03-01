#!/usr/bin/env python


from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
def get_dict_id_seq_exchange(f):
	rec_seq, file_dict = [], {}
	handle = open(f,"rU")
	for record in SeqIO.parse(handle, "fasta"):
		rec_seq.append(record.seq)
		file_dict[record.seq] = record.id
	handle.close()
	return rec_seq, file_dict
	

def dict_merge(dict1, dict2):#For same dict key,It will choose dict2's value. 
	dictMerged2=dict(dict1, **dict2)
	return dictMerged2
def main():
	merged_file = open("macaque-merged-trimmed-Variable-region-min-sequence-unique.fasta", "w")
	infile_411 = "macaca-411_get_trimmed_Variable_region-min-sequence_unique.fasta"
	infile_426 = "macaca-426_get_trimmed_Variable_region-min-sequence_unique.fasta"
	infile_55 = "macaca-55_get_trimmed_Variable_region-min-sequence_unique.fasta"
	infile_520 = "macaca-520_get_trimmed_Variable_region-min-sequence_unique.fasta"
	infile_625 = "macaca-625_get_trimmed_Variable_region-min-sequence_unique.fasta"
	
	print "Step: 1"
	id_411, dict_411 = get_dict_id_seq_exchange(infile_411)
	id_426, dict_426 = get_dict_id_seq_exchange(infile_426)
	id_55, dict_55 = get_dict_id_seq_exchange(infile_55)
	id_520, dict_520 = get_dict_id_seq_exchange(infile_520)
	id_625, dict_625 = get_dict_id_seq_exchange(infile_625)
	
	print "merge dict"
	dict_merge1 = dict_merge(dict_426, dict_411)
	dict_merge2 = dict_merge(dict_55, dict_merge1)
	dict_merge3 = dict_merge(dict_520, dict_merge2)
	dict_merge4 = dict_merge(dict_625, dict_merge3)
	
	#print "set new id"
	#for index, (key, value) in enumerate(dict_merge4.items()):
	#	dict_merge4[key] = ">"+"macaca"+"_"+"%s"%value.split('-')[-1]+"_"+"%s"%index
	print "Write file"
	for key, value in dict_merge4.items():
		record = SeqRecord(key, id = value, description = "")
		SeqIO.write(record, merged_file, "fasta")

main()