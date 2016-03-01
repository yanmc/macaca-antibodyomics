#!/usr/bin/python
import csv
from Bio import SeqIO
#f = open('human.rna.gbff','r')
f = open('/zzh_gpfs/data/human.rna.gbff','r')
outfile = csv.writer(open('./reference_annotation_ncbi.txt','w'),delimiter = '\t')
outfile.writerow(['NM_**','Definition','Cytogenetic','Gene_length','geneName(ID)','Num of exon','Comment','Reference'])

count_exon = 0
for record in SeqIO.parse(f, "genbank"):
	rec_id = record.id
	rec_difi = record.description
	comm = record.annotations['comment']
	summary_start = comm.find('Summary:')
	if summary_start != -1:
		summary_end = comm.find(']')
		rec_comment = comm[summary_start:summary_end+1]
	else:
		rec_comment = '-'
	rec_reference = record.annotations['references']
	for rec in record.features:
		if rec.type == 'gene':
			gene_len = len(rec.location)
			gene_name = (''.join(rec.qualifiers['gene'])+'['+''.join(rec.qualifiers['db_xref'][0])+']')
			gene_name_1 = ''.join(rec.qualifiers['gene'])
		if rec.type == 'source':
			cytog = (''.join(rec.qualifiers['map']))
		if rec.type == 'exon':
			count_exon += 1
		#if rec.type == 'variation':
			#var_loc = int(rec.location.end)
		#if rec.type == 'CDS':
			#if rec.location.strand == -1:
				#cds_seq = record.seq.reverse_complement()[rec.location.start-1:rec.location.end]
			#else:
				#cds_seq = record.seq[rec.location.start-1:rec.location.end]
		
	
	outfile.writerow([gene_name_1,rec_id, rec_difi, cytog, gene_len, gene_name, count_exon,rec_comment,rec_reference])
	