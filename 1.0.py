#!/usr/bin/env python
# encoding: utf-8
"""
1.0.py -p project

Created by Mingchen on 2015-05-04.
Copyright (c) 2015 __MyCompanyName__. All rights reserved

"""
import glob, os, sys, subprocess
from Bio import SeqIO
from mytools_ymc import *

def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch

def trim_fastq_by_quality(the_file,prj_folder):
	handle = open(the_file, "rU")
	fname = retrieve_name_body(the_file)
	print "Triming...",fname
	writer = open("%s/1.1-trimed-fastq-file/%s-trimed.fastq"%(prj_folder,fname), "w")
	for record in SeqIO.parse(handle, "fastq") :
		quality_type = list(record.letter_annotations)[0]
		quality_list = record.letter_annotations[quality_type]
		position_list = []
		for index in range(0,len(quality_list)):
			if quality_list[index] > 20:
				position_list.append(index)
		new_record = record[position_list[0] : position_list[-1]+1]
		SeqIO.write(new_record, writer, "fastq")
	handle.close()
	writer.close()
def unique_fasta(prj_folder, project):
	handle = "%s/1.2-merged-fastq-file/%s.assembled.fasta"%(prj_folder, project)
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
	print "Quality contorl..."
	infiles = glob.glob("%s/1.0-origin/*.fastq"%(prj_folder))
	if len(infiles) != 2:
		print "The %s be loaded in error, are they two?"%infiles
	for the_file in infiles:
		trim_fastq_by_quality(the_file,prj_folder)	
	print "Merging..."	
	infiles = glob.glob("%s/1.1-trimed-fastq-file/*.fastq"%(prj_folder))
	
	os.chdir("%s/1.2-merged-fastq-file/"%(prj_folder))
	merge = subprocess.call("pear -f %s -r %s -o %s"%(infiles[0],infiles[1],project),shell=True)
	
	print "Convert fastq to fasta..."
	merged_file = "%s/1.2-merged-fastq-file/%s.assembled.fastq"%(prj_folder, project)
	fname, suffix = os.path.splitext(merged_file)
	count = SeqIO.convert(merged_file, "fastq","%s.fasta"%fname, "fasta")
	print count
	print "There are  %i records have been Converted!" %(count)
	
	print "Unique fasta file..."
	unique_fasta(prj_folder, project)
	
	print "Split large file to small..."
	os.chdir("%s/1.3-splited-fasta-file/"%(prj_folder))
	record_iter = SeqIO.parse(open("%s_unique.fasta"%fname), "fasta")
	for i, batch in enumerate(batch_iterator(record_iter, 10000)) :
		filename = "%s-%i.fasta" % (project, i+1)
		handle = open(filename, "w")
		count = SeqIO.write(batch, handle, "fasta")
		handle.close()
		print "Wrote %i records to %s" % (count, filename)
	print "Begin IgBLAST..."
	os.chdir("%s/1.4-IgBLAST-output/"%(prj_folder))
	IgBLAST_result = open("%s-igblast-output.txt "%(project),"w")	
	mv_database = subprocess.call("cp -r /zzh_gpfs/home/zzhgroup/yanmingchen/IgBLAST_database/ ./",shell = True)
	#IgBLAST_run = subprocess.call("igblastn -germline_db_V ./IgBLAST_database/20150429-human-gl-v -germline_db_J ./IgBLAST_database/20150429-human-gl-j -germline_db_D ./IgBLAST_database/20150429-human-gl-d -organism human -domain_system imgt -query %s -auxiliary_data optional_file/human_gl.aux -outfmt '7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq score frames qframe sframe positive ppos btop staxids stitle sstrand qcovs qcovhsp' -num_alignments_V 10 -num_alignments_D 10 -num_alignments_J 10 -out IgBLAST_result"%merged_file,shell=True)
	IGBLAST_infiles = glob.glob("%s/1.3-splited-fasta-file/*.fasta"%(prj_folder))
	for index, the_file in enumerate(IGBLAST_infiles):
		print index,the_file
		IgBLAST_run = subprocess.call("igblastn -germline_db_V ./IgBLAST_database/20150429-human-gl-v -germline_db_J ./IgBLAST_database/20150429-human-gl-j -germline_db_D ./IgBLAST_database/20150429-human-gl-d -organism human -domain_system imgt -query %s -auxiliary_data optional_file/human_gl.aux -outfmt '7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq score frames qframe sframe positive ppos btop staxids stitle sstrand qcovs qcovhsp' -num_alignments_V 10 -num_alignments_D 10 -num_alignments_J 10 -out IgBLAST_result_%i &"%(the_file,index+1),shell=True)
		#IGBLAST_runner = subprocess.call("igblastn -germline_db_V /zzh_gpfs/home/zzhgroup/yanmingchen/IgBLAST_database/20150429-human_gl_v -germline_db_J /zzh_gpfs/home/zzhgroup/yanmingchen/IgBLAST_database/20150429-human_gl_j -germline_db_D /zzh_gpfs/home/zzhgroup/yanmingchen/IgBLAST_database/20150429-human_gl_d -organism human -domain_system imgt -query %s -auxiliary_data optional_file/human_gl.aux -outfmt '7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qlen slen qseq sseq score frames qframe sframe positive ppos btop staxids stitle sstrand qcovs qcovhsp' -num_alignments_V 10 -num_alignments_D 10 -num_alignments_J 10 -out ./%s-%i.txt & "%(the_file,project,index+1),shell=True)
		#fname, suffix = os.path.splitext(the_file)
		#the_count = SeqIO.convert(the_file, 'fastq', '%s.fasta'%fname, "fasta")
		#print "%s Converted %i records" % (fname,the_count)
#main()


if __name__ == '__main__':

	# get parameters from input
	if len(sys.argv) < 3 :
		print __doc__
		sys.exit(0)

	dict_args = processParas(sys.argv, p = "project")
	project = getParas(dict_args, "project")
	
	# create 1st and 2nd subfolders
	create_folders(os.getcwd())
	main()
	print "Finished"
