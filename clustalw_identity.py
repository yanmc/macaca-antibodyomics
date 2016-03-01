#!/usr/bin/python
import Bio
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline
from mytools import *
import sys,os,glob
from Bio.Align.Applications import ClustalwCommandline

def caculate_identity(clustalw_result, result):
	
	fs = glob.glob(clustalw_result)
	for infile in fs:
		print "processing %s"%infile
		iden = 0
		align = AlignIO.read(infile, 'clustal')
		summary_align = AlignInfo.SummaryInfo(c_align)
		dict_align = zip(align[0].seq, align[1].seq)
		for index, (i,j) in enumerate(dict_align):
			if i==j and  i!='-':
				iden+=1
		identity = round(float(iden)/float(len(align[0].seq))*100,2)
		result.writerow([align[0].id, align[1].id, identity])

def do_clustalw(file_for_clustalw):
	infiles = glob.glob(file_for_clustalw)

	#clustalw_exe = r"/Applications/clustalw-2.1-macosx/clustalw2"
	clustalw_exe = r"/zzh_gpfs/apps/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"
	assert os.path.isfile(clustalw_exe), "Clustal W executable missing"
	for in_file in infiles:
		print "Processing %s ......."%in_file
		in_file = in_file.replace('&','\&')
		in_file = in_file.replace('*','\*')
		clustalw_cline = ClustalwCommandline(clustalw_exe, infile=in_file)
		stdout, stderr = clustalw_cline()



def main():
	result = csv.writer(open('%s.txt'%timepoint,'a+'),delimiter = '\t')
	row_germ = SeqIO.index(infile_a,'fasta')
	rank_germ = SeqIO.index(infile_b,'fasta')
	for index1, i in enumerate(row_germ):
		for index2,j in enumerate(rank_germ):
			out = open('%s_%s_pair.fasta'%(index1+1,index2+1),'w')
			SeqIO.write(row_germ[i],out,'fasta')
			SeqIO.write(rank_germ[j],out,'fasta')
			out.close()
			file_for_clustalw = '%s_%s_pair.fasta'%(index1+1,index2+1)
			do_clustalw(file_for_clustalw)
			clustalw_result = '%s_%s_pair.aln'%(index1+1,index2+1)
			caculate_identity(clustalw_result, result)
			os.system(rm "%s_%s_pair.*"%(index1+1,index2+1))
if __name__ == '__main__':
	# check parameters
	if len(sys.argv) < 7 :
		print __doc__
		sys.exit(0)

		# get parameters from input
		dict_args = processParas(sys.argv, a="infile_a", b="infile_b", t="timepoint")
		infile_a, infile_b, timepoint = getParas(dict_args, "infile_a", "infile_b", "timepoint")

		main()