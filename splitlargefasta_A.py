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

import os, sys, glob
from Bio import SeqIO
from mytools import *
def spliter(infile):
	fname, suffix = os.path.splitext(infile)
	record_iter = SeqIO.parse(open(infile),"fasta")
	for i, batch in enumerate(batch_iterator(record_iter, 10000)) :
		filename = "%s-%s-split-%s.fasta" % (fname,chain,(i+1))
	   	handle = open(filename, "w")
	   	count = SeqIO.write(batch, handle, "fasta")
 	  	handle.close()
 	  	print "Wrote %i records to %s" % (count, filename)
def main():
	infiles = glob.glob(infile_model)
	for infile in infiles:
		
		spliter(infile)
if __name__ == '__main__':
        # check parameters
        if len(sys.argv) < 7 :
                print __doc__
                sys.exit(0)

        # get parameters from input
        dict_args = processParas(sys.argv, org = "organism",c = "chain",i="infile")
        infile, chain,organism= getParas(dict_args,"infile","chain","organism")

        main()
