#!/usr/bin/env python
"""
Deblur sequences
Created on Thu Sep 19 17:43:09 2013
@author: amnon
"""

__version__ = "1.995"

import argparse

from os import listdir
from os.path import isfile,join
import sys
import numpy as np
from cogent.parse.fasta import MinimalFastaParser
# for the pairwise alignment
from cogent.align.align import global_pairwise
from cogent.align.align import make_dna_scoring_dict
from cogent import DNA
from cogent.align.algorithm import nw_align
#from Bio import pairwise2

import LogMe

def SeqToArray(seq):
    """ convert a string sequence to a numpy array"""
    seqa=np.zeros(len(seq),dtype=np.int8)
    for ind,base in enumerate(seq):
        if base=='A':
            seqa[ind]=0
        elif base=='a':
            seqa[ind]=0
        elif base=='C':
            seqa[ind]=1
        elif base=='c':
            seqa[ind]=1
        elif base=='G':
            seqa[ind]=2
        elif base=='g':
            seqa[ind]=2
        elif base=='T':
            seqa[ind]=3
        elif base=='t':
            seqa[ind]=3
        elif base=='-':
            seqa[ind]=4
        else:
            seqa[ind]=5
    return(seqa)


def RemoveError(log,seqs,seqsnp,sfreq,readerror,meanerror,ofracerr,indelprob,indelmax,pyroseq):
    """ Deblur the reads
    Input:
        log - a LogMe log module to write the debluring info
        seqs - the list of sequences
        seqsnp - a list of numpy arrays of the sequences (for faster comparison) - from SeqToArray()
        sfreq - dictionary (based on the sequence) of the number of reads for each sequence
        readerror - the maximal read error expected (fraction - typically 0.01)
        meanerror - the mean read error used for peak spread normalization - typically 0.01
        ofracerr - the error distribution array, or 0 if use default
        indelprob - the probability for an indel (currently constant for number of indels until max is reached)
        indelmax - the maximal number of indels expected by errors (error cutoff)
        pyroseq - if set, use pairwise alignment for pyrosequencing data
    Output:
        sfreq - the deblurred number of reads for each sequence (0 if not present)
        debugdata - a list of strings
    Notes:
        meanerror is used only for normalizing the peak height before deblurring, whereas readerror
        is used for calculating the expected number of errors for each position
        error distribution array X should be of length >10, where Xi = max frequency of error hamming i
        if it is 0, we use the default distribution
    """
    # take the list values so it won't change
    fracerr=list(ofracerr)
    
    # we assume all sequences are of equal length
    commonlen=len(seqs[0])
    for cseq in seqs:
        if not(commonlen==len(cseq)):
            print("Not all sequences are same length!!!!")
            print(commonlen)
            print(len(cseq))
            print(cseq)
    print ("processing",len(seqs),"sequences")

    numreal=0
    for cchar in seqs[0]:
        if not (cchar=='-'):
            numreal+=1
    modfactor=pow((1-meanerror),numreal)

    # create the error profile from the read error
    # exponential independent
    #   fracerr=[]
    #   for a in range(10):
    #       fracerr.append(pow(readerror,a)/modfactor)

    # empirical
    #    fracerr=[1.0/modfactor,pow(readerror,1)/modfactor,2*pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor,pow(readerror,2)/modfactor]

    # used for the 22 mock mixture
    #    fracerr=[1.0/modfactor,pow(readerror,1)/modfactor,0.01,0.01,0.01,0.005,0.005,0.005,0.005,0.005,0.005,0.001,0.001,0.001,0.001,0.001,0.001,0.0005,0.0001,0.0001]

    # used for the 44 mock mixture
    #   e1=pow(readerror,1)/modfactor
    #   fracerr=[1.0/modfactor,e1,e1/4,e1/5,e1/6,e1/8,e1/10,e1/15,e1/20,e1/30,e1/40,e1/50,e1/50,e1/50,e1/50,e1/50,e1/50,e1/100,e1/500,e1/500]

    # if fracerr not supplied, use the default (22 mock mixture setup)
    log.log("original fracer parameter:",fracerr)
    if fracerr==0:
        fracerr=[1.0/modfactor,pow(readerror,1)/modfactor,0.01,0.01,0.01,0.005,0.005,0.005,0.005,0.005,0.005,0.001,0.001,0.001,0.001]
        log.log("modified fracerr because it was 0")
    else:
        for idx,val in enumerate(fracerr):
            fracerr[idx]=fracerr[idx]/modfactor

    maxhdist=len(fracerr)-1

    print "fracerr"
    print fracerr
    print "readerror"
    print readerror
    print "modfactor"
    print modfactor   

    log.log("indel prob:",indelprob)
    log.log("indel max:",indelmax)
    log.log("readerror:",readerror)
    log.log("meanerror:",meanerror)
    log.log("mod factor:",modfactor)
    log.log("fracerr:",fracerr)     

    # for pairwise alignment:
    DNAm = make_dna_scoring_dict(10, -8, -8)

    for idx,cseq in enumerate(seqs):
        csfreq=sfreq[cseq]
        # no need to remove neighbors if freq. is <=0
        if csfreq<=0:
            continue
        # correct for the fact that many reads are expected to be mutated
        numerr=[] 
        for a in range(len(fracerr)):
            numerr.append(fracerr[a]*csfreq)

        # if it's low level, just continue
        if numerr[1]<0.1:
            continue

        # compare to all other sequences and calculate hamming dist
        cseqnp=seqsnp[idx]
        oseqlen=len(seqs[idx].rstrip('-'))
        for idxtmp,seqnptmp in enumerate(seqsnp):
            # don't compare to ourselves (dist=0)
            if idxtmp==idx:
                continue
            # calculate the hamming distance
            hdist=np.count_nonzero(np.not_equal(seqnptmp,cseqnp))
            # if far away, don't need to correct
            if hdist>maxhdist:
                continue
            # close, so lets calculate exact distance

            numsub=0
            numindel=0
            # experimental try 2
            # s1=seqs[idx].replace('-','')
            # s2=seqs[idxtmp].replace('-','')
            # cseq1,cseq2=nw_align(s1,s2)

            # experimental: pairwise align the sequences
            if pyroseq:
                s0=DNA.makeSequence(seqs[idx])
                s0=s0.degap()
                s1=DNA.makeSequence(seqs[idxtmp])
                s1=s1.degap()
                print s0._seq
                print s1._seq
                align = global_pairwise(s0, s1, DNAm, 10, 9)
                a0=align.getGappedSeq('seq_0')
                a1=align.getGappedSeq('seq_1')            
                cseq1=a0._seq
                cseq2=a1._seq

                len1=len(cseq1.rstrip('-'))
                len2=len(cseq2.rstrip('-'))
                oseqlen=len(cseq1)
                for cpos in range(oseqlen):
                    if not (cseq1[cpos]==cseq2[cpos]):
                        if cseq1[cpos]=='-':
                            if cpos<len1:
                                numindel+=1
                        else:
                            if cseq2[cpos]=='-':
                                if cpos<len2:
                                    numindel+=1
                            else:
                                numsub+=1

            # not pyrosequencing so use the faster global alignment
            else:
                for cpos in range(oseqlen):
                    if not (cseqnp[cpos]==seqnptmp[cpos]):
                        # 4 is '-'
                        if seqnptmp[cpos]==4:
                            numindel+=1
                        else:
                            if cseqnp[cpos]==4:
                                numindel+=1
                            else:
                                numsub+=1

            nerr=numerr[numsub]

            # remove errors due to (PCR?) indels (saw in 22 mock mixture)
            if numindel>0:
                nerr=nerr*indelprob
            if numindel>indelmax:
                nerr=0

            # if the effect is small - don't do anything
            if nerr<0.1:
                continue
            # met all the criteria - so correct the frequency of the neighbor
            sfreq[seqs[idxtmp]]-=nerr
            # if sfreq[seqs[idxtmp]]<=0:
            #     if sfreq[seqs[idxtmp]]+nerr>0:
            #         log.log("Removed sequence ",idxtmp," due to sequence ",idx)
            #         log.log("seq:",idx," and ",idxtmp," have ",numindel," indels and ",numsub,"substitutions")
            #         log.log(cseq1)
            #         log.log(cseq2)
            #         log.log("true seq freq:",csfreq)
            #         log.log("freq from ",sfreq[seqs[idxtmp]]+nerr," to ",sfreq[seqs[idxtmp]])
            # else:
            #     if numindel>0:
            #         log.log("====indels but no delete!!!!")
            #         log.log("seq:",idx," and ",idxtmp," have ",numindel," indels and ",numsub,"substitutions")
            #         log.log(cseq1)
            #         log.log(cseq2)
            #         log.log("true seq freq:",csfreq)
            #         log.log("freq from ",sfreq[seqs[idxtmp]]+nerr," to ",sfreq[seqs[idxtmp]])
    return(sfreq)




def CleanSeqs(dirname,readerror,meanerror,errordist,indelprob,indelmax,pyroseq):
    print "Cleaning"
# prepare the file list in the directory
    filelist=[f for f in listdir(dirname) if isfile(join(dirname,f))]

# loop over the files
    for cfile in filelist:
        if cfile[-5:]!='.tuni':
            continue
        sfreq={}
        seqs=[]
        seqnames=[]
        seqsnp=[]
        sortfreq=[]
        log=LogMe.LogMe(join(dirname,cfile)+'.log')
        log.log("CleanSeqsIndel.py version",__version__)
        log.log("Original dir name:",dirname)
        log.log("Original file name:",cfile)
        print cfile
        fafile=open(join(dirname,cfile))
        for seqid,seq in MinimalFastaParser(fafile):
            # get the number or reads from the header string
            numseqs=float(seqid[seqid.find(';size=')+6:-1])
            
            # convert sequence to a numpy array
            seqa=SeqToArray(seq)
            # hash the number of reads
            sfreq[seq]=numseqs
            # and store the list of sequences (needs to be sorted)
            # we use a hash for the frequencies and a numpy array for the hamming comparisons
            seqs.append(seq)
            seqnames.append(seqid)
            seqsnp.append(seqa)
            # store the ordered list of frequencies for sorting
            sortfreq.append(numseqs)
        
        if len(seqs)==0:
            print("No sequences in file: "+cfile)
        else:
            # sort the sequences according to frequencies (descending)
            sortorder=sorted(range(len(sortfreq)),key=sortfreq.__getitem__, reverse=True)
            seqs=[seqs[i] for i in sortorder]
            seqnames=[seqnames[i] for i in sortorder]
            seqsnp=[seqsnp[i] for i in sortorder]

 #           for cseq in seqs:
 #               log.log(cseq," size=",sfreq[cseq])
            # after loading the file - remove the read errors (MAIN FUNCTION)
            print "orig readerror",readerror
            cfreq=RemoveError(log,seqs,seqsnp,sfreq,readerror,meanerror,errordist,indelprob,indelmax,pyroseq)
#            log.log("=================================")
#            for cseq in seqs:
#                log.log(cseq," size=",sfreq[cseq])
            # and finally save the new fasta as a '.clean' file
            ofile=open(join(dirname,cfile+'.clean'),'w')
            cid=1
            for cseq in seqs:
                if round(cfreq[cseq])>0:
                    ofile.write('>aa'+str(cid)+';size='+str(int(round(cfreq[cseq])))+';\n')
                    cid+=1
                    csequp=cseq.upper()
                    for a in range(len(csequp)):
                        if not (csequp[a]=='-'):
                            ofile.write(csequp[a])
                    ofile.write('\n')
            ofile.close()
    

def main(argv):
    parser=argparse.ArgumentParser(description='Clean read errors from illumina reads. version '+__version__)
    parser.add_argument('dirname',help='input dir (containing .tuni files unique and truncated)')
    parser.add_argument('-e','--readerror',help='readerror rate',default=0.05)
    parser.add_argument('--indelmax',help='maximal indel number',default=3)
    parser.add_argument('-i','--indelprob',help='indel probability (same for N indels)',default=0.01)
    parser.add_argument('-p','--pyroseq',help='Use pairwise alignment for pyrosequencing (slower)',action='store_true')
    parser.add_argument('-m','--meanerror',help='the mean error, used for original sequence estimate (default same as readerror)',default=-1)
    parser.add_argument('-d','--errordist',help='a comma separated list of error probabilities for each edit distance (min length=10)',default=0)
    args=parser.parse_args(argv)
    if args.meanerror==-1:
        args.meanerror=args.readerror

    # cast the error profile to a list if not 0 (default error profile)
    errordist=args.errordist
    if not errordist==0:
        errorstr=errordist
        errordist=errorstr.split(',')
        errordist=list(map(float,errordist))
    
    CleanSeqs(args.dirname,float(args.readerror),float(args.meanerror),errordist,float(args.indelprob),float(args.indelmax),args.pyroseq)
    
if __name__ == "__main__":
    main(sys.argv[1:])                