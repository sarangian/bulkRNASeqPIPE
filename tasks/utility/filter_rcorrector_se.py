#!/usr/bin/env python

"""
author: adam h freedman
afreedman405 at gmail.com
data: Fri Aug 26 10:55:18 EDT 2016
This script takes as an input Rcorrector error corrected Illumina paired-reads
in fastq format and:
1. Removes any reads that Rcorrector indentifes as containing an error,
but can't be corrected, typically low complexity sequences. For these,
the header contains 'unfixable'.
2. Strips the ' cor' from headers of reads that Rcorrector fixed, to avoid
issues created by certain header formats for downstream tools.
3. Write a log with counts of (a) read pairs that were removed because one end
was unfixable, (b) corrected left and right reads, (c) total number of
read pairs containing at least one corrected read.
Currently, this script handles single-end data, and handle either unzipped
or gzipped files on the fly, so long as the gzipped files end with 'gz'.
"""

import sys        
import gzip
#from itertools import izip,izip_longest
import argparse
from os.path import basename

def get_input_streams(file):
    if file[-2:]=='gz':
        handle=gzip.open(file,'rb')

    else:
        handle=open(file,'r')

    return handle

        
def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)  
    

if __name__=="__main__": 
    parser = argparse.ArgumentParser(description="options for filtering and logging rCorrector fastq outputs")
    parser.add_argument('-in','--single_reads',dest='singlereads',type=str,help='Single Reads fastq file')
    parser.add_argument('-o','--out_prefix',dest='outprefix',type=str,help="prefix for filtered fastq output")
    opts = parser.parse_args()

    s_out=open(opts.outprefix+'_'+basename(opts.singlereads).replace('.gz',''),'w')
    s_cor_count=0
    s_unfix_count=0
    s_stream=get_input_streams(opts.singlereads)

    with s_stream as f1:
        S=grouper(f1, 4)
        s_counter=0
        for s_entry in S:
            s_counter+=1
            if s_counter%100000==0:
                print "%s reads processed" % s_counter
        
            head1,seq1,placeholder1,qual1=[i.strip() for i in s_entry]
            if 'unfixable' in head1:
                s_unfix_count+=1
            else:
                if 'cor' in head1:
                    s_cor_count+=1

                head1=head1.split('l:')[0][:-1] # keeps all before the low kmer count statistic and removes the trailing whitespace character
                #head1=head1.replace(' cor','')
                s_out.write('%s\n' % '\n'.join([head1,seq1,placeholder1,qual1]))

    unfix_log=open('rmunfixable.log','w')
    unfix_log.write('total SE reads:%s\nremoved SE reads:%s\n''corrected:%s' % (s_counter,s_unfix_count,s_cor_count))

    s_out.close()
