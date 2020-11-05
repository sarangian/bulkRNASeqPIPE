#!/usr/bin/env python3
import luigi
import sys



from tasks.readCleaning.rawReadQC import readqc
from tasks.readCleaning.rawReadQC import rawReadsQC
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.preProcessReads import cleanReads
from tasks.readCleaning.reFormatReads import reformat


#ALIGNMENT FREE
from tasks.rnaSeq.alignmentFree import indexTransctriptome
from tasks.rnaSeq.alignmentFree import generate_transcript_count_file
from tasks.rnaSeq.alignmentFree import alignment_free_differential_expression_analysis


#ALIGNMENT BASED 
from tasks.rnaSeq.alignmentBased import alignment_based_differential_expression_analysis
from tasks.rnaSeq.alignmentBased import align_rnaseq_reads_with_genome
from tasks.rnaSeq.alignmentBased import generate_gene_count_file
from tasks.rnaSeq.alignmentBased import index_genome

#from tasks.rnaSeq.deNovo import transcriptome_assembly
#from tasks.rnaSeq.deNovo import index_quanify_cluster_da_transcripts
from tasks.rnaSeq.deNovo import denovo_differential_expression_analysis


from tasks.rnaSeq.annotation import prokaryotic_annotation
from tasks.rnaSeq.annotation import make_tx_to_gene



if __name__ == '__main__':

    luigi.run()
