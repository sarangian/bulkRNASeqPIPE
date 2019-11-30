from __future__ import print_function
import pandas as pd
import math
import gzip
#from itertools import izip,izip_longest
import argparse
from os.path import basename
import sys
from Bio import SeqIO


#Function to check sample metadata file
def check_metadata(input_file):
    df = pd.read_table(input_file, sep='\t+', engine='python')

    labels = df.ix[:,0]
    samples = df.ix[:,1]
    groups = df.ix[:,2]
    print('''_______________________________Sample Metadata Sanity Check__________________________''')
    if labels.is_unique == True:
        print ("                    Labels OK")
    else:
        print ('''                  Labels should be unique.
                   Please provide a proper sample metadata file with unique labels. Delete the current project and 
                   prepare a new project with proper metadata file
                  
                   Example
                   -----------------------------------
                   label     samples     group
                   sample1   sample1     condition1
                   sample2   sample2     condition1
                   sample3   sample3     condition2
                   sample4   sample4     condition4
                  
                  ''')

    if samples.is_unique == True:
        print ("                   Sample Names OK")
    else:
        print('''                  SampleNames should be unique.
                    Please provide a proper sample metadata file with unique labels. Delete the current project and 
                    prepare a new project with proper tab separated metadata file.

                    Example
                     -----------------------------------
                    label     samples     group
                    sample1   sample1     condition1
                    sample2   sample2     condition1
                    sample3   sample3     condition2
                    sample4   sample4     condition4     

                    ''')

    if groups.is_unique == True:

        print('''                  SampleNames should be unique.
                    Please provide a tab separated sample metadata file with proper grouping. Delete the current 
                    project and prepare a new project with proper metadata file.

                    Example
                    -----------------------------------
                    label     samples     group
                    sample1   sample1     condition1
                    sample2   sample2     condition1
                    sample3   sample3     condition2
                    sample4   sample4     condition4     

            ''')
    else:
        print ("                   Sample Groups OK.")

    print('''__________________________________________________________________________''')



#Function to calucale the value of genomeSAindexNbases for STAR to index genome
def prepare_sample_list(input_file):
    df = pd.read_table(input_file, sep='\t+', engine='python')
    df1 = df[["samples"]]
    samples = df1.to_csv('samples.txt',header=False,index=False)

    return samples

def prepare_trinity_sample(input_file):
    df = pd.read_table(input_file, sep='\t+', engine='python', header=None)
    df.columns = ['samples','conditions']
    df['count'] = df.groupby('conditions').cumcount() + 1
    df['replicate'] = df["conditions"]+ "_Rep" + df["count"].map(str)
    df['R1'] = df["samples"]+ "_R1.fastq"
    df['R2'] = df["samples"]+ "_R2.fastq"
    df1 = df[['conditions','replicate','R1','R2']]
    trinity_samples = df1.to_csv('trinity_samples.txt',header=False,index=False)
    return trinity_samples

def prepare_corset_target(input_file):
    df = pd.read_table(input_file, sep='\t+', engine='python')
    df1 = df[["samples", "conditions"]]
    df2 = df1.set_index("samples", drop = False)
    df_groups = df2.groupby('conditions')

    condition_group = ",".join(df1['conditions'])
    sample_group = ",".join(df1['samples'])
    CorSet_Input = '-g {} -n {}'.format(condition_group, sample_group)
    return CorSet_Input

'''
def genomeSAindexNbases(genome_fasta):
    with open(genome_fasta) as genome:
        seqs={}
        for line in genome:
            line=line.rstrip()

            if line.startswith('>'):
                words=line.split()
                name=words[0][1:]
                seqs[name]=''
            else:
                seqs[name]=seqs[name]+line
        size=0
        for scaffold in seqs.keys():
            size=size+len(seqs[scaffold])
        #print size
        log2 = math.log(size, 2.0)
    #print log2
        index = (log2/2) - 1
        gsanb= int(index)
        return gsanb'''

def genomeSAindexNbases(genome_fasta):
    with open(genome_fasta) as genome:
        totalsize=0
        for rec in SeqIO.parse(genome, 'fasta'):
            totalsize = totalsize + len(rec)
        log2 = math.log(totalsize, 2.0)
        index = (log2/2) - 1
        gsanb = int(index)
        return gsanb

#Function to prepare rockhopper paired-end input
def prepare_rockhopper_pe_input(input_file):
    df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
    df1.columns = ["SampleName", "Condition"]
    df2 = df1.set_index("SampleName", drop=False)
    df_groups = df2.groupby('Condition')

    grp_list = []
    for name, group in df_groups:
        grp_list.append(name)
    grp_str = ','.join(grp_list)
    parsed_string = "-L "

    group_name_list = []
    for name,group in df_groups:
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'
        left_read_name_list = [x + left_read_name_suffix for x in group.SampleName.tolist()]
        right_read_name_list = [x + right_read_name_suffix for x in group.SampleName.tolist()]

        result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]

        result1 = [x+"%" +y for x, y in result]
        group_string = ','.join(result1)
        group_name_list.append(group_string)

    parse_group_string = " ".join(group_name_list)

    parsed_string = "{} {} {}".format(parsed_string, grp_str, parse_group_string)
    return parsed_string


def prepare_rockhopper_se_input(input_file):

    df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
    df1.columns = ["SampleName", "Condition"]
    df2 = df1.set_index("SampleName", drop=False)
    df_groups = df2.groupby('Condition')

    grp_list = []
    for name, group in df_groups:
        grp_list.append(name)
    grp_str = ','.join(grp_list)
    parsed_string = "-L "

    group_name_list = []
    for name,group in df_groups:
        read_name_suffix = '.fastq'
        read_name_list = [x + read_name_suffix for x in group.SampleName.tolist()]

        result1 = [x for x in read_name_list]
        group_string = ','.join(result1)
        group_name_list.append(group_string)

    parse_group_string = " ".join(group_name_list)

    parsed_string = "{} {} {}".format(parsed_string, grp_str, parse_group_string)
    return parsed_string

#Function to prepare rna-spades paired-end input
def prepare_rnaspades_pe_input(input_file):
    with open(input_file) as ifh:
        sample_name_list = ifh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        left_read_name_prefix = '--pe1-1 '
        right_read_name_suffix = '_R2.fastq'
        right_read_name_prefix = '--pe1-2 '

        left_read_name_list = [x + left_read_name_suffix for x in sample_name_list]
        left_read_name_list = [left_read_name_prefix + x for x in left_read_name_list]

        right_read_name_list = [x + right_read_name_suffix for x in sample_name_list]
        right_read_name_list = [right_read_name_prefix + x for x in right_read_name_list]

        result = [item for sublist in zip(left_read_name_list, right_read_name_list) for item in sublist]
        RNASpades_PE_Input = " ".join(result)
        return RNASpades_PE_Input

#Function to prepare rna-spades single-end input
def prepare_rnaspades_se_input(input_file):
    with open(input_file) as ifh:
        sample_name_list = ifh.read().splitlines()
        read_name_suffix = '.fastq'
        read_name_prefix = '--s1 '

        read_name_list = [x + read_name_suffix for x in sample_name_list]
        read_name_list = [read_name_prefix + x for x in read_name_list]

        RNASpades_SE_Input = " ".join(read_name_list)
        return RNASpades_SE_Input


#Function to prepare Trinity paired-end input
def prepare_trinity_pe_input(input_file):
    with open(input_file) as ifh:
        sample_name_list = ifh.read().splitlines()
        left_read_name_suffix = '_R1.fastq'
        right_read_name_suffix = '_R2.fastq'

        left_read_name_list = [ x + left_read_name_suffix for x in sample_name_list]
        right_read_name_list =[ x + right_read_name_suffix for x in sample_name_list]

        left_reads = ','.join(left_read_name_list)
        right_reads = ','.join(right_read_name_list)

        Trinity_PE_Input = "--left " + left_reads + " --right " + right_reads
        return Trinity_PE_Input

#Function to prepare Trinity Single-end input
def prepare_trinity_se_input(input_file):
    with open(input_file) as ifh:
        sample_name_list = ifh.read().splitlines()
        read_name_suffix = '.fastq'
        read_name_list = [ x + read_name_suffix for x in sample_name_list]
        reads = ','.join(read_name_list)
        Trinity_SE_Input = "--single " + reads + " "
        return Trinity_SE_Input


#Function to prepare Corset Input
def prepare_corset_input(input_file):
    df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
    df1.columns = ["SampleName", "Condition"]
    df2 = df1.set_index("SampleName", drop = False)
    df_groups = df2.groupby('Condition')

    condition_group = ",".join(df1['Condition'])
    sample_group = ",".join(df1['SampleName'])
    CorSet_Input = '-g {} -n {}'.format(condition_group, sample_group)
    return CorSet_Input


def format_rockhopper_transcripts(input_file,length):
    data = pd.read_table(input_file, sep='\t+', engine='python')
    transcripts = data[data.Length >= int(length)]

    df = transcripts.assign(Transcript_ID=[1 + i for i in xrange(len(transcripts))])[
        ['Transcript_ID'] + transcripts.columns.tolist()]

    fasta_file = []
    for index, row in df.iterrows():
        header = ">transcript_" + str(row['Transcript_ID']) + "_length_" + str(row['Length'])
        nucleotides = row['Sequence']
        fasta = "{}\n{}\n".format(header,nucleotides)
        fasta_file.append(fasta)
        out_fasta = "".join(fasta_file)
    return out_fasta

