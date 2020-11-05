#!/usr/bin/env python
import pandas as pd
import optparse
import sys

parser = optparse.OptionParser()

parser.add_option('-i', '--infile',
                  help="input file name (.tsv)",
                  default=False,
                  )

parser.add_option('-l', '--minlen',
                  help="input file name (.tsv)",
                  default=False,
                 )

#parser.add_option('-o', '--outfile',
#                  help="output file name",
#                  default=False,
#                  )


options,args = parser.parse_args()
option_dict = vars(options)

inf = option_dict.get('infile')
leng = option_dict.get('minlen')

#print 'input file  :', inf
#print 'output file :',  outf
#print 'minimum transcript length: ', leng
data = pd.read_table(inf, sep='\t+', engine='python')
transcripts = data[data.Length >= int(leng)]

df = transcripts.assign(Transcript_ID=[1 + i for i in xrange(len(transcripts))])[['Transcript_ID'] + transcripts.columns.tolist()]

for index, row in df.iterrows():
    print ">transcript_" + str(row['Transcript_ID']) + "_length_"  + str(row['Length'])
    print row['Sequence']


