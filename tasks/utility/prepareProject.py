#!/usr/bin/env python
import optparse
import re
import os
from sys import exit
import pandas as pd
import psutil
from collections import OrderedDict

########################################################################################################################

required="rnaseqDir fileExtention readType domain projectName " \
		 "".split()
parser = optparse.OptionParser()

totalcpus = psutil.cpu_count()
threads = (totalcpus-1)

mem = psutil.virtual_memory()
maxMemory= int((mem.available/1073741824) -1)


parser.add_option('-f', '--fastqDir',
				  help="Path to Directory containing raw RNASeq reads",
				  dest='rnaseqDir'
				  )

parser.add_option('-x', '--fileExtention',
				  help="extensions for a FASTQ files. Allowed extensions are [fq / fq.gz / fastq / fastq.gz]",
				  dest='fileExtention'
				  )


parser.add_option('-r', '--readType',
				  help="RNASeq read type =string[single or paired]",
				  dest='readType'
				  )

parser.add_option('-d', '--domain',
				  help="Organism Domain =string[prokaryote OR eukaryote]",
				  dest='domain'
				  )
parser.add_option('-a', '--adapterFile',
				  help="Path to adapter file",
				  dest='adapterFile'
				  )

parser.add_option('-g', '--genomeDir',
				  help="Path to Directory containing genome.fna and /or genome.gff files ",
				  dest='genomeDir'
				 )

parser.add_option('-p', '--projectName',
				  help="Path to the Project Directory ",
				  dest='projectName',
				  default='MyProject'
				 )

parser.add_option('-s', '--schedulerPort',
				  help="[Optional Parameter] Scheduler Port Number. default =int[8888] ",
				  type="int",
				  default=8082
				 )

parser.add_option('-e', '--emailAddress',
				  help="Provide your email address =email[abc@xyz.com]",
				  dest='emailAddress'
				 )
parser.add_option('-c', '--threads',
				  help="[Optional Parameter, ] Number of threads. Default = (total threads -1)",
				  dest='threads',
				  type="int",
				  default = threads)

parser.add_option('-m', '--maxMemory',
				  help="[Optional Parameter] Maximum Memory in GB. Default = (available memory in GB -1)",
				  dest='maxMemory',
				  type="int",
				  default=maxMemory
				 )

options,args = parser.parse_args()
for r in required:
	if options.__dict__[r] is None:
		parser.error("parameter %s required" %r)


option_dict = vars(options)
rnaseqDir = option_dict.get('rnaseqDir')
suffix = option_dict.get('fileExtention')
domain = option_dict.get('domain')
adapterPath = option_dict.get('adapterFile')
genomeDir = option_dict.get('genomeDir')
genomeName = option_dict.get('organismName')
projectName = option_dict.get('projectName')
readType = option_dict.get('readType')
email=option_dict.get('emailAddress')
port = option_dict.get('schedulerPort')
cpu = option_dict.get('threads')
memory = option_dict.get('maxMemory')

def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)

# Example
createFolder(projectName)

if not (os.path.isdir(rnaseqDir)):
	print ("Path provided for RNASeq Directory does not exists.")
	exit()

if not (os.path.isdir(genomeDir)):
	print ("Path Provided for Genome Directory does not exists.")
	exit()




if os.path.isdir(genomeDir):
	files = [f for f in os.listdir(genomeDir) if os.path.isfile(os.path.join(genomeDir, f))]

	file_list=[]
	for file in files:
		bname = file.split('.')
		genome = bname[0]
		file_list.append(genome)

	genome_set=set(file_list)
	number=len(genome_set)
	if number != 1:
		print("the annotation file (.gff) and genome fasta file (.fna) must have identical prefix and the the "
			  "{genomefolder} must contain only these two files".format(genomefolder=genomeDir))
		exit()

	genome_list=sorted(genome_set)
	genomefile=genome_list[0]

	with open('luigi.cfg', 'w') as fh:
		fh.write('[core]\n')
		fh.write('default-scheduler-port:{port}\n'.format(port=port))
		fh.write('error-email={email}\n\n'.format(email=email))
		fh.write('[GlobalParameter]\n')
		fh.write('rnaSeqDir={rnaseqDir}\n'.format(rnaseqDir=rnaseqDir))
		fh.write('suffix={suffix}\n'.format(suffix=suffix))
		fh.write('domain={domain}\n'.format(domain=domain))
		fh.write('adapter={adapterPath}\n'.format(adapterPath=adapterPath))
		fh.write('genomeDir={genomeDir}\n'.format(genomeDir=genomeDir))
		if domain == "prokaryote":
			fh.write('featureType=CDS\n')
		if domain == "eukaryote":
			fh.write('featureType=exon\n')	
	#fh.write('genomeName={genomeName}\n'.format(genomeName=genomeName))
		fh.write('genomeName={genomeName}\n'.format(genomeName=genomefile))
		fh.write('projectName={projectName}\n'.format(projectName=projectName))
		fh.write('readType={readType}\n'.format(readType=readType))
		fh.write('threads={cpu}\n'.format(cpu=cpu))
		fh.write('maxMemory={memory}\n'.format(memory=memory))
		fh.close()

	'''
	if genome_list[0] != genomeName:
		print("the parameter passed with -o i.e. {genomename} and prefix of genome fasta file i.e. {genomefile} "
			  "present at {genomefolder} must be same".format(genomename=genomeName,genomefile=genomefile ,genomefolder=genomeDir))
		exit()
	'''

if os.path.isdir(rnaseqDir):
	files = [f for f in os.listdir(rnaseqDir) if os.path.isfile(os.path.join(rnaseqDir, f))]
	final_list = []

	if readType=="paired":
		pairedEndList = re.compile(r'^(.+?)_R[12]\.(fastq|fq|fastq\.gz|fq\.gz)?$')

		#pairedEndList = re.compile(r'^(.+?)_R[12]\.(fastq|fq(?:\.gz)?)$')


		for file in files:
			if pairedEndList.search(file):
				if suffix in file:
					bname = file.split('_R')
					sample = bname[0]
					final_list.append(sample)
				else:
					print("\nFile Extension Error\n")
					print("The file extention \"{suffix}\" you provided with -x argument, is not matching with files at {rnaseqDir}\n".format(suffix=suffix,rnaseqDir=rnaseqDir))
					print("List of files at {rnaseqDir}".format(rnaseqDir=rnaseqDir))
					for file in files:
						print (file)
					exit ()

	if readType=="single":
		singlEndList = re.compile(r'^(.+?)\.(fastq|fq|fastq\.gz|fq\.gz)?$')
		for file in files:
			if singlEndList.search(file):
				if suffix in file:
					bname = file.split('.')
					sample = bname[0]
					final_list.append(sample)
				else:
					print("\nFile Extension Error\n")
					print("The file extention \"{suffix}\" you provided with -x argument, is not matching with files at {rnaseqDir}\n".format(suffix=suffix,rnaseqDir=rnaseqDir))
					print("List of files at {rnaseqDir}".format(rnaseqDir=rnaseqDir))
					for file in files:
						print (file)
					exit ()

	sample_set = set(final_list)
	sample_list = sorted(sample_set)
	number=len(sample_list)

	print ("\nfound a total of {number} samples\n".format(number=number))

	for sample in range(len(sample_list)):
		print (sample +1,":" ,sample_list[sample])

	print ('''
	Now you have to enter the biological conditions associated with each sample
	
	Example:
	
	samples         conditions
	sample_s1_rep1  control
	sample_s1_rep2  control
	sample_s2_rep1  treated
	sample_s2_rep2  treated

	Note: (1) condition name is case sensitive.
		  (2) each sample must have atleast one replicate

		   ''')
		
	dictionary=OrderedDict()
	for sample in range(len(sample_list)):
		condition=input('Enter Condition for %s: ' %(sample_list[sample]))
		dictionary.update({sample_list[sample]:condition})
	#print (dictionary)

	df = pd.DataFrame()
	df['lable'] = dictionary.keys()
	df['samples'] = dictionary.keys()
	df['conditions'] = dictionary.values()
	sampledf = df[["samples"]]
	groupdf = df[["samples","conditions"]]


	if (len(df.conditions.unique()) >2) or ((df.groupby('conditions').size().min()) == 1):
		print ('\nError in generating target file\n.................NOTE.................... \n1. Name of the condition is case sensitive. \n2. Each sample must have atleast one replicate.\n3. There must be two different conditions\n\nPlease Check The Values you entered\n')
		
	else :
		print('Target file generated\n')
		df.to_csv(projectName +"/"+"target.tsv", sep='\t', index=False)
		df.to_csv(projectName +"/"+"target.tsv", sep='\t', index=False)
		sampledf.to_csv(projectName + "/" + "samples.txt", index=False, header=False)
		groupdf.to_csv(projectName + "/" + "group.tsv", sep='\t', index=False, header=False)

else:
	print ("Directory not exists.")
