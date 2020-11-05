#! /usr/bin/env python3
import os
import re
import shutil
from collections import OrderedDict
import optparse
from sys import exit
import psutil
import subprocess
import pandas as pd

totalcpus = psutil.cpu_count()
threads = (totalcpus-1)
mem = psutil.virtual_memory()
maxMemory= int((mem.available/1073741824) -1)


required="dataDir readType organismDomain projectDir".split()

parser = optparse.OptionParser()
parser.add_option('-i', '--inputDir',
				  help="[OPTIONAL] Path to Input Data Directory, Default 'raw_data'",
				  dest='dataDir',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"raw_data")))
				  )

parser.add_option('-p', '--projectDir',
				  help="Name of the project directory, Default 'alignment_free_dea'",
				  dest='projectDir',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"alignment_free_dea")))
				  )

parser.add_option('-r', '--readType',
				  help="[OPTIONAL] rnaSeq read type. Choose for 'se' or 'pe' ",
				  dest='readType',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"raw_data")))
				  )


parser.add_option('-d', '--organismDomain',
				  help="[OPTIONAL] rnaSeq read type. Choose for 'prokaryote' or 'eukaryote' ",
				  dest='organismDomain'
				  )


parser.add_option('-o', '--symLinkDir',
				  help="[Optional] Path to symbolic link Data Directory, Default 'raw_data_symlink'",
				  dest='symLinkDir',
				  default=(os.path.abspath(os.path.join((os.getcwd()),"raw_data_symlink")))
				  )


parser.add_option('-s', '--schedulerPort',
				  help="[Optional Parameter] Scheduler Port Number. default =int[8888] ",
				  type="int",
				  default=8082
				 )

parser.add_option('-e', '--emailAddress',
				  help="Provide your email address =email[abc@xyz.com]",
				  dest='emailAddress',
				  default="a.n.sarangi@gmail.com"
				 )

parser.add_option('-t', '--threads',
				  help="[Optional Parameter, ] Number of threads. Default = (total threads -1)",
				  dest='threads',
				  type="int",
				  default = threads)

parser.add_option('-x', '--maxMemory',
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

dataDir = option_dict.get('dataDir')
outDir = option_dict.get('symLinkDir')
projectDir=option_dict.get('projectDir')

readType=option_dict.get('readType')
domain=option_dict.get('organismDomain')

adapter=os.path.join(os.getcwd(),"tasks","utility",'adapters.fasta.gz')



def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)


def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output


excludeDir=os.path.abspath(os.path.join(outDir, "ex" ))
pairedendDir=os.path.abspath(os.path.join(outDir, "pe" ))
singleendDir=os.path.abspath(os.path.join(outDir, "se"))
genomeDir=os.path.abspath(os.path.join(outDir, "geno" ))
transcriptomeDir=os.path.abspath(os.path.join(outDir, "tran"))
annotationDir=os.path.abspath(os.path.join(outDir, "anot"))
projectDir=os.path.abspath(os.path.join(projectDir))




email=option_dict.get('emailAddress')
port = int(option_dict.get('schedulerPort'))
cpu = int(option_dict.get('threads'))
memory = int(option_dict.get('maxMemory'))

createFolder(excludeDir)
createFolder(outDir)
createFolder(pairedendDir)
createFolder(singleendDir)
createFolder(transcriptomeDir)
createFolder(genomeDir)
createFolder(annotationDir)
createFolder(projectDir)
createFolder("config")

if os.path.isdir(dataDir):
	files = [f for f in os.listdir(dataDir) if os.path.isfile(os.path.join(dataDir, f))]
	keys = []
	fileList = re.compile(r'^(.+?).(gff|gtf|fasta|fna|ffn|fa|fastq|fq|fastq\.gz|fq\.gz)?$')
	for file in files:
		if fileList.search(file):
			keys.append(file)

dicts = OrderedDict ()
#keys = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

for i in keys:
		dicts[i] = str(input("Enter Data Type of {data}: \tchoose from [pe, se, geno, tran, anot, ex]: ".format(data=i)))
print(dicts)

for key, val in dicts.items():
	if not os.path.exists(os.path.join(outDir, val)):
		os.mkdir(os.path.join(outDir, val))

##ln -nsf method
for key, val in dicts.items():
	dest = (os.path.join(outDir,val,key))
	src = (os.path.join(dataDir,key))
	source = os.path.abspath(src)
	destination = os.path.abspath(dest)
	escape="\'"
	print("Source:\t {skip}{src}{skip}".format(skip=escape,src=source))
	print("Destination:\t {skip}{dest}{skip}".format(skip=escape,dest=destination))
	#print("Desitnation:", '\'destination\'')

	link_cmd = "ln -nsf "
	create_symlink = "{link_cmd} {source} {destination}".format(link_cmd=link_cmd,source=source,destination=destination)
	print("****** NOW RUNNING COMMAND ******: " + create_symlink)
	print (run_cmd(create_symlink))

###########################################
def peSamples(pe_dir):
	pe_read_list=os.listdir(pairedendDir)

	sample_list=[]
	for read in pe_read_list:
		pe_allowed_extn=["_R1.fq","_R1.fastq","_R1.fq.gz","_R1.fastq.gz"]
		if any (read.endswith(ext) for ext in pe_allowed_extn):
			
			sample_name=read.split("_R1.f",1)[0]
			sample_list.append(sample_name)

			file_extn=read.split('.',1)[1]

	with open ((os.path.join(os.getcwd(),"config",'pe_samples.lst')),'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()

	return file_extn

def seSamples(se_dir):
	se_read_list=os.listdir(singleendDir)
	sample_list=[]
	for read in se_read_list:
		se_allowed_extn=[".fq",".fastq",".fq.gz",".fastq.gz"]
		
		if any (read.endswith(ext) for ext in se_allowed_extn):
			
			sample_name=read.split(".",1)[0]
			sample_list.append(sample_name)

			file_extn=read.split('.',1)[1]


	with open ((os.path.join(os.getcwd(),"config",'se_samples.lst')),'w') as file:
		for sample in sample_list:
			file.write("%s\n" % sample)
	file.close()

	return file_extn

################################################################################
#Generate target file

def peTrget(pe_dir):
	pe_read_list=os.listdir(pairedendDir)

	sample_list=[]
	for read in pe_read_list:
		pe_allowed_extn=["_R1.fq","_R1.fastq","_R1.fq.gz","_R1.fastq.gz"]
		if any (read.endswith(ext) for ext in pe_allowed_extn):
			
			sample_name=read.split("_R1.f",1)[0]
			sample_list.append(sample_name)

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
		
		
	else:
		print('Target file generated\n')
		df.to_csv("config" +"/"+"target.tsv", sep='\t', index=False)
		df.to_csv("config" +"/"+"target.tsv", sep='\t', index=False)
		sampledf.to_csv("config" + "/" + "samples.txt", index=False, header=False)
		groupdf.to_csv("config" + "/" + "group.tsv", sep='\t', index=False, header=False)

	

def seTarget(se_dir):
	se_read_list=os.listdir(singleendDir)
	sample_list=[]
	for read in se_read_list:
		se_allowed_extn=[".fq",".fastq",".fq.gz",".fastq.gz"]
		
		if any (read.endswith(ext) for ext in se_allowed_extn):
			
			sample_name=read.split(".",1)[0]
			sample_list.append(sample_name)

			file_extn=read.split('.',1)[1]


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
		
		else:
			print('Target file generated\n')
			df.to_csv("config" +"/"+"target.tsv", sep='\t', index=False)
			df.to_csv("config" +"/"+"target.tsv", sep='\t', index=False)
			sampledf.to_csv("config" + "/" + "samples.txt", index=False, header=False)
			groupdf.to_csv("config" + "/" + "group.tsv", sep='\t', index=False, header=False)
################################################################################

def genome(genome_dir):
	genome_file_list=os.listdir(genomeDir)
	for genome in genome_file_list:
		genome_allowed_extn=[".fasta",".fa",".fna",".fas"]
		
		if any (genome.endswith(ext) for ext in genome_allowed_extn):
			
			genome_name,file_extn=genome.split(".",1)[0],genome.split('.',1)[1]
	
	return genome_name,file_extn

def transcriptome(transcriptome_dir):
	transcriptome_file_list=os.listdir(transcriptomeDir)
	for transcriptome in transcriptome_file_list:
		transcriptome_allowed_extn=[".fasta",".fa",".ffn",".fas"]
		
		if any (transcriptome.endswith(ext) for ext in transcriptome_allowed_extn):
			
			transcriptome_name,file_extn=transcriptome.split(".",1)[0],transcriptome.split('.',1)[1]
	
	return transcriptome_name,file_extn


def annotation(annotation_dir):
	annotation_file_list=os.listdir(annotationDir)
	for annotation in annotation_file_list:
		annotation_allowed_extn=[".gff",".gtf"]
		
		if any (annotation.endswith(ext) for ext in annotation_allowed_extn):
			
			annotation_name,file_extn=annotation.split(".",1)[0],annotation.split('.',1)[1]
	
	return annotation_name,file_extn



genomeDir=os.path.abspath(os.path.join(outDir, "geno"))
transcriptomeDir=os.path.abspath(os.path.join(outDir, "tran"))
annotationDir=os.path.abspath(os.path.join(outDir, "anot"))
	
if readType=="pe":
	readExtn=peSamples(pairedendDir)
	#print('pe ext',pe_ext)
	peTrget(pairedendDir)
if readType=="se":
	readExtn=seSamples(singleendDir)
	#print('se ext',se_ext)
	seTarget(singleendDir)

if not (len(os.listdir(genomeDir))==0):
	genome=genome(genomeDir)
	genomeName=genome[0]
	genomeExtn=genome[1]
	print('genome name: ',genomeName)
	print('genome extension: ',genomeExtn)


if (len(os.listdir(genomeDir))==0):
	genomeName="NA"
	genomeExtn="NA"
	print('genome name: ',genomeName)
	print('genome extension: ',genomeExtn)



if not (len(os.listdir(transcriptomeDir))==0):
	transcriptome=transcriptome(transcriptomeDir)
	transcriptomeName=transcriptome[0]
	transcriptomeExtn=transcriptome[1]
	print('transcriptome name: ',transcriptomeName)
	print('transcriptome extension: ', transcriptomeExtn)

if (len(os.listdir(transcriptomeDir))==0):
	transcriptomeName="NA"
	transcriptomeExtn="NA"
	print('transcriptome name: ',transcriptomeName)
	print('transcriptome extension: ', transcriptomeExtn)


if not (len(os.listdir(annotationDir))==0):
	annotation=annotation(annotationDir)
	annotationName=annotation[0]
	annotationExtn=annotation[1]
	print('annotation name: ', annotationName)
	print('annotation extension: ', annotationExtn)


if (len(os.listdir(annotationDir))==0):
	annotationName="NA"
	annotationExtn="NA"
	print('annotation name: ', annotationName)
	print('annotation extension: ', annotationExtn)





with open('luigi.cfg', 'w') as config:
	config.write('[core]\n')
	config.write('default-scheduler-port:{port}\n'.format(port=port))
	config.write('error-email={email}\n\n'.format(email=email))
	config.write('[GlobalParameter]\n')

	config.write('project_name={projectDir}/\n'.format(projectDir=projectDir))
	config.write('adapter={adapter}\n'.format(adapter=adapter))

	config.write('domain={domain}\n'.format(domain=domain))

	if domain=="prokaryote":
		config.write('feature_type=CDS\n')
	if domain=="eukaryote":
		config.write('feature_type=exon\n')

	
	if not (len(os.listdir(genomeDir))==0):	
		config.write('genome_dir={genomeDir}/\n'.format(genomeDir=genomeDir))
		config.write('genome_name={genomeName}\n'.format(genomeName=genomeName))
		config.write('genome_suffix={genomeExtn}\n'.format(genomeExtn=genomeExtn))

	if (len(os.listdir(genomeDir))==0):	
		config.write('genome_dir=NA\n')
		config.write('genome_name=NA\n')
		config.write('genome_suffix=NA\n')	
	

	
	config.write('read_suffix={readExtn}\n'.format(readExtn=readExtn))


	if not (len(os.listdir(singleendDir))==0):
		config.write('rnaseq_dir={singleendDir}/\n'.format(singleendDir=singleendDir))
		config.write('read_library_type=se\n')

	if not (len(os.listdir(pairedendDir))==0):
		config.write('rnaseq_dir={pairedendDir}/\n'.format(pairedendDir=pairedendDir))
		config.write('read_library_type=pe\n')

	if (len(os.listdir(transcriptomeDir))==0):
		config.write('transcriptome_dir={transcriptomeDir}/\n'.format(transcriptomeDir=transcriptomeDir))
		config.write('transcriptome_name={genomeName}\n'.format(genomeName=genomeName))
		config.write('transcriptome_suffix=ffn\n')
	else:
		config.write('transcriptome_dir={transcriptomeDir}/\n'.format(transcriptomeDir=transcriptomeDir))
		config.write('transcriptome_name={transcriptomeName}\n'.format(transcriptomeName=transcriptomeName))
		config.write('transcriptome_suffix={transcriptomeExtn}\n'.format(transcriptomeExtn=transcriptomeExtn))
	
	if (len(os.listdir(annotationDir))==0):
		config.write('annotation_dir={annotationDir}/\n'.format(annotationDir=annotationDir))
		config.write('annotation_name={genomeName}\n'.format(genomeName=genomeName))
		config.write('annotation_suffix=NA\n')

	else:
		config.write('annotation_dir={annotationDir}/\n'.format(annotationDir=annotationDir))
		config.write('annotation_name={annotationName}\n'.format(annotationName=annotationName))
		config.write('annotation_suffix={annotationExtn}\n'.format(annotationExtn=annotationExtn))
	
	config.write('threads={cpu}\n'.format(cpu=threads)) 
	config.write('maxMemory={memory}\n'.format(memory=maxMemory))
	config.close()



if not (len(os.listdir(pairedendDir))==0):
	read_library_type="pe"

if not (len(os.listdir(singleendDir))==0):
	read_library_type="se"



with open('run_alignment_based_dea_with_read_cleaning.sh', 'w') as script:
	script.write('python rnaseq.py alignmentBasedDEA \\\n--pre-process-reads yes \\\n--bbduk-quality-trim f \\\n--bbduk-trim-quality 6 \\\n--bbduk-min-length 50 \\\n--bbduk-trim-front 0 \\\n--bbduk-trim-tail 0 \\\n--bbduk-min-average-quality 10 \\\n--bbduk-mingc 0.0 \\\n--bbduk-maxgc  1.0 \\\n--rnaseq-aligner star \\\n--dea-method deseq2 \\\n--reference-condition control --local-scheduler')
	script.close()

with open('run_alignment_based_dea_without_read_cleaning.sh', 'w') as script:
	script.write('python rnaseq.py alignmentBasedDEA \\\n--pre-process-reads no \\\n--rnaseq-aligner star \\\n--dea-method deseq2 \\\n--reference-condition control --local-scheduler')
	script.close()

with open('run_alignment_free_dea_with_read_cleaning.sh', 'w') as script:
	script.write('python rnaseq.py alignmentFreeDEA \\\n--pre-process-reads yes \\\n--bbduk-quality-trim f \\\n--bbduk-trim-quality 6 \\\n--bbduk-min-length 50 \\\n--bbduk-trim-front 0 \\\n--bbduk-trim-tail 0 \\\n--bbduk-min-average-quality 10 \\\n--bbduk-mingc 0.0 \\\n--bbduk-maxgc  1.0 \\\n--quant-method salmon \\\n--dea-method deseq2 \\\n--reference-condition control --local-scheduler')
	script.close()

with open('run_alignment_free_dea_without_read_cleaning.sh', 'w') as script:
	script.write('python rnaseq.py alignmentFreeDEA \\\n--pre-process-reads no \\\n--quant-method salmon \\\n--dea-method deseq2 \\\n--reference-condition control --local-scheduler')
	script.close()
