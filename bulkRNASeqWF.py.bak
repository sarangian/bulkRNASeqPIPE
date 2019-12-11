#!/usr/bin/env python3
from utility import utility
import time
import pandas as pd
import json
import re
import luigi
import os
import subprocess
import logging
import enum
import shutil


########################################################################################################################
#logger = logging.getLogger('luigi-interface')

class GlobalParameter(luigi.Config):
	basefolder = luigi.Parameter(default=(os.getcwd()) + "/")
	logfolder = luigi.Parameter(default=(os.getcwd()) + '/log/')
	utilityFolder = luigi.Parameter(default=(os.getcwd()) + '/utility/')


def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
						 shell=True,
						 universal_newlines=True,
						 stdout=subprocess.PIPE,
						 executable='/bin/bash')
	output = p.communicate()[0]
	return output

def prepare_sample_list(input_file):
	df = pd.read_table(input_file, sep='\t+', engine='python')
	df1 = df[["samples"]]
	samples = df1.to_csv('samples.txt',header=False,index=False)

	return samples

########################################################################################################################
class Download(luigi.Task):
	accession = luigi.Parameter()

	def output(self):
		folder = GlobalParameter().basefolder + "/raw_data/genome/"

		return {'out_R1': luigi.LocalTarget(folder + self.accession + ".gff"),
				'out_R2': luigi.LocalTarget(folder + self.accession + ".genome.fna")
				}

	def run(self):
		folder = GlobalParameter().basefolder + "/raw_data/genome/"
		cmd_download = "[ -d  {folder} ] || mkdir -p {folder}; cd {folder}; wget https://dfast.nig.ac.jp/genomes/{accession}/{accession}.genome.fna && wget https://dfast.nig.ac.jp/genomes/{accession}/{accession}.gff".format(
			accession=self.accession, folder=folder)
		print("****** NOW RUNNING COMMAND ******: " + cmd_download)
		print (run_cmd(cmd_download))


########################################################################################################################
class TimeTask(object):
	'''
	A mixin that when added to a luigi task, will print out
	the tasks execution time to standard out, when the task is
	finished

	@luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
	def print_execution_time(self, processing_time):
		print('### PROCESSING TIME ###: ' + str(processing_time))'''


	@luigi.Task.event_handler(luigi.Event.PROCESSING_TIME)
	def processing_time(task, processing_time):
		print ('\nProcessing time\t task {0} : {1} seconds\n'.format(task, processing_time))
########################################################################################################################
class prepareProject(luigi.Task):
	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string [=MyProject]")
	sampleMetadataFile = luigi.Parameter(description="""The user has to supply a tab delimited file which describes the experiment, 
	i.e. which contains the name of the biological condition associated with each sample.This file has one row per sample 
	and is composed of at least three columns with headers. column 1 : sample labels --mst be unique. 
	column 2 : name of the samples --must be unique. column 3 : biological conditions associated with the samples.
	NOTE: each sample must have two replicates""")

	sampleHeader = luigi.Parameter(default="samples",
								   description='''header name of the sample column. string([=samples])''')
	groupHeader = luigi.Parameter(default="group",
								   description='''header name of the sample group column . string([=group])''')
	refCondName = luigi.Parameter(default="control",
								  description='''Name of the reference condition . string([=control])''')


	def output(self):
		projectFolder = os.path.join(GlobalParameter().basefolder, self.projectName + "/")
		return {'out': luigi.LocalTarget(projectFolder + "samples.txt"),
				'out2': luigi.LocalTarget(projectFolder + "group.tsv")}

	def run(self):
		basefolder = GlobalParameter().basefolder
		projectFolder = os.path.join(GlobalParameter().basefolder, self.projectName + "/")
		metadataFile = "{basefolder}{metadataFile}".format(basefolder=basefolder, metadataFile=self.sampleMetadataFile)

		cmd_sanity_check = "[ -d  {projectFolder} ] || mkdir -p {projectFolder}; cd {projectFolder}; " \
								   "checkTarget.R -t {metadataFile} " \
								   "-s {sampleHeader} " \
								   "-v {groupHeader} " \
								   "-c {refCondName} ".\
			format(metadataFile=metadataFile,
				   sampleHeader=self.sampleHeader,
				   groupHeader=self.groupHeader,
				   refCondName=self.refCondName,
				   projectFolder=projectFolder)

		print("****** NOW RUNNING COMMAND ******: " + cmd_sanity_check)
		print (run_cmd(cmd_sanity_check))


########################################################################################################################
class rawSamplesQC(luigi.Task, TimeTask):
	sampleName = luigi.Parameter()
	readType = luigi.Parameter()
	projectName = luigi.Parameter()
	parameter = luigi.Parameter(default="-t 8")

	def output(self):

		raw_read_qc_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "QCReports", self.readType
										  + "_" + "FastQC-RawReads" + "/")

		if self.readType == "paired":
			return {'out_R1': luigi.LocalTarget(raw_read_qc_folder + self.sampleName + "_R1_fastqc.html"),
					'out_R2': luigi.LocalTarget(raw_read_qc_folder + self.sampleName + "_R2_fastqc.html")
					}
		if self.readType == "single":
			return {'out': luigi.LocalTarget(raw_read_qc_folder + self.sampleName + "_fastqc.html")}

	def run(self):
		raw_read_qc_folder = os.path.join(GlobalParameter().basefolder, self.projectName,"QCReports", self.readType
										  + "_" + "FastQC-RawReads")
		raw_rnaseq_read_folder = os.path.join(GlobalParameter().basefolder, "raw_data", "RNASeq/")

		cmd_raw_pe_qc = "[ -d  {raw_read_qc_folder} ] || mkdir -p {raw_read_qc_folder}; cd {raw_read_qc_folder}; fastqc " \
						"--threads 16 " \
						"{raw_rnaseq_read_folder}{sampleName}_R1.fastq.gz " \
						"{raw_rnaseq_read_folder}{sampleName}_R2.fastq.gz " \
						"-o {raw_read_qc_folder}".format(sampleName=self.sampleName,
														 raw_read_qc_folder=raw_read_qc_folder,
														 raw_rnaseq_read_folder=raw_rnaseq_read_folder,
														 parameter=self.parameter)

		cmd_raw_se_qc = "[ -d  {raw_read_qc_folder} ] || mkdir -p {raw_read_qc_folder}; cd {raw_read_qc_folder}; fastqc " \
						"--threads 16 " \
						"{raw_rnaseq_read_folder}{sampleName}.fastq.gz " \
						"-o {raw_read_qc_folder}".format(sampleName=self.sampleName,
														 raw_read_qc_folder=raw_read_qc_folder,
														 raw_rnaseq_read_folder=raw_rnaseq_read_folder,
														 parameter=self.parameter)

		if self.readType == "paired":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pe_qc)
			print (run_cmd(cmd_raw_pe_qc))


		else:
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_se_qc)
			print (run_cmd(cmd_raw_se_qc))

#######################################################################################################################
class rawReadsQC(luigi.Task):
	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string [=MyProject]")

	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")




	def requires(self):

		return [rawSamplesQC(sampleName=i,
							readType=self.readType,
							projectName=self.projectName)
					for i in [line.strip()
						for line in open((os.path.join(GlobalParameter().basefolder, self.projectName,"samples.txt")))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))



########################################################################################################################
class bbduk(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string [=MyProject]")
	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	parameter = luigi.Parameter(default="ktrim=r k=21 tbo=t tpe=t hdist=1 mink=11 qtrim=w trimq=20")
	Xms = luigi.Parameter(default="2",description="JVM will be started with Xms amount (in GB) of memory. Default: -Xms 2")
	Xmx = luigi.Parameter(default="20",description="JVM will be able to use a maximum of Xmx (in GB) amount of memory. Default -Xmx 20")
	cpu = luigi.Parameter(default="2", description="worker thread number. Default: -n_cpu 2")
	minL = luigi.Parameter(default="50",description="reads shorter than min_length will be discarded. Default: "
													"min_length=50")
	minAQ = luigi.Parameter(default="10", description="Reads with average quality (after trimming) below "
																"this will be discarded. Default: min_avg_quality=10")
	minBQ = luigi.Parameter(default="0",description="Reads with any base below this quality (after trimming) will be "
												   "discarded. Default: min_base_quality=0")
	minGC = luigi.Parameter(default="0.0", description="Discard reads with GC content below this. Default: min_gc=0.001")
	maxGC = luigi.Parameter(default="1.0", description="Discard reads with GC content below this. Default: max_gc=0.999")
	kmer = luigi.Parameter(default="13", description="Kmer length used for finding contaminants. Default: kmer=13")
	minQ = luigi.Parameter(default="10", description="Regions with average quality below minQ will be trimmed. "
													 "Default: minQ=10")
	trimF = luigi.Parameter(default="0", description="trimming how many bases in front for read. Default: "
														  "trim_front=0")
	trimT = luigi.Parameter(default="0", description="trimming how many bases in tail for read. Default: "
														  "trim_tail=0")



	def output(self):

		clean_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "CleanedReads", "Cleaned_" +
										 self.readType + "-end_Reads/")

		if self.readType == "paired":

			return {'out1': luigi.LocalTarget(clean_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(clean_read_folder + self.sampleName + "_R2.fastq")}
		else:
			return {'out': luigi.LocalTarget(clean_read_folder + self.sampleName + ".fastq")}

	def run(self):

		bbduk_read_stat_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "QCReports",
											  "BBDuk-CleanStats", "Cleaned_" + self.readType + "-end_Reads/")
		clean_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "CleanedReads", "Cleaned_" +
										 self.readType + "-end_Reads/")
		raw_rnaseq_read_folder = os.path.join(GlobalParameter().basefolder, "raw_data", "RNASeq/")

		cmd_clean_pe = "[ -d  {clean_read_folder} ] || mkdir -p {clean_read_folder}; " \
					   "mkdir -p {bbduk_read_stat_folder}; bbduk.sh " \
					   "-Xms{Xms}g " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "minlength={minL} " \
					   "minavgquality={minAQ} " \
					   "minbasequality={minBQ} " \
					   "trimq={minQ} " \
					   "ftl={trimF} " \
					   "ftr2={trimT} " \
					   "mingc={minGC} " \
					   "maxgc={maxGC} " \
					   "in1={raw_rnaseq_read_folder}{sampleName}_R1.fastq.gz " \
					   "in2={raw_rnaseq_read_folder}{sampleName}_R2.fastq.gz " \
					   "out={clean_read_folder}{sampleName}_R1.fastq " \
					   "out2={clean_read_folder}{sampleName}_R2.fastq " \
					   "ziplevel=9 ref={adapter} " \
					   "stats={bbduk_read_stat_folder}{sampleName}.stat " \
					   "bqhist={bbduk_read_stat_folder}{sampleName}_Qual.hist " \
					   "gchist={bbduk_read_stat_folder}{sampleName}_GC.hist " \
					   "{parameter} ".format(sampleName=self.sampleName,
											 projectName=self.projectName,
											 adapter=self.adapter,
											 clean_read_folder=clean_read_folder,
											 Xms=self.Xms,
											 Xmx=self.Xmx,
											 cpu=self.cpu,
											 minL=self.minL,
											 minAQ=self.minAQ,
											 minBQ=self.minBQ,
											 minQ=self.minQ,
											 trimF=self.trimF,
											 trimT=self.trimT,
											 minGC=self.minGC,
											 maxGC=self.maxGC,
											 kmer=self.kmer,
											 raw_rnaseq_read_folder=raw_rnaseq_read_folder,
											 bbduk_read_stat_folder=bbduk_read_stat_folder,
											 parameter=self.parameter)

		cmd_clean_se = "[ -d  {clean_read_folder} ] || mkdir -p {clean_read_folder};bbduk.sh " \
					   "-Xms{Xms}g " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "minlength={minL} " \
					   "minavgquality={minAQ} " \
					   "minbasequality={minBQ} " \
					   "mingc={minGC} " \
					   "maxgc={maxGC} " \
					   "trimq={minQ} " \
					   "ftl={trimF} " \
					   "ftr2={trimT} " \
					   "in={raw_rnaseq_read_folder}{sampleName}.fastq.gz " \
					   "out={clean_read_folder}{sampleName}.fastq " \
					   "ziplevel=9 ref={adapter} " \
					   "stats={bbduk_read_stat_folder}{sampleName}.stat " \
					   "bqhist={bbduk_read_stat_folder}{sampleName}_Qual.hist " \
					   "gchist={bbduk_read_stat_folder}{sampleName}_GC.hist " \
					   "{parameter} ".format(sampleName=self.sampleName,
											 projectName=self.projectName,
											 adapter=self.adapter,
											 clean_read_folder=clean_read_folder,
											 Xms=self.Xms,
											 Xmx=self.Xmx,
											 minQ=self.minQ,
											 trimF=self.trimF,
											 trimT=self.trimT,
											 cpu=self.cpu,
											 minL=self.minL,
											 minAQ=self.minAQ,
											 minBQ=self.minBQ,
											 minGC=self.minGC,
											 maxGC=self.maxGC,
											 kmer=self.kmer,
											 raw_rnaseq_read_folder=raw_rnaseq_read_folder,
											 bbduk_read_stat_folder=bbduk_read_stat_folder,
											 parameter=self.parameter)

		if self.readType == "paired":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pe)
			print (run_cmd(cmd_clean_pe))

		else:
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_se)
			print (run_cmd(cmd_clean_se))

########################################################################################################################
class cleanedSampleQC(luigi.Task):
	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string [=MyProject]")
	sampleName = luigi.Parameter(description="name of the sampleName to be analyzed. (string)")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")
	parameter = luigi.Parameter(default="-t 16")

	def requires(self):

			return [rawSamplesQC(sampleName=self.sampleName,
							   readType=self.readType,
							   projectName=self.projectName),

						bbduk(sampleName=self.sampleName,
							   adapter=self.adapter,
							   readType=self.readType,
							   projectName=self.projectName)
					]


	def output(self):
		clean_read_qc_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "QCReports",self.readType + "_" + "FastQC-ProcessedReads/")


		if self.readType == "paired":
			return {'out_pe_R1': luigi.LocalTarget(clean_read_qc_folder + self.sampleName + "_R1_fastqc.html"),
					'out_pe_R2': luigi.LocalTarget(clean_read_qc_folder + self.sampleName + "_R2_fastqc.html")}
		if self.readType == "single":
			return {'out': luigi.LocalTarget(clean_read_qc_folder + self.sampleName + "_fastqc.html")}

	def run(self):
		clean_read_qc_folder = os.path.join(GlobalParameter().basefolder,self.projectName,"QCReports",self.readType +
											"_" + "FastQC-ProcessedReads/")
		clean_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "CleanedReads", "Cleaned_" + self.readType + "-end_Reads/")

		cmd_clean_pe_qc = "[ -d  {clean_read_qc_folder} ] || mkdir -p {clean_read_qc_folder}; " \
						  "cd {clean_read_qc_folder}; fastqc " \
						  "-t 2 " \
						  "{clean_read_folder}{sampleName}_R1.fastq " \
						  "{clean_read_folder}{sampleName}_R2.fastq " \
						  "-o {clean_read_qc_folder}".format(sampleName=self.sampleName,
													   clean_read_qc_folder=clean_read_qc_folder,
													   clean_read_folder=clean_read_folder,
													   parameter=self.parameter)

		cmd_clean_se_qc = "[ -d  {clean_read_qc_folder} ] || mkdir -p {clean_read_qc_folder}; " \
						  "cd {clean_read_qc_folder}; " \
						  "fastqc -t 16 " \
						  "{clean_read_folder}{sampleName}.fastq " \
						  "-o {clean_read_qc_folder}".format(sampleName=self.sampleName,
															 clean_read_qc_folder=clean_read_qc_folder,
															 clean_read_folder=clean_read_folder,
															 parameter=self.parameter)

		if self.readType == "paired":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pe_qc)
			print (run_cmd(cmd_clean_pe_qc))


		else:
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_se_qc)
			print (run_cmd(cmd_clean_se_qc))


##########################################################################################################################
class preProcessSamples(luigi.Task, TimeTask):
	'''sampleListFile = luigi.Parameter(default="samples.txt",
									 description="path to the file containing name of the samples to be analyzed. "
												 "one sampleName name per row. (string [=samples.txt])")'''
	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string [=MyProject]")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	def requires(self):
		return [cleanedSampleQC(sampleName=i,
								readType=self.readType,
								adapter=self.adapter,
								projectName=self.projectName)
				for i in [line.strip()
						  for line in
						  open((os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt")))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('pre_processing_workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))


#########################################################################################################################
# Correct Reads using RCorrector

class RCorrector(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string [=MyProject]")
	sampleName = luigi.Parameter(description="name of the sampleName to be analyzed. (string)")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")
	kmerLength = luigi.Parameter(default="23",description="Length of K-mer. Must be <=32. (int [=23]")
	cpu = luigi.Parameter(default="4",description="Number of threads to be used. (int [=4]")
	parameter = luigi.Parameter(default="-t 2")

	def requires(self):
		if self.readType == "paired":
			return [bbduk(sampleName=self.sampleName,
							   adapter=self.adapter,
							   readType="paired",
							   projectName=self.projectName)]
		else:
			return [bbduk(sampleName=self.sampleName,
							   adapter=self.adapter,
							   readType="single",
							   projectName=self.projectName)]

	def output(self):
		corrected_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName,
											 "CorrectedReads", "Corrected_" + self.readType + "-end_Reads/")


		if self.readType == "paired":

			return {'out_pe_R1': luigi.LocalTarget(corrected_read_folder +"/" + self.sampleName + "_R1.fastq"),
					'out_pe_R2': luigi.LocalTarget(corrected_read_folder +"/" + self.sampleName + "_R2.fastq")
					}

		else:
			return {'out_se': luigi.LocalTarget(corrected_read_folder + self.sampleName + ".fastq")}

	def run(self):
		corrected_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName,
											 "CorrectedReads", "Corrected_" + self.readType + "-end_Reads/")
		clean_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "CleanedReads", "Cleaned_" + self.readType + "-end_Reads/")



		cmd_correct_pe = "[ -d  {corrected_read_folder} ] || mkdir -p {corrected_read_folder}; " \
						 "cd {corrected_read_folder}; run_rcorrector.pl " \
						 "-1 {clean_read_folder}{sampleName}_R1.fastq " \
						 "-2 {clean_read_folder}{sampleName}_R2.fastq " \
						 "-od {corrected_read_folder} " \
						 "-t 2".format(sampleName=self.sampleName,
									   corrected_read_folder=corrected_read_folder,
									   clean_read_folder=clean_read_folder,
									   parameter=self.parameter)

		cmd_filter_uncorrected_pe = "cd {corrected_read_folder}; filter_rcorrector_pe.py " \
									"-1 {corrected_read_folder}{sampleName}_R1.cor.fq " \
									"-2 {corrected_read_folder}{sampleName}_R2.cor.fq " \
									"-o filter".format(sampleName=self.sampleName,
													   corrected_read_folder=corrected_read_folder,
													   parameter=self.parameter)

		cmd_correct_se = "[ -d  {corrected_read_folder} ] || mkdir -p {corrected_read_folder}; " \
						 "cd {corrected_read_folder}; run_rcorrector.pl " \
						 "-s {clean_read_folder}{sampleName}.fastq " \
						 "-od {corrected_read_folder} " \
						 "-t 2".format(sampleName=self.sampleName,
									   corrected_read_folder=corrected_read_folder,
									   clean_read_folder=clean_read_folder,
									   parameter=self.parameter)

		cmd_filter_uncorrected_se = "cd {corrected_read_folder}; filter_rcorrector_se.py " \
									"-in {corrected_read_folder}{sampleName}.cor.fq " \
									"-o filter".format(sampleName=self.sampleName,
													   corrected_read_folder=corrected_read_folder,
													   parameter=self.parameter)

		if self.readType == "paired":

			print("****** NOW RUNNING COMMAND ******: " + cmd_correct_pe)
			print (run_cmd(cmd_correct_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_filter_uncorrected_pe)
			print (run_cmd(cmd_filter_uncorrected_pe))

			cmd_rename_forward = "cd {corrected_read_folder}; " \
								 "mv {corrected_read_folder}filter_{sampleName}_R1.cor.fq {corrected_read_folder}{sampleName}_R1.fastq" \
				.format(sampleName=self.sampleName,
						corrected_read_folder=corrected_read_folder)

			print("****** NOW RUNNING COMMAND ******: " + cmd_rename_forward)
			print (run_cmd(cmd_rename_forward))

			cmd_rename_reverse = "cd {corrected_read_folder}; " \
								 "mv {corrected_read_folder}filter_{sampleName}_R2.cor.fq {corrected_read_folder}{sampleName}_R2.fastq" \
				.format(sampleName=self.sampleName,
						corrected_read_folder=corrected_read_folder)

			print("****** NOW RUNNING COMMAND ******: " + cmd_rename_reverse)
			print (run_cmd(cmd_rename_reverse))

			cmd_clean = "cd {corrected_read_folder}; rm *.cor.fq ".format(
				corrected_read_folder=corrected_read_folder)
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean)
			print (run_cmd(cmd_clean))


		if self.readType == "single":

			print("****** NOW RUNNING COMMAND ******: " + cmd_correct_se)
			print (run_cmd(cmd_correct_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_filter_uncorrected_se)
			print (run_cmd(cmd_filter_uncorrected_se))

			cmd_rename_single = "cd {corrected_read_folder}; " \
								"mv {corrected_read_folder}filter_{sampleName}.cor.fq {corrected_read_folder}{sampleName}.fastq" \
				.format(sampleName=self.sampleName,
						corrected_read_folder=corrected_read_folder)

			print("****** NOW RUNNING COMMAND ******: " + cmd_rename_single)
			print (run_cmd(cmd_rename_single))

			cmd_clean = "cd {corrected_read_folder}; rm *.cor.qc".format(corrected_read_folder=corrected_read_folder)

			print("****** NOW RUNNING COMMAND ******: " + cmd_clean)
			print (run_cmd(cmd_clean))



##########################################################################################################################
##########################################################################################################################
class correctReads(luigi.Task):
	'''sampleListFile = luigi.Parameter(default="samples.txt",
									 description="path to the file containing name of the samples to be analyzed. "
												 "one sampleName name per row. (string [=samples.txt])")'''
	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string [=MyProject]")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")



	def requires(self):
		if self.readType == "paired":
			return [RCorrector(sampleName=i,
								readType="paired",
								adapter=self.adapter,
								projectName=self.projectName)
					for i in [line.strip()
							  for line in
							  open((os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt")))]]
		else:
			return [RCorrector(sampleName=i,
								readType="single",
								adapter=self.adapter,
								projectName=self.projectName)
					for i in [line.strip()
							  for line in
							  open((os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt")))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))

##########################################################################################################################
# Run Rockhopper De-novo Assembler
class denovoTransAssemble(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string ["
											  "=MyProject])")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	parameter = luigi.Parameter(default="java -Xmx40G -u 500")

	domain = luigi.Parameter(description="domain of the organism . (string [=prokaryote] OR [=eukaryote])")

	assembler = luigi.Parameter(description="transcript assembly tool to be used. (string (if domain =prokaryote: ["
											"=rockhopper2] if domain=eukaryote [=trinity] OR [=spades])")

	def requires(self):
		return [correctReads(projectName=self.projectName,
							readType=self.readType,
							adapter=self.adapter)]

	def output(self):

		assembled_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "DenovoAssembly",
												   self.assembler + "_DenovoAssembly_" + self.readType + "/")
		salmon_index_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   self.assembler + "_DenovoSalmonIndex_" + self.readType + "/")


		if all ([self.readType == "paired" , self.domain == "prokaryote", self.assembler == "rockhopper2"]):
			return {'out_1': luigi.LocalTarget(assembled_transcript_folder + "transcripts.txt"),
					'out_2': luigi.LocalTarget(assembled_transcript_folder + "transcripts.fna")
					}
		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2"]):
			return {'out_1': luigi.LocalTarget(assembled_transcript_folder + "transcripts.fna"),
					'out_2': luigi.LocalTarget(assembled_transcript_folder + "transcripts.fna")
					}

		if all ([self.readType == "paired" , self.domain == "eukaryote", self.assembler == "trinity"]):
			return {'out': luigi.LocalTarget(assembled_transcript_folder + "Trinity.fasta")}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity"]):
			return {'out': luigi.LocalTarget(assembled_transcript_folder + "Trinity.fasta")}


		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades"]):
			return {'out_1': luigi.LocalTarget(assembled_transcript_folder + "soft_filtered_transcripts.fasta"),
					'out_2': luigi.LocalTarget(assembled_transcript_folder + "hard_filtered_transcripts.fasta"),
					'out_3': luigi.LocalTarget(salmon_index_folder + "hash.bin")
					}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades"]):
			return {'out_1': luigi.LocalTarget(assembled_transcript_folder + "soft_filtered_transcripts.fasta"),
					'out_2': luigi.LocalTarget(assembled_transcript_folder + "hard_filtered_transcripts.fasta")
					}

	def run(self):

		basefol = GlobalParameter().basefolder
		corrected_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName,
											 "CorrectedReads", "Corrected_" + self.readType + "-end_Reads/")
		salmon_index_folder = os.path.join(GlobalParameter().basefolder,self.projectName,"TransQuant", self.assembler
										   + "_DenovoSalmonIndex_" + self.readType + "/")
		assembled_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName,"DenovoAssembly",
												   self.assembler + "_DenovoAssembly_" + self.readType + "/")

		#Prokaryotic Denovo assembly using Rockhopper 2
		input_group_file = ((os.path.join(GlobalParameter().basefolder, self.projectName, "group.tsv")))

		#Command to generate Rockhopper Paired-end input
		cmd_rockhopper_pe_read_input = utility.prepare_rockhopper_pe_input(input_group_file)

		#Command to run Rockhopper
		cmd_run_rockhopper_rb_pe = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
								   "cd {corrected_read_folder};" \
								   "java -Xmx12G " \
								   "-cp $(which Rockhopper.jar) Rockhopper  {rockhopper_pe_read_input} " \
								   "-s false " \
								   "-p 2 " \
								   "-o {assembled_transcript_folder}" \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					corrected_read_folder=corrected_read_folder,
					rockhopper_pe_read_input=cmd_rockhopper_pe_read_input,
					parameter=self.parameter)

		#Command to generate transcripts in fasta format from tsv file
		cmd_format_rockhopper_transcript_pe = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
											  "cd {assembled_transcript_folder};" \
											  "FormatRockhopTrans.py " \
											  "-i {assembled_transcript_folder}transcripts.txt " \
											  "-l 500 > {assembled_transcript_folder}transcripts.fna"\
			.format(assembled_transcript_folder=assembled_transcript_folder)



		cmd_run_salmom_index_pe = "[ -d  {salmon_index_folder} ] || mkdir -p {salmon_index_folder}; " \
								  "salmon index " \
								  "-t {assembled_transcript_folder}transcripts.fna " \
								  "-i {salmon_index_folder} " \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					salmon_index_folder=salmon_index_folder)



		# Code for single end data
		# Command to generate Rockhopper Paired-end input
		cmd_rockhopper_se_read_input = utility.prepare_rockhopper_se_input(input_group_file)

		# Command to run Rockhopper
		cmd_run_rockhopper_rb_se = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
								   "cd {corrected_read_folder};" \
								   "java -Xmx12G " \
								   "-cp $(which Rockhopper.jar) Rockhopper  {rockhopper_se_read_input} " \
								   "-s false " \
								   "-p 2 " \
								   "-o {assembled_transcript_folder}" \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					corrected_read_folder=corrected_read_folder,
					rockhopper_se_read_input=cmd_rockhopper_se_read_input,
					parameter=self.parameter)

		cmd_format_rockhopper_transcript_se = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
											  "cd {assembled_transcript_folder};" \
											  "FormatRockhopTrans.py " \
											  "-i {assembled_transcript_folder}transcripts.txt " \
											  "-l 500 > {assembled_transcript_folder}transcripts.fna".\
			format(assembled_transcript_folder=assembled_transcript_folder)


		cmd_run_salmom_index_se = "[ -d  {salmon_index_folder} ] || mkdir -p {salmon_index_folder}; " \
								  "salmon index " \
								  "-t {assembled_transcript_folder}transcripts.fna " \
								  "-i {salmon_index_folder} " \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					salmon_index_folder=salmon_index_folder)


 #Eukaryotic Denovo Transcriptome Assembly using Rockhopper2

		if all([self.readType == "paired", self.domain == "prokaryote", self.assembler == "rockhopper2"]):


			#print("****** NOW RUNNING COMMAND ******: " + cmd_rockhopper_pe_read_input)
			#print run_cmd(cmd_rockhopper_pe_read_input)

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_rockhopper_rb_pe)
			print (run_cmd(cmd_run_rockhopper_rb_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_format_rockhopper_transcript_pe)
			print (run_cmd(cmd_format_rockhopper_transcript_pe))

			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmom_index_pe)
			print (run_cmd(cmd_run_salmom_index_pe))


		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2"]):
			#print("****** NOW RUNNING COMMAND ******: " + cmd_rockhopper_se_read_input)
			#print run_cmd(cmd_rockhopper_se_read_input)

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_rockhopper_rb_se)
			print (run_cmd(cmd_run_rockhopper_rb_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_format_rockhopper_transcript_se)
			print (run_cmd(cmd_format_rockhopper_transcript_se))

			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmom_index_se)
			print (run_cmd(cmd_run_salmom_index_se))




########################################################################################################################
#Eukaryotic Denovo Transcriptome Assembly using Trinity                                                                #
########################################################################################################################
		input_sample_file = (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))
		cmd_trinity_pe_input = utility.prepare_trinity_pe_input(input_sample_file)
		cmd_trinity_se_input = utility.prepare_trinity_se_input(input_sample_file)


		cmd_run_trinity_pe = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
						  "cd {corrected_read_folder}; " \
							 "Trinity --seqType fq " \
							 "--normalize_by_read_set " \
							 "--min_contig_length 500 " \
							 "--max_memory 40G " \
							 " {trinity_pe} " \
							 "--output {assembled_transcript_folder} " \
							 "--CPU 4 "\
			.format(assembled_transcript_folder=assembled_transcript_folder,
					corrected_read_folder=corrected_read_folder,trinity_pe=cmd_trinity_pe_input)


		cmd_run_trinity_se = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
							 "cd {corrected_read_folder};" \
							 "Trinity --seqType fq " \
							 "--normalize_by_read_set " \
							 "--min_contig_length 500 " \
							 "--max_memory 40G " \
							 " {trinity_se}" \
							 "--output {assembled_transcript_folder} " \
							 "--CPU 4 " \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					corrected_read_folder=corrected_read_folder, trinity_se=cmd_trinity_se_input)



		cmd_run_salmom_index_trinity = "[ -d  {salmon_index_folder} ] || mkdir -p {salmon_index_folder}; " \
								  "salmon index -t {assembled_transcript_folder}Trinity.fasta -i {salmon_index_folder} " \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					salmon_index_folder=salmon_index_folder)

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "trinity"]):
			#print("****** NOW RUNNING COMMAND ******: " + cmd_trinity_pe_input)
			#print run_cmd(cmd_trinity_pe_input)

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmom_index_trinity)
			print (run_cmd(cmd_run_salmom_index_trinity))

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity"]):
			#print("****** NOW RUNNING COMMAND ******: " + cmd_trinity_se_input)
			#print run_cmd(cmd_trinity_se_input)

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmom_index_trinity)
			print (run_cmd(cmd_run_salmom_index_trinity))

########################################################################################################################
# Eukaryotic Denovo Transcriptome Assembly using RNASpades                                                             #
########################################################################################################################
		input_sample_file = ((os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt")))
		cmd_spades_pe_input = utility.prepare_rnaspades_pe_input(input_sample_file)
		cmd_spades_se_input = utility.prepare_rnaspades_se_input(input_sample_file)

		cmd_run_rnaspades_pe = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
							   "cd {corrected_read_folder};" \
							   "spades.py --rna " \
							   "--only-assembler " \
							   "-m 10 " \
							   " {spades_pe_input} " \
							   "-o {assembled_transcript_folder} " \
							   "-t 1 " \
			.format(corrected_read_folder=corrected_read_folder,
					assembled_transcript_folder=assembled_transcript_folder,spades_pe_input=cmd_spades_pe_input)

		cmd_run_rnaspades_se = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
							   "cd {corrected_read_folder}; " \
							   "spades.py --rna " \
							   "--only-assembler " \
							   "-m 10 " \
							   " {spades_pe_input} " \
							   "-o {assembled_transcript_folder} " \
							   "-t 1 " \
			.format(corrected_read_folder=corrected_read_folder,
					assembled_transcript_folder=assembled_transcript_folder, spades_pe_input=cmd_spades_se_input)

		cmd_run_salmom_index_spades = "[ -d  {salmon_index_folder} ] || mkdir -p {salmon_index_folder}; " \
									  "salmon index " \
									  "-t {assembled_transcript_folder}hard_filtered_transcripts.fasta " \
									  "-i {salmon_index_folder} " \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					salmon_index_folder=salmon_index_folder)

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades"]):

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_rnaspades_pe)
			print (run_cmd(cmd_run_rnaspades_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmom_index_spades)
			print (run_cmd(cmd_run_salmom_index_spades))

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades"]):

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_rnaspades_se)
			print (run_cmd(cmd_run_rnaspades_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmom_index_spades)
			print (run_cmd(cmd_run_salmom_index_spades))


#######################################################################################################################################################################################
##Run Salmom Quantification
class salmonQuant(luigi.Task):
	sampleName = luigi.Parameter()
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	adapter = luigi.Parameter(default="/opt/RNASeqPIPE/tools/bbmap/resources/adapters.fasta.gz")
	assembler = luigi.Parameter(description="""name of the assembler to be used to assemble corrected RNASeq reads . (
											 string [=rockhopper2] OR [=spades] OR [=trinity].
											 NOTE: rockhopper2 is exclusively used to assemble bacterial RNASeq reads and must
											 not be used for domain eukaryote """)
	quantParameter = luigi.Parameter(default="-p 1")

	def requires(self):
		'''return [denovoAssemble(projectName=self.projectName,
							   sampleName=self.sampleName,
							   adapter=self.adapter,
							   readType=self.readType,
							   assembler=self.assembler)]'''
		return []

	def output(self):

		salmon_quant_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   self.assembler + "_DenovoSalmonQuant_" + self.readType + "/")

		if all([self.readType == "paired", self.domain == "prokaryote", self.assembler == "rockhopper2"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "trinity"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

	def run(self):

		salmon_quant_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   self.assembler + "_DenovoSalmonQuant_" + self.readType + "/")
		salmon_map_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",self.assembler
										 + "_DenovoSalmonMAP_" + self.readType + "/")
		corrected_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName,
											 "CorrectedReads", "Corrected_" + self.readType + "-end_Reads/")
		salmon_index_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   self.assembler + "_DenovoSalmonIndex_" + self.readType + "/")

		cmd_run_salmon_quant_pe = "[ -d {salmon_quant_folder} ] || mkdir -p {salmon_quant_folder}; " \
								  "[ -d {salmon_map_folder} ] || mkdir -p {salmon_map_folder}; " \
								  "cd {salmon_map_folder}; " \
								  "salmon quant --no-version-check -p 1 " \
								  "-i {salmon_index_folder} " \
								  "-l A " \
								  "-1 {corrected_read_folder}{sampleName}_R1.fastq " \
								  "-2 {corrected_read_folder}{sampleName}_R2.fastq " \
								  "--dumpEq " \
								  "--output {salmon_quant_folder}{sampleName} " \
								  "--validateMappings " \
								  "--writeMappings | samtools view -bS - | samtools sort -m 30G " \
								  "-o {salmon_map_folder}{sampleName}.bam" \
			.format(salmon_quant_folder=salmon_quant_folder,
					salmon_map_folder=salmon_map_folder,
					sampleName=self.sampleName,
					corrected_read_folder=corrected_read_folder,
					salmon_index_folder=salmon_index_folder,
					quantParameter=self.quantParameter)

		cmd_run_salmon_quant_se = "[ -d {salmon_quant_folder} ] || mkdir -p {salmon_quant_folder}; " \
								  "mkdir -p {salmon_map_folder}; " \
								  "cd {salmon_map_folder};" \
								  "salmon quant --no-version-check -p 1 " \
								  "--seqBias " \
								  "--gcBias " \
								  "-i {salmon_index_folder} " \
								  "-l A " \
								  "-r {corrected_read_folder}{sampleName}.fastq " \
								  "--dumpEq " \
								  "--output {salmon_quant_folder}{sampleName} " \
								  "--validateMappings " \
								  "--writeMappings | samtools view -bS - | samtools sort -m 30G " \
								  "-o {salmon_map_folder}{sampleName}.bam" \
			.format(salmon_quant_folder=salmon_quant_folder,
					salmon_map_folder=salmon_map_folder,
					sampleName=self.sampleName,
					corrected_read_folder=corrected_read_folder,
					salmon_index_folder=salmon_index_folder,
					quantParameter=self.quantParameter)

		if all([self.readType == "paired", self.domain == "prokaryote", self.assembler == "rockhopper2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))



		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))



		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "trinity"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))


		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))



		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))



		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))


##########################################################################################################################

class quantifyAssembledTranscripts(luigi.Task, TimeTask):
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")
	assembler = luigi.Parameter(description="""name of the assembler to be used to assemble corrected RNASeq reads . (
											 string [=rockhopper2] OR [=spades] OR [=trinity].
											 NOTE: rockhopper2 is exclusively used to assemble bacterial RNASeq reads and must
											 not be used for domain eukaryote """)
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	def requires(self):

		return [salmonQuant(assembler=self.assembler,
							domain=self.domain,
							readType=self.readType,
							projectName=self.projectName,
							adapter=self.adapter,
							sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt")))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))

#####################################################################################################################################################3
# Run Corsetanalysis
class clusterContigs(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	assembler = luigi.Parameter(description="""name of the assembler to be used to assemble corrected RNASeq reads . (
											 string [=rockhopper2] OR [=spades] OR [=trinity].
											 NOTE: rockhopper2 is exclusively used to assemble bacterial RNASeq reads and must
											 not be used for domain eukaryote """)

	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	def requires(self):

		return [quantifyAssembledTranscripts(domain=self.domain,
									adapter=self.adapter,
									assembler=self.assembler,
									readType=self.readType,
									projectName=self.projectName),

				denovoTransAssemble(domain=self.domain,
							   assembler=self.assembler,
							   projectName=self.projectName,
							   readType=self.readType,
							   adapter=self.adapter)
				]


	def output(self):

		corset_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.assembler +
									 "_DenovoSalmonQuant_" + self.readType, "Corset/")

		return {'out1': luigi.LocalTarget(corset_folder + "clusters.txt"),
				'out2': luigi.LocalTarget(corset_folder + "counts.txt")
				}


	def run(self):
		corset_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.assembler +
									 "_DenovoSalmonQuant_" + self.readType, "Corset/")
		salmon_quant_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   self.assembler + "_DenovoSalmonQuant_" + self.readType + "/")
		assembled_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "DenovoAssembly",
												   self.assembler + "_DenovoAssembly_" + self.readType + "/")

		super_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName,"DenovoAssembly",
								self.assembler + "_DenovoAssembly_" + self.readType, "SuperTranscript/")


		# Paired end reads

		input_group_file = ((os.path.join(GlobalParameter().basefolder, self.projectName, "group.tsv")))

		# Command to generate Rockhopper Paired-end input
		cmd_corset_read_input = utility.prepare_corset_input(input_group_file)


		cmd_run_corset = "[ -d  {corset_folder} ] || mkdir -p {corset_folder}; cd {corset_folder}; " \
							"corset -D 99999999999  {corset_read_input} " \
							"-i salmon_eq_classes {salmon_quant_folder}*/aux_info/eq_classes.txt " \
			.format(corset_folder=corset_folder,
					corset_read_input=cmd_corset_read_input,
					salmon_quant_folder=salmon_quant_folder)

		cmd_run_supertrans = "python $(which Lace.py) " \
								"--cores 1 " \
								"{assembled_transcript_folder}transcripts.fna " \
								"{corset_folder}clusters.txt " \
								"--outputDir {super_transcript_folder} " \
			.format(super_transcript_folder=super_transcript_folder,corset_folder=corset_folder,
					assembled_transcript_folder=assembled_transcript_folder)



		print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset)
		print (run_cmd(cmd_run_corset))

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_supertrans)
		print (run_cmd(cmd_run_supertrans))
####################################################################################################
class mapReadsToGenome(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string ["
											  "=MyProject])")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	genomeName = luigi.Parameter()

	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	parameter = luigi.Parameter(default="java -Xmx40G -u 500")

	domain = luigi.Parameter(description="domain of the organism . (string [=prokaryote] OR [=eukaryote])")

	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)

	def requires(self):

		return [alignReads(readType=self.readType,
						   aligner=self.aligner,
						   genomeName=self.genomeName,
						   projectName=self.projectName,
						   domain = self.domain,
						   adapter=self.adapter,
						   sampleName=i)
				for i in [line.strip()
						  for line in
								open (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))]]

	def output(self):
		mapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.aligner + "_map_" + self.readType + "/")

		#STAR OUTPUT
		if all([self.readType == "paired", self.aligner == "star",self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "star",self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}


		#HISAT2 OUTPUT
		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}


		#BOWTIE2 OUTPUT
		if all([self.readType == "paired", self.aligner == "bowtie2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "bowtie2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}


		#DART OUTPUT
		if all([self.readType == "paired", self.aligner == "dart", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "dart", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "dart", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}



		# SEGEMEHL OUTPUT
		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genomeName + ".bam")}

	def run(self):
		mapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.aligner + "_map_" + self.readType + "/")

		cmd_merge_bam = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "samtools merge {genomeName}.bam *.bam " \
								.format(mapFolder=mapFolder,genomeName=self.genomeName)
		print("****** NOW RUNNING COMMAND ******: " + cmd_merge_bam)
		print (run_cmd(cmd_merge_bam))

####################################################################################################
#Genome Guided Transcript Assembly
####################################################################################################
class genomeGuidedTransAssembly(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string ["
											  "=MyProject])")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	genomeName = luigi.Parameter()

	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	domain = luigi.Parameter(description="domain of the organism . (string [=prokaryote] OR [=eukaryote])")

	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)

	maxMemory = luigi.Parameter(description="Maximum Memory in GB",default="20")
	maxIntron = luigi.Parameter(description="Maximum Intron Length",default="2000")
	minContigLength = luigi.Parameter(description="Minimum Contig Length",default="200")
	threads = luigi.Parameter(description="Number of threads to be used", default="20")




	def requires(self):

		return [mapReadsToGenome(readType=self.readType,
						   aligner=self.aligner,
						   genomeName=self.genomeName,
						   projectName=self.projectName,
						   domain = self.domain,
						   adapter=self.adapter)]
						   

#Trinity --genome_guided_bam hexcentricum.bam --genome_guided_max_intron 10000 --max_memory 20G --CPU 10

	def output(self):
		gg_assembled_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "GGAssembly",
													 "trinity_GGAssembly_" + self.aligner +"_" + self.readType + "/")


		if all([self.readType == "paired", self.aligner == "star",self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "star",self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}


		#HISAT2 OUTPUT
		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}


		#BOWTIE2 OUTPUT
		if all([self.readType == "paired", self.aligner == "bowtie2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "bowtie22", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}


		#DART OUTPUT
		if all([self.readType == "paired", self.aligner == "dart", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "dart", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "dart", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}



		# SEGEMEHL OUTPUT
		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(gg_assembled_transcript_folder + "/" + "Trinity-GG.fasta")}


	def run(self):
		mapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.aligner + "_map_" + self.readType + "/")
		gg_assembled_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "GGAssembly",
													 "trinity_GGAssembly_" + self.aligner +"_" + self.readType + "/")

		gg_transcript_index_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_Bowtie2_Index_" + self.readType + "/")



		cmd_run_trinity_pe = "[ -d  {gg_assembled_transcript_folder} ] || mkdir -p {gg_assembled_transcript_folder}; " \
							 "Trinity " \
							 "--min_contig_length {minContigLength} " \
							 "--max_memory {maxMemory}G " \
							 "--genome_guided_bam {mapFolder}{genomeName}.bam " \
							 "--genome_guided_max_intron {maxIntron} " \
							 "--output {gg_assembled_transcript_folder} " \
							 "--CPU {threads} "\
							 "--include_supertranscripts "\
							 "--full_cleanup "\
			.format(gg_assembled_transcript_folder=gg_assembled_transcript_folder,
					mapFolder=mapFolder,
					minContigLength=self.minContigLength,
					maxMemory=self.maxMemory,
					genomeName=self.genomeName,
					maxIntron=self.maxIntron,
					threads=self.threads)


		cmd_run_trinity_se = "[ -d  {gg_assembled_transcript_folder} ] || mkdir -p {gg_assembled_transcript_folder}; " \
							 "Trinity " \
							 "--min_contig_length {minContigLength} " \
							 "--max_memory {maxMemory}G " \
							 "--genome_guided_bam {mapFolder}{genomeName}.bam " \
							 "--genome_guided_max_intron {maxIntron} " \
							 "--output {gg_assembled_transcript_folder} " \
							 "--CPU {threads} "\
							 "--include_supertranscripts "\
							 "--full_cleanup "\
			.format(gg_assembled_transcript_folder=gg_assembled_transcript_folder,
					mapFolder=mapFolder,
					minContigLength=self.minContigLength,
					maxMemory=self.maxMemory,
					genomeName=self.genomeName,
					maxIntron=self.maxIntron,
					threads=self.threads)


		cmd_run_gg_bowtie2_index = "[ -d  {gg_transcript_index_folder} ] || mkdir -p {gg_transcript_index_folder}; " \
							 "cd {gg_transcript_index_folder}; " \
							 "bowtie2-build {gg_assembled_transcript_folder}Trinity-GG.fasta {genomeName} " \
			.format(gg_assembled_transcript_folder=gg_assembled_transcript_folder,
				    gg_transcript_index_folder=gg_transcript_index_folder,
					genomeName=self.genomeName)




		if all([self.readType == "paired", self.aligner == "star",self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "star",self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))


		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))


		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))


		#HISAT2 OUTPUT
		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))


		#BOWTIE2 OUTPUT
		if all([self.readType == "paired", self.aligner == "bowtie2", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "bowtie22", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))


		#DART OUTPUT
		if all([self.readType == "paired", self.aligner == "dart", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "dart", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "dart", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))



		# SEGEMEHL OUTPUT
		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_pe)
			print (run_cmd(cmd_run_trinity_pe))
			print (run_cmd(cmd_run_gg_bowtie2_index))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_trinity_se)
			print (run_cmd(cmd_run_trinity_se))
			print (run_cmd(cmd_run_gg_bowtie2_index))


####################################################################################################
class mapReadToTans(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string ["
											  "=MyProject])")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	genomeName = luigi.Parameter()
	sampleName = luigi.Parameter()
	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")
	domain = luigi.Parameter(description="domain of the organism . (string [=prokaryote] OR [=eukaryote])")
	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)
	
	def requires(self):

		return [genomeGuidedTransAssembly(projectName=self.projectName,
								readType=self.readType,
								domain = self.domain,
								aligner = self.aligner,
								genomeName=self.genomeName,
								adapter=self.adapter)]


	def output(self):
		gg_transcript_map_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_bowtie2_map_" + self.readType + "/")

		return {'out1': luigi.LocalTarget(gg_transcript_map_folder + "/" + self.sampleName + ".bam")}

	def run(self):

		gg_transcript_map_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_bowtie2_map_" + self.readType + "/")

		gg_transcript_index_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_Bowtie2_Index_" + self.readType + "/")

		clean_read_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "CleanedReads", "Cleaned_" +
										 self.readType + "-end_Reads/")




		################################################################################################################
		#4 Bowtie2 Aligner
		################################################################################################################

		cmd_run_bowtie2_map_pe = "[ -d {gg_transcript_map_folder} ] || mkdir -p {gg_transcript_map_folder}; " \
								 "cd {gg_transcript_map_folder}; " \
								 "bowtie2 -x {gg_transcript_index_folder}{genomeName} " \
								 "-p 2 " \
								 "-a " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq |samtools view -bS - | samtools " \
								 "sort " \
								 "-o {gg_transcript_map_folder}{sampleName}.bam " \
			.format(gg_transcript_map_folder=gg_transcript_map_folder,
					sampleName=self.sampleName,
					genomeName=self.genomeName,
					cleanedReadFolder=clean_read_folder,
					gg_transcript_index_folder=gg_transcript_index_folder)


		cmd_run_bowtie2_map_se = "[ -d {gg_transcript_map_folder} ] || mkdir -p {gg_transcript_map_folder}; " \
								 "cd {gg_transcript_map_folder}; " \
								 "bowtie2 -x {gg_transcript_index_folder}{genomeName} " \
								 "-p 2 " \
								 "-a " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "|samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {gg_transcript_map_folder}{sampleName}.bam " \
			.format(gg_transcript_map_folder=gg_transcript_map_folder,
					sampleName=self.sampleName,
					genomeName=self.genomeName,
					cleanedReadFolder=clean_read_folder,
					gg_transcript_index_folder=gg_transcript_index_folder)


		if self.readType == "paired":
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_pe)
			print (run_cmd(cmd_run_bowtie2_map_pe))

		if self.readType == "single":
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_se)
			print (run_cmd(cmd_run_bowtie2_map_se))

####################################################################################################

class mapReadToGGTansript(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="name of the project used for storing the analysis results."
											  "name should not contain any spaces or special characters (string ["
											  "=MyProject])")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	genomeName = luigi.Parameter()

	#sampleName = luigi.Parameter()

	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	domain = luigi.Parameter(description="domain of the organism . (string [=prokaryote] OR [=eukaryote])")

	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)

	
	
	def requires(self):

		return [mapReadToTans(projectName=self.projectName,
								readType=self.readType,
								domain = self.domain,
								genomeName=self.genomeName,
								aligner = self.aligner,
								adapter=self.adapter,
								sampleName=i)
				for i in [line.strip()
						  for line in
								open (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))]]


	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))

########################################################################################################
# Run Corsetanalysis
class clusterGGTranscripts(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	genomeName = luigi.Parameter()
	
	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")

	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)

	def requires(self):
		return [mapReadToTans(projectName=self.projectName,
								readType=self.readType,
								domain = self.domain,
								aligner=self.aligner,
								genomeName=self.genomeName,
								adapter=self.adapter,
								sampleName=i)
				for i in [line.strip()
						  for line in
								open (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))]]


	def output(self):

		corset_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_bowtie2_Corset_" + self.readType + "/")

		return {'out1': luigi.LocalTarget(corset_folder + "clusters.txt"),
				'out2': luigi.LocalTarget(corset_folder + "counts.txt")
				}


	def run(self):
		corset_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_bowtie2_Corset_" + self.readType + "/")
		
		gg_transcript_map_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_bowtie2_map_" + self.readType + "/")

		super_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_bowtie2_Corset_" + self.readType, "SuperTranscript/")
		gg_assembled_transcript_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "GGAssembly",
													 "trinity_GGAssembly_" + self.aligner +"_" + self.readType + "/")


		# Paired end reads

		input_group_file = ((os.path.join(GlobalParameter().basefolder, self.projectName, "group.tsv")))

		# Command to generate Rockhopper Paired-end input
		cmd_corset_read_input = utility.prepare_corset_input(input_group_file)


		cmd_run_corset = "[ -d  {corset_folder} ] || mkdir -p {corset_folder}; cd {corset_folder}; " \
							"corset -D 99999999999  {corset_read_input} {gg_transcript_map_folder}*.bam " \
							.format(corset_folder=corset_folder,
								 	corset_read_input=cmd_corset_read_input,
									gg_transcript_map_folder=gg_transcript_map_folder)

		cmd_run_supertrans = "python $(which Lace.py) " \
								"--cores 1 " \
								"{gg_assembled_transcript_folder}Trinity-GG.fasta " \
								"{corset_folder}clusters.txt " \
								"--outputDir {super_transcript_folder} " \
			.format(super_transcript_folder=super_transcript_folder,corset_folder=corset_folder,
					gg_assembled_transcript_folder=gg_assembled_transcript_folder)



		print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset)
		print (run_cmd(cmd_run_corset))

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_supertrans)
		print (run_cmd(cmd_run_supertrans))
########################################################################################################
# EDGER Analysis on corsetCounts
# Run EDGER
class denovoDEA(luigi.Task, TimeTask):
	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")
	assembler = luigi.Parameter(description="""name of the assembler to be used to assemble corrected RNASeq reads . (
											 string [=rockhopper2] OR [=spades] OR [=trinity].
											 NOTE: rockhopper2 is exclusively used to assemble bacterial RNASeq reads and must
											 not be used for domain eukaryote """)
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")

	deaMethod = luigi.Parameter(description='''name of the method to be used for differential expression analysis
																				=string([=deseq2] OR =[edger])''')
	resultTag = luigi.Parameter(default="cond2VScond1", description="""tag to be appended to output folder (string=[
										 cond2vscond1])""")
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	reportName = luigi.Parameter(default="Corset_DESeq2_HTML_Report")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")


	#batchCond = luigi.Parameter(default="",significant=False, description="blocking factor.  (string [=batch] OR [=NULL]")
	#coresetQuant = luigi.Parameter()

	factorInt = luigi.Parameter(default="group",description="factor of intrest column of the target file (string ["
															"=group]). ")
	refCond = luigi.Parameter(default="WT", description="reference biological condition.  (string [=WT]")
	targetFile = luigi.Parameter(default="target.txt", description="path to the design/target file. (string ["
																   "=target.txt]")
	alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
	pAdjustMethod = luigi.Parameter(default="BH", description="p-value adjustment method.  (string ["
																	  "=BH] OR [=BY]")

	fitType = luigi.Parameter(default="mean", description="mean-variance relationship.  (string ["
																	  "=parametric] OR [=local] OR [=mean]")

	sizeFactor = luigi.Parameter(default="median", description="method to estimate the size factors.  (string ["
															   "=median] OR =[shorth]")




	def requires(self):
		return [clusterContigs(domain=self.domain,
					   assembler=self.assembler,
					   projectName=self.projectName,
					   readType=self.readType,
					   adapter=self.adapter)]

	def output(self):

		#edgeResultFolter = GlobalParameter().basefolder + "/" + self.projectName + "_DEAnalysis/" + self.assembler + "_" +  self.sampleGroupFile + "_" + self.readType + "/" + "edgeR/"
		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis", self.deaMethod +
									"_" + self.assembler + "_" +  self.resultTag + "_" + self.readType + "/")

		if all([self.readType == "paired", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}


		if all([self.readType == "paired", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}


	def run(self):
		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis", self.deaMethod +
									"_" + self.assembler + "_" + self.resultTag + "_" + self.readType + "/")

		corset_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.assembler +
									 "_DenovoSalmonQuant_" + self.readType, "Corset/")

		basefol = GlobalParameter().basefolder + "/"
		rmd_DESeq2File = GlobalParameter().utilityFolder + "PlotDESEQ2.Rmd"
		rmd_edgeRFile = GlobalParameter().utilityFolder + "PlotEDGER.Rmd"


		cmd_run_corset_deseq = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
						  "cd {resultFolder};" \
						  "corset_DESeq2.r " \
						  "-t {basefol}{targetFile} " \
						  "-q {corset_folder}counts.txt " \
						   "-v {factorInt} " \
						   "-c {refCond} " \
						   "-f {fitType} " \
						   "-a {alpha} " \
						   "-p {pAdjustMethod} " \
						   "-l {sizeFactor} " \
							   "-T {rmd_DESeq2File}" \
						.format(resultFolder=resultFolder,
								basefol=basefol,
								targetFile=self.targetFile,
								corset_folder=corset_folder,
								factorInt=self.factorInt,
								refCond=self.refCond,
								fitType=self.fitType,
								alpha=self.alpha,
								pAdjustMethod=self.pAdjustMethod,
								sizeFactor=self.sizeFactor,
								rmd_DESeq2File=rmd_DESeq2File)

		cmd_run_corset_edger = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
							   "cd {resultFolder};" \
							   "corset_edgeR.r " \
							   "-t {basefol}{targetFile} " \
							   "-q {corset_folder}counts.txt " \
							   "-v {factorInt} " \
							   "-c {refCond} " \
							   "-f {fitType} " \
							   "-a {alpha} " \
							   "-p {pAdjustMethod} " \
							   "-l {sizeFactor} " \
							   "-T {rmd_edgeRFile}" \
			.format(resultFolder=resultFolder,
					basefol=basefol,
					targetFile=self.targetFile,
					corset_folder=corset_folder,
					factorInt=self.factorInt,
					refCond=self.refCond,
					fitType=self.fitType,
					alpha=self.alpha,
					pAdjustMethod=self.pAdjustMethod,
					sizeFactor=self.sizeFactor,
					rmd_edgeRFile=rmd_edgeRFile)



		if all([self.readType == "paired", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

##############################################################################################################
#Run EDGER
		if all([self.readType == "paired", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.readType == "single", self.domain == "prokaryote", self.assembler == "rockhopper2",
				self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "trinity",
				self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.readType == "paired", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.readType == "single", self.domain == "eukaryote", self.assembler == "spades",
				self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

########################################################################################################################
# EDGER Analysis on corsetCounts
# Run EDGER
class genomeGuidedDEA(luigi.Task, TimeTask):
	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")
	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")

	deaMethod = luigi.Parameter(description='''name of the method to be used for differential expression analysis
																				=string([=deseq2] OR =[edger])''')
	resultTag = luigi.Parameter(default="cond2VScond1", description="""tag to be appended to output folder (string=[
										 cond2vscond1])""")
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	reportName = luigi.Parameter(default="Corset_DESeq2_HTML_Report")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	genomeName = luigi.Parameter()
	#batchCond = luigi.Parameter(default="",significant=False, description="blocking factor.  (string [=batch] OR [=NULL]")
	#coresetQuant = luigi.Parameter()

	factorInt = luigi.Parameter(default="group",description="factor of intrest column of the target file (string ["
															"=group]). ")
	refCond = luigi.Parameter(default="WT", description="reference biological condition.  (string [=WT]")
	targetFile = luigi.Parameter(default="target.txt", description="path to the design/target file. (string ["
																   "=target.txt]")
	alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
	pAdjustMethod = luigi.Parameter(default="BH", description="p-value adjustment method.  (string ["
																	  "=BH] OR [=BY]")

	fitType = luigi.Parameter(default="mean", description="mean-variance relationship.  (string ["
																	  "=parametric] OR [=local] OR [=mean]")

	sizeFactor = luigi.Parameter(default="median", description="method to estimate the size factors.  (string ["
															   "=median] OR =[shorth]")




	def requires(self):
		return [clusterGGTranscripts(domain=self.domain,
					   aligner=self.aligner,
					   projectName=self.projectName,
					   readType=self.readType,
					   genomeName=self.genomeName,
					   adapter=self.adapter)]

	def output(self):

		#edgeResultFolter = GlobalParameter().basefolder + "/" + self.projectName + "_DEAnalysis/" + self.assembler + "_" +  self.sampleGroupFile + "_" + self.readType + "/" + "edgeR/"
		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis", self.deaMethod +
									"_" + "GG_Trinity" + "_" +  self.resultTag + "_" + self.readType + "/")

	

		if all([self.readType == "paired", self.domain == "eukaryote", self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.deaMethod == "deseq2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

				

		if all([self.readType == "paired", self.domain == "eukaryote", self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.deaMethod == "edger"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		

	def run(self):
		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis", self.deaMethod +
									"_" + "GG_Trinity" + "_" +  self.resultTag + "_" + self.readType + "/")

		corset_folder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
										   "GG_Transcript_bowtie2_Corset_" + self.readType + "/")

		basefol = GlobalParameter().basefolder + "/"

		rmd_DESeq2File = GlobalParameter().utilityFolder + "PlotDESEQ2.Rmd"
		rmd_edgeRFile = GlobalParameter().utilityFolder + "PlotEDGER.Rmd"


		cmd_run_corset_deseq = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
						  "cd {resultFolder};" \
						  "corset_DESeq2.r " \
						  "-t {basefol}{targetFile} " \
						  "-q {corset_folder}counts.txt " \
						   "-v {factorInt} " \
						   "-c {refCond} " \
						   "-f {fitType} " \
						   "-a {alpha} " \
						   "-p {pAdjustMethod} " \
						   "-l {sizeFactor} " \
							   "-T {rmd_DESeq2File}" \
						.format(resultFolder=resultFolder,
								basefol=basefol,
								targetFile=self.targetFile,
								corset_folder=corset_folder,
								factorInt=self.factorInt,
								refCond=self.refCond,
								fitType=self.fitType,
								alpha=self.alpha,
								pAdjustMethod=self.pAdjustMethod,
								sizeFactor=self.sizeFactor,
								rmd_DESeq2File=rmd_DESeq2File)

		cmd_run_corset_edger = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
							   "cd {resultFolder};" \
							   "corset_edgeR.r " \
							   "-t {basefol}{targetFile} " \
							   "-q {corset_folder}counts.txt " \
							   "-v {factorInt} " \
							   "-c {refCond} " \
							   "-f {fitType} " \
							   "-a {alpha} " \
							   "-p {pAdjustMethod} " \
							   "-l {sizeFactor} " \
							   "-T {rmd_edgeRFile}" \
			.format(resultFolder=resultFolder,
					basefol=basefol,
					targetFile=self.targetFile,
					corset_folder=corset_folder,
					factorInt=self.factorInt,
					refCond=self.refCond,
					fitType=self.fitType,
					alpha=self.alpha,
					pAdjustMethod=self.pAdjustMethod,
					sizeFactor=self.sizeFactor,
					rmd_edgeRFile=rmd_edgeRFile)



		if all([self.readType == "paired", self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))

		if all([self.readType == "single",self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
			print (run_cmd(cmd_run_corset_deseq))


##############################################################################################################
#Run EDGER
		if all([self.readType == "paired", self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))

		if all([self.readType == "single", self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
			print (run_cmd(cmd_run_corset_edger))


##########################################################################################################################
##Draft Genome based Differential Expression Analysis

# Run Prokka
class annotateGenome(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
												 name should not contain any spaces or special characters 
												 (string [=MyProject])""")

	genomeName = luigi.Parameter(description="""name of the draft or assembled genome file [with out .fna extension] 
											present in folder /raw_data/genome string(=[my_genome]""")


	transcriptName = luigi.Parameter(description="""name of the predicted transcript, which must be the genome name with
												out file extension. string[=my_genome]""")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	readType = luigi.Parameter(description="sequencing read type. (string [=single] OR [=paired]")

	n_cpu = luigi.Parameter(default="4", description="number of CPU to be used. (int [=8]")

	minContigLength = luigi.Parameter(default="200", description="minimum length of contigs. (int [=200]" )

	parameter = luigi.Parameter(default="--cpu 4")

	def requires(self):
		return []

	def output(self):

		genomeAnnotationFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome", self.transcriptName + "/")
		return {'out1': luigi.LocalTarget(genomeAnnotationFolder +"/"  + self.transcriptName + ".ffn"),
				'out2': luigi.LocalTarget(genomeAnnotationFolder + "/" + "tx2gene.csv")
				}

	def run(self):
		genomeAnnotationFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome", self.transcriptName + "/")
		draftGenomeFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "genome", self.genomeName + "/")


		cmd_run_prokka = "prokka {draftGenomeFolder}{genomeName}.fna " \
						 "--cpu {n_cpu} " \
						 "--prefix {genomeName} " \
						 "--mincontiglen {minContigLength} " \
						 "--outdir {genomeAnnotationFolder} --force --rfam" \
			.format(genomeName=self.genomeName, transcriptName=self.transcriptName,
					genomeAnnotationFolder=genomeAnnotationFolder,
					minContigLength=self.minContigLength,
					draftGenomeFolder=draftGenomeFolder,
					n_cpu=self.n_cpu,
					parameter=self.parameter)

		print ("****** NOW RUNNING COMMAND ******: " + cmd_run_prokka)
		print (run_cmd(cmd_run_prokka))

		print("Generating tx2gene")
		awk_cmd = 'BEGIN{FS="\\t"}{print ""$1"," ""$1"*" $7}'
		cmd_run_tx2gene = "cd {genomeAnnotationFolder}; " \
						  "awk '{awk}' {genomeName}.tsv > tx2gene.csv" \
			.format(genomeAnnotationFolder=genomeAnnotationFolder,
					genomeName=self.genomeName,
					awk=awk_cmd)

		print ("****** NOW RUNNING COMMAND ******: " + cmd_run_tx2gene)
		print (run_cmd(cmd_run_tx2gene))

		print("Generating GTF from GFF")
		sed_cmd2 = ''' sed 's/ID//g' '''
		awk_cmd2 = 'BEGIN{OFS="\\t"}{print $1,"PROKKA","CDS",$2,$3,".",$4,".","gene_id " $5}'
		cmd_run_gff2gtf = "cd {genomeAnnotationFolder}; " \
						  "grep -v '#' {genomeName}.gff | grep 'ID=' | {sed_cmd2} | " \
						  "cut -f1,4,5,7,9 | awk '{awk_cmd2}' > {genomeName}.gtf " \
			.format(genomeAnnotationFolder=genomeAnnotationFolder,
					genomeName=self.genomeName,
					sed_cmd2=sed_cmd2,
					awk_cmd2=awk_cmd2)
		print ("****** NOW RUNNING COMMAND ******: " + cmd_run_gff2gtf)
		print (run_cmd(cmd_run_gff2gtf))

###########################################################################################################################
# Run Salmon Index
class indexTranscript(luigi.Task):

	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
												name should not contain any spaces or special characters 
												(string [=MyProject])""")
	sampleName = luigi.Parameter()
	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in 
											folder /raw_data/genome/ string(=[my_genome]""")

	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")
	transcriptName = luigi.Parameter(description="""name of the predicted transcript, which must be the genome name with
											out file extension. string[=my_genome]""")

	predTranscript = luigi.Parameter(description="""for prokaryote, if it is required to predict transcript from genome
											  using PROKKA. (string [=yes] OR [=no])""")

	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired]""")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	quantMethod = luigi.Parameter(default="salmon", description="""method to be used for quantification of transcripts. string (=[salmon] OR =[kallisto])""" )
	Salmon_Index_Parameter = luigi.Parameter(default="--type quasi -k 31")
	Kallisto_Index_Parameter = luigi.Parameter(default="-k 31")

	def requires(self):
		if (self.predTranscript == "yes"):
			return [bbduk(projectName=self.projectName,
							   sampleName=self.sampleName,
							   readType=self.readType,
							   adapter=self.adapter),

					annotateGenome(projectName=self.projectName,
								   genomeName=self.genomeName,
								   transcriptName=self.transcriptName,
								   readType=self.readType,
								   adapter=self.adapter)]

		if (self.predTranscript == "no"):
			return [bbduk(projectName=self.projectName,
							   sampleName=self.sampleName,
							   readType=self.readType,
							   adapter=self.adapter)]


	def output(self):
		TranscriptIndexFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.quantMethod + "_index_" + self.readType + "/")

		if (self.quantMethod == "salmon") and (self.readType == "paired"):
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
					'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}

		if (self.quantMethod == "salmon") and (self.readType == "single"):
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
					'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}

		if (self.quantMethod == "kallisto") and (self.readType == "paired"):
			return {'out': luigi.LocalTarget(TranscriptIndexFolder + "/" + "kallisto.index")}

		if (self.quantMethod == "kallisto") and (self.readType == "single"):
			return {'out': luigi.LocalTarget(TranscriptIndexFolder + "/" + "kallisto.index")}

	def run(self):
		TranscriptIndexFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.quantMethod + "_index_" + self.readType + "/")
		transcriptFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome",self.transcriptName + "/")


		cmd_run_salmon_index_pe = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {transcriptFolder}{transcriptName}.ffn " \
								  "-i {TranscriptIndexFolder}" \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					transcriptName=self.transcriptName,
					transcriptFolder=transcriptFolder,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_salmon_index_se = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {transcriptFolder}{transcriptName}.ffn " \
								  "-i {TranscriptIndexFolder}" \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					transcriptName=self.transcriptName,
					transcriptFolder=transcriptFolder,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_kallisto_index_pe = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
									"cd {TranscriptIndexFolder}; " \
									"kallisto index " \
									"--index=kallisto.index {transcriptFolder}{transcriptName}.ffn " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					transcriptName=self.transcriptName,
					transcriptFolder=transcriptFolder,
					Kallisto_Index_Parameter=self.Kallisto_Index_Parameter)

		cmd_run_kallisto_index_se = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
									"cd {TranscriptIndexFolder}; " \
									"kallisto index " \
									"--index=kallisto.index {transcriptFolder}{transcriptName}.ffn " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					transcriptName=self.transcriptName,
					transcriptFolder=transcriptFolder,
					Kallisto_Index_Parameter=self.Kallisto_Index_Parameter)

		if (self.quantMethod == "salmon") and (self.readType == "paired"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_pe)
			print (run_cmd(cmd_run_salmon_index_pe))


		if (self.quantMethod == "salmon") and (self.readType == "single"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_se)
			print (run_cmd(cmd_run_salmon_index_se))



		if (self.quantMethod == "kallisto") and (self.readType == "paired"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_index_pe)
			print (run_cmd(cmd_run_kallisto_index_pe))


		if (self.quantMethod == "kallisto") and (self.readType == "single"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_index_se)
			print (run_cmd(cmd_run_kallisto_index_se))


###########################################################################################################################
class transQuant(luigi.Task, TimeTask):

	projectName=luigi.Parameter(description="""name of the project used for storing the analysis results.
								 name should not contain any spaces or special characters (string [=MyProject])""")

	sampleName = luigi.Parameter(description="""name of a single sample.""")


	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired]""")

	quantMethod = luigi.Parameter(description="""method to be used for quantification of transcripts. string (=[salmon] OR =[kallisto])""")

	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in 
											folder /raw_data/genome/ string(=[my_genome]""")

	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")

	transcriptName = luigi.Parameter(description="""name of the predicted transcript, which must be the genome name with
											out file extension. string[=my_genome]""")

	predTranscript=luigi.Parameter(description="""for prokaryote, if it is required to predict transcript from genomeName
											  using PROKKA. (string [=yes] OR [=no])""")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	quantParameter = luigi.Parameter(default="-p 2")

	def requires(self):

		return [indexTranscript(projectName=self.projectName,
								readType=self.readType,
								quantMethod=self.quantMethod,
								predTranscript=self.predTranscript,
								domain = self.domain,
								transcriptName = self.transcriptName,
								genomeName=self.genomeName,
								adapter=self.adapter,
								sampleName=i)
				for i in [line.strip()
						  for line in
								open (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))]]

	def output(self):

		TranscriptQuantFolder = os.path.join(GlobalParameter().basefolder, self.projectName,"TransQuant",
											 self.quantMethod + "_quant_" + self.readType, self.quantMethod + "_" + "quant" + "/")

		if (self.readType == "paired") and (self.quantMethod == "salmon"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "cmd_info.json")
					}
		if (self.readType == "single") and (self.quantMethod == "salmon"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "cmd_info.json")
					}

		if (self.readType == "paired") and (self.quantMethod == "kallisto"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "abundance.tsv"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "run_info.json")
					}

		if (self.readType == "single") and (self.quantMethod == "kallisto"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "abundance.tsv"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "run_info.json")
					}

	def run(self):

		transcriptFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome",self.transcriptName + "/")

		TranscriptIndexFolder = os.path.join(GlobalParameter().basefolder,
								 self.projectName,
								     "TransQuant",
								 self.quantMethod + "_index_" + self.readType + "/")

		TranscriptQuantFolder = os.path.join(GlobalParameter().basefolder,
								 self.projectName,
								     "TransQuant",
						                 self.quantMethod + "_quant_" + self.readType,
								 self.quantMethod + "_" + "quant" + "/")

		TranscriptQuantSampleFolder = os.path.join(GlobalParameter().basefolder,
								       self.projectName,
									   "TransQuant",
					   self.quantMethod + "_quant_" + self.readType,
			                               self.quantMethod + "_" + "quant",
								  self.sampleName + "/")

		TranscriptMapFolder = os.path.join(GlobalParameter().basefolder,
										   self.projectName,
										   "TransQuant",
										   self.quantMethod + "_map_" + self.readType + "/")

		cleanedReadFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "CleanedReads", "Cleaned_" +
										 self.readType + "-end_Reads/")
		#genome_folder = GlobalParameter().basefolder + "/" + self.projectName + "_annotated_" + self.genome + "/"
		'''genomeAnnotationFolder = os.path.join(GlobalParameter().basefolder,
											  self.projectName,
											  "annotated_" + self.genome + "/")'''

		cmd_run_salmon_quant_pe = "[ -d {TranscriptQuantFolder} ] || mkdir -p {TranscriptQuantFolder}; " \
								  "[ -d {TranscriptMapFolder} ] || mkdir -p {TranscriptMapFolder}; " \
								  "cd {TranscriptMapFolder}; " \
								  "salmon quant --no-version-check -p 12 " \
								  "-i {TranscriptIndexFolder} " \
								  "-l A " \
								  "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								  "-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								  "--dumpEq " \
								  "--output {TranscriptQuantFolder}{sampleName} " \
								  "--validateMappings " \
								  "--writeMappings | samtools view -bS - | samtools sort -m 10G " \
								  "-o {TranscriptMapFolder}{sampleName}.bam" \
			.format(TranscriptQuantFolder=TranscriptQuantFolder,
					TranscriptMapFolder=TranscriptMapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					TranscriptIndexFolder=TranscriptIndexFolder,
					quantParameter=self.quantParameter)

		cmd_run_salmon_quant_se = "[ -d {TranscriptQuantFolder} ] || mkdir -p {TranscriptQuantFolder}; " \
								  "mkdir -p {TranscriptMapFolder}; " \
								  "cd {TranscriptQuantFolder};" \
								  "salmon quant --no-version-check -p 12 " \
								  "--seqBias " \
								  "--gcBias " \
								  "-i {TranscriptIndexFolder} " \
								  "-l A " \
								  "-r {cleanedReadFolder}{sampleName}.fastq " \
								  "--dumpEq " \
								  "--output {TranscriptQuantFolder}{sampleName} " \
								  "--validateMappings " \
								  "--writeMappings | samtools view -bS - | samtools sort -m 10G " \
								  "-o {TranscriptMapFolder}{sampleName}.bam" \
			.format(TranscriptQuantFolder=TranscriptQuantFolder,
					TranscriptMapFolder=TranscriptMapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					TranscriptIndexFolder=TranscriptIndexFolder,
					quantParameter=self.quantParameter)


		cmd_run_kallisto_quant_pe = "[ -d {TranscriptQuantFolder} ] || mkdir -p {TranscriptQuantFolder}; " \
									"mkdir -p {TranscriptMapFolder}; mkdir -p {TranscriptQuantFolder}; " \
									"cd {TranscriptQuantFolder}; " \
									"kallisto quant " \
									"--threads=16 " \
									"--index={TranscriptIndexFolder}kallisto.index " \
									"--output-dir={TranscriptQuantFolder}{sampleName} " \
									"{cleanedReadFolder}{sampleName}_R1.fastq " \
									"{cleanedReadFolder}{sampleName}_R2.fastq " \
									"--pseudobam " \
									"--genomebam " \
									"-g {transcriptFolder}{transcriptName}.gtf" \
			.format(TranscriptQuantFolder=TranscriptQuantFolder,
					TranscriptIndexFolder=TranscriptIndexFolder,
					transcriptFolder=transcriptFolder,
					TranscriptMapFolder=TranscriptMapFolder,
					sampleName=self.sampleName,
					transcriptName=self.transcriptName,
					cleanedReadFolder=cleanedReadFolder,
					quantParameter=self.quantParameter)

		mean_sd = ''' awk 'BEGIN { t=0.0;sq=0.0; n=0;} ;NR%4==2 {n++;L=length($0);t+=L;sq+=L*L;}END{
				m=t/n;printf("-l %f -s %f",m,(sq/n-m*m)+0.000001);}' '''

		cmd_run_mean_sd_se = "cd {cleanedReadFolder}; {mean_sd} {cleanedReadFolder}{sampleName}.fastq > {cleanedReadFolder}{sampleName}.txt" \
			.format(cleanedReadFolder=cleanedReadFolder,
					mean_sd=mean_sd,
					sampleName=self.sampleName)

		cmd_run_kallisto_quant_se = "[ -d {TranscriptQuantFolder} ] || mkdir -p {TranscriptQuantFolder}; " \
									"mkdir -p {TranscriptMapFolder}; mkdir -p {TranscriptQuantSampleFolder}; " \
									"cd {TranscriptQuantSampleFolder}; " \
									"kallisto quant " \
									"--threads=12 " \
									"--index={TranscriptIndexFolder}kallisto.index " \
									"--output-dir={TranscriptQuantFolder}{sampleName} " \
									"--single " \
									"{cleanedReadFolder}{sampleName}.fastq " \
									"$(<{cleanedReadFolder}{sampleName}.txt) " \
									"--pseudobam " \
									"--genomebam " \
									"-g {transcriptFolder}{transcriptName}.gtf" \
			.format(TranscriptQuantFolder=TranscriptQuantFolder,
					TranscriptQuantSampleFolder=TranscriptQuantSampleFolder,
					TranscriptIndexFolder=TranscriptIndexFolder,
					transcriptFolder=transcriptFolder,
					TranscriptMapFolder=TranscriptMapFolder,
					sampleName=self.sampleName,
					transcriptName=self.transcriptName,
					cleanedReadFolder=cleanedReadFolder,
					quantParameter=self.quantParameter)

		cmd_move_kallisto_bam_se = "mv {TranscriptQuantSampleFolder}pseudoalignments.bam {TranscriptMapFolder}{sampleName}.bam " \
								   "&& " \
								   "mv {TranscriptQuantSampleFolder}pseudoalignments.bam.bai {TranscriptMapFolder}{sampleName}.bam.bai" \
			.format(TranscriptQuantSampleFolder=TranscriptQuantSampleFolder,
					sampleName=self.sampleName,
					TranscriptMapFolder=TranscriptMapFolder)

		cmd_move_kallisto_bam_pe = "mv {TranscriptQuantSampleFolder}pseudoalignments.bam {TranscriptMapFolder}{sampleName}.bam " \
								   "&& " \
								   "mv {TranscriptQuantSampleFolder}pseudoalignments.bam.bai {TranscriptMapFolder}{sampleName}.bam.bai" \
			.format(TranscriptQuantSampleFolder=TranscriptQuantSampleFolder,
					sampleName=self.sampleName,
					TranscriptMapFolder=TranscriptMapFolder)

		if (self.readType == "paired") and (self.quantMethod == "salmon"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))


		if (self.readType == "single") and (self.quantMethod == "salmon"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))


		if (self.readType == "paired") and (self.quantMethod == "kallisto"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_quant_pe)
			print (run_cmd(cmd_run_kallisto_quant_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_move_kallisto_bam_pe)
			print (run_cmd(cmd_move_kallisto_bam_pe))

		if (self.readType == "single") and (self.quantMethod == "kallisto"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_mean_sd_se)
			print (run_cmd(cmd_run_mean_sd_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_quant_se)
			print (run_cmd(cmd_run_kallisto_quant_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_move_kallisto_bam_se)
			print (run_cmd(cmd_move_kallisto_bam_se))


##########################################################################################################################
class quantifyTranscripts(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
							 name should not contain any spaces or special characters (string [=MyProject])""")


	quantMethod = luigi.Parameter(default="salmon", description="""method to be used for quantification of 
																transcripts. string (=[salmon] OR =[kallisto])""")

	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")

	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in 
											folder /raw_data/genome/ string(=[my_genome]""")

	transcriptName=luigi.Parameter(description="""name of the predicted transcript, which must be the genome name with
											out file extension. string[=my_genome]""")

	predTranscript=luigi.Parameter(description="""for prokaryote, if it is required to predict transcript from genomeName
											  using PROKKA. (string [=yes] OR [=no])""")

	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired]""")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	def requires(self):
		return [transQuant(projectName=self.projectName,
							   readType=self.readType,
							   quantMethod=self.quantMethod,
							   transcriptName=self.transcriptName,
							   genomeName=self.genomeName,
							   domain=self.domain,
							   predTranscript=self.predTranscript,
							   adapter=self.adapter,
							   sampleName=i)
				for i in [line.strip()
						  for line in
								open (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))]]

	'''def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))'''
	def output(self):
		tx2genefolder = os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome",self.transcriptName + "/")
		return {'out': luigi.LocalTarget(tx2genefolder +"/"  + "tx2gene.csv")}
	
	def run(self):
		transcriptFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome",self.transcriptName + "/")
		cmd_tx2gene_from_gtf = "[ -d  {transcriptFolder} ] || mkdir -p {transcriptFolder}; " \
						 "cd {transcriptFolder}; tx2gene.R " \
						 "-a gtf " \
						 "-p {transcriptFolder}{transcriptName}.gtf " \
						 "-o tx2gene.csv".format(transcriptFolder=self.transcriptFolder,
									   transcriptName=self.transcriptName)
		
		print("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_from_gtf)
		print (run_cmd(cmd_tx2gene_from_gtf))


#################################################################################################################################
class transcriptomeBasedDEA(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
												 name should not contain any spaces or special characters 
												 (string [=MyProject])""")

	domain = luigi.Parameter(description="""domain of the organism. 
											(string [=prokaryote] OR [=eukaryote])""")

	predTranscript=luigi.Parameter(description="""for prokaryote, if it is required to predict transcript from the draft
												  OR assembled genome using PROKKA.            
												  (string [=yes] OR [=no])""")

	readType = luigi.Parameter(description="""sequencing read type. 
											  (string [=single] OR [=paired]""")

	quantMethod = luigi.Parameter(description="""method to be used for quantification of transcripts. 
												  string (=[salmon] OR =[kallisto])""")

	deaMethod = luigi.Parameter(description="""method to be used for differential expression analysis. 
												  string (=[deseq2] OR =[edger])""")

	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in 
											folder /raw_data/genome/ string(=[my_genome]""")

	transcriptName=luigi.Parameter(description="""name of the predicted transcript, which must be the genome name with
											out file extension. string[=my_genome]""")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")



	resultTag = luigi.Parameter()
	reportName = luigi.Parameter(default="featureCount_DESeq2_HTML_Report")

	factorInt = luigi.Parameter(default="group", description="factor of intrest column of the target file (string ["
															 "=group]). ")
	refCond = luigi.Parameter(default="WT", description="reference biological condition.  (string [=WT]")
	targetFile = luigi.Parameter(default="target.txt", description="path to the design/target file. (string ["
																   "=target.txt]")
	alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
	pAdjustMethod = luigi.Parameter(default="BH", description="p-value adjustment method.  (string ["
															  "=BH] OR [=BY]")

	fitType = luigi.Parameter(default="mean", description="mean-variance relationship.  (string ["
														  "=parametric] OR [=local] OR [=mean]")

	sizeFactor = luigi.Parameter(default="median", description="method to estimate the size factors.  (string ["
															   "=median] OR =[shorth]")



	def requires(self):
		return [quantifyTranscripts(projectName=self.projectName,
						   readType=self.readType,
						   quantMethod=self.quantMethod,
						   predTranscript=self.predTranscript,
						   genomeName=self.genomeName,
						   domain=self.domain,
						   transcriptName=self.transcriptName,
						   adapter=self.adapter)]


	def output(self):

		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis",
									self.deaMethod + "_" + self.quantMethod
									+ "_" + self.resultTag + "_" + self.readType + "/")

		if all([self.readType == "paired", self.quantMethod == "salmon", self.deaMethod == "deseq2"]):
			return {'out': luigi.LocalTarget(resultFolder + "/"  + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.quantMethod == "salmon", self.deaMethod == "deseq2"]):
			return {'out': luigi.LocalTarget(resultFolder + "/" + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.quantMethod == "kallisto", self.deaMethod == "deseq2"]):
			return {'out': luigi.LocalTarget(resultFolder + "/" + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.quantMethod == "kallisto", self.deaMethod == "deseq2"]):
			return {'out': luigi.LocalTarget(resultFolder + "/" + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.quantMethod == "salmon", self.deaMethod == "edger"]):
			return {'out': luigi.LocalTarget(resultFolder + "/" + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.quantMethod == "salmon", self.deaMethod == "edger"]):
			return {'out': luigi.LocalTarget(resultFolder + "/" + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.quantMethod == "kallisto", self.deaMethod == "edger"]):
			return {'out': luigi.LocalTarget(resultFolder + "/" + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.quantMethod == "kallisto", self.deaMethod == "edger"]):
			return {'out': luigi.LocalTarget(resultFolder + "/" + self.reportName + "/" + "index.html")}


	def run(self):
		baseFolder = GlobalParameter().basefolder
		tx2genefolder = os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome",self.transcriptName + "/")

		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis", self.deaMethod +
									"_" + self.quantMethod + "_" + self.resultTag + "_" + self.readType + "/")

		transcriptQuantFolder = os.path.join(GlobalParameter().basefolder,self.projectName,"TransQuant",
									 self.quantMethod + "_quant_" + self.readType + "/")

		rmd_DESeq2File = GlobalParameter().utilityFolder + "PlotDESEQ2.Rmd"
		rmd_edgeRFile = GlobalParameter().utilityFolder + "PlotEDGER.Rmd"

		cmd_run_salmon_DESeq2 = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
									  "cd {resultFolder};" \
									  "salmon_DESeq2.r " \
									  "-t {baseFolder}{targetFile} " \
									  "-q {transcriptQuantFolder} " \
									  "-G {tx2genefolder} " \
									  "-v {factorInt} " \
									  "-c {refCond} " \
									  "-f {fitType} " \
									  "-a {alpha} " \
									  "-p {pAdjustMethod} " \
									  "-l {sizeFactor} " \
									  "-T {rmd_DESeq2File}" \
			.format(resultFolder=resultFolder,
					baseFolder=baseFolder,
					targetFile=self.targetFile,
					transcriptQuantFolder=transcriptQuantFolder,
					tx2genefolder=tx2genefolder,
					factorInt=self.factorInt,
					refCond=self.refCond,
					fitType=self.fitType,
					alpha=self.alpha,
					pAdjustMethod=self.pAdjustMethod,
					sizeFactor=self.sizeFactor,
					rmd_DESeq2File=rmd_DESeq2File)

		cmd_run_kallisto_DESeq2 = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
								"cd {resultFolder};" \
								"kallisto_DESeq2.r " \
								"-t {baseFolder}{targetFile} " \
								"-q {transcriptQuantFolder} " \
								"-G {tx2genefolder} " \
								"-v {factorInt} " \
								"-c {refCond} " \
								"-f {fitType} " \
								"-a {alpha} " \
								"-p {pAdjustMethod} " \
								"-l {sizeFactor} " \
								"-T {rmd_DESeq2File}" \
			.format(resultFolder=resultFolder,
					baseFolder=baseFolder,
					targetFile=self.targetFile,
					transcriptQuantFolder=transcriptQuantFolder,
					tx2genefolder=tx2genefolder,
					factorInt=self.factorInt,
					refCond=self.refCond,
					fitType=self.fitType,
					alpha=self.alpha,
					pAdjustMethod=self.pAdjustMethod,
					sizeFactor=self.sizeFactor,
					rmd_DESeq2File=rmd_DESeq2File)

		cmd_run_salmon_edgeR = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
								"cd {resultFolder};" \
								"salmon_edgeR.r " \
								"-t {baseFolder}{targetFile} " \
								"-q {transcriptQuantFolder} " \
								"-G {tx2genefolder} " \
								"-v {factorInt} " \
								"-c {refCond} " \
								"-a {alpha} " \
								"-p {pAdjustMethod} " \
								"-T {rmd_edgeRFile}" \
			.format(resultFolder=resultFolder,
					baseFolder=baseFolder,
					targetFile=self.targetFile,
					transcriptQuantFolder=transcriptQuantFolder,
					tx2genefolder=tx2genefolder,
					factorInt=self.factorInt,
					refCond=self.refCond,
					fitType=self.fitType,
					alpha=self.alpha,
					pAdjustMethod=self.pAdjustMethod,
					sizeFactor=self.sizeFactor,
					rmd_edgeRFile=rmd_edgeRFile)

		cmd_run_kallisto_edgeR = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
								  "cd {resultFolder};" \
								  "kallisto_edgeR.r " \
								  "-t {baseFolder}{targetFile} " \
								  "-q {transcriptQuantFolder} " \
								  "-G {tx2genefolder} " \
								  "-v {factorInt} " \
								  "-c {refCond} " \
								  "-a {alpha} " \
								  "-p {pAdjustMethod} " \
								  "-T {rmd_edgeRFile}" \
			.format(resultFolder=resultFolder,
					baseFolder=baseFolder,
					targetFile=self.targetFile,
					transcriptQuantFolder=transcriptQuantFolder,
					tx2genefolder=tx2genefolder,
					factorInt=self.factorInt,
					refCond=self.refCond,
					fitType=self.fitType,
					alpha=self.alpha,
					pAdjustMethod=self.pAdjustMethod,
					sizeFactor=self.sizeFactor,
					rmd_edgeRFile=rmd_edgeRFile)



		##########Commands
		#DESEQ2
		if all([self.readType == "paired", self.quantMethod == "salmon", self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_DESeq2)
			print (run_cmd(cmd_run_salmon_DESeq2))

		if all([self.readType == "single", self.quantMethod == "salmon", self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_DESeq2)
			print (run_cmd(cmd_run_salmon_DESeq2))

		if all([self.readType == "paired", self.quantMethod == "kallisto", self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_DESeq2)
			print (run_cmd(cmd_run_kallisto_DESeq2))

		if all([self.readType == "single", self.quantMethod == "kallisto", self.deaMethod == "deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_DESeq2)
			print (run_cmd(cmd_run_kallisto_DESeq2))

		#EDGER

		if all([self.readType == "paired", self.quantMethod == "salmon", self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_edgeR)
			print (run_cmd(cmd_run_salmon_edgeR))

		if all([self.readType == "single", self.quantMethod == "salmon", self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_edgeR)
			print (run_cmd(cmd_run_salmon_edgeR))

		if all([self.readType == "paired", self.quantMethod == "kallisto", self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_edgeR)
			print (run_cmd(cmd_run_kallisto_edgeR))

		if all([self.readType == "single", self.quantMethod == "kallisto", self.deaMethod == "edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_edgeR)
			print (run_cmd(cmd_run_kallisto_edgeR))








##############################################################################################
#Reference based expression Analysis
# ###################################################################################################################################################################################################
# Run STAR Aligner  Followed by HTSeq Count and DESeq2 #Aligner Gives CoreDump Error
# Run STAR Index
class indexGenome(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")

	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in 
											folder /raw_data/genome/ string(=[my_genome])""")

	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")

	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)

	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")

	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")



	def requires(self):

		return []

	def output(self):
		genomeIndexFolder = os.path.join(GlobalParameter().basefolder, "raw_data","genome",
								   self.genomeName + "_genome_index", self.aligner + "_index" + "/")

		if all([self.readType == "paired", self.aligner == "star", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}

		if all([self.readType == "single", self.aligner == "star", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}

		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}

		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}


		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genomeName + ".1.ht2")}
		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genomeName + ".1.ht2")}
		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genomeName + ".1.ht2")}
		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genomeName + ".1.ht2")}




		if all([((self.readType == "paired") or (self.readType == "single")),
				((self.domain == "prokaryote") or (self.domain == "eukaryote")),
				(self.aligner == "segemehl")]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "index.idx")}


		if all([((self.readType == "paired") or (self.readType == "single")),
				((self.domain == "prokaryote") or (self.domain == "eukaryote")),
				(self.aligner == "dart")]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "index.bwt")}


		if all([((self.readType == "paired") or (self.readType == "single")),
				((self.domain == "prokaryote") or (self.domain == "eukaryote")),
				(self.aligner == "bowtie2")]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "index.1.bt2")}




	def run(self):
		refGenomeFolder = os.path.join(GlobalParameter().basefolder,"raw_data","genome", self.genomeName + "/")

		genomeIndexFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "genome",
										 self.genomeName + "_genome_index", self.aligner + "_index" + "/")

		genomeFastaFile= "{refGenomeFolder}{genomeName}.fna".format(refGenomeFolder=refGenomeFolder,
																	genomeName=self.genomeName)

		gsan = utility.genomeSAindexNbases(genomeFastaFile)

		cmd_run_star_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder}; STAR " \
								"--runMode genomeGenerate " \
								"--genomeSAindexNbases {gsan} " \
								"--genomeFastaFiles {refGenomeFolder}{genomeName}.fna " \
								"--genomeDir {genomeIndexFolder} " \
			.format(genomeIndexFolder=genomeIndexFolder,
					genomeName=self.genomeName,
					gsan=gsan,
					refGenomeFolder=refGenomeFolder)



		cmd_run_bowtie2_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
								   "bowtie2-build -f {refGenomeFolder}{genomeName}.fna {genomeIndexFolder}index " \
			.format(genomeIndexFolder=genomeIndexFolder,
					genomeName=self.genomeName,
					refGenomeFolder=refGenomeFolder)



		cmd_run_dart_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
								"bwt_index {refGenomeFolder}{genomeName}.fna {genomeIndexFolder}index " \
			.format(genomeIndexFolder=genomeIndexFolder,
					genomeName=self.genomeName,
					refGenomeFolder=refGenomeFolder)


		cmd_run_hisat2_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
								"hisat2-build {refGenomeFolder}{genomeName}.fna {genomeIndexFolder}{genomeName} " \
			.format(genomeIndexFolder=genomeIndexFolder,
					genomeName=self.genomeName,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
							   "segemehl.x -x {genomeIndexFolder}index.idx -d {refGenomeFolder}{genomeName}.fna  " \
			.format(genomeIndexFolder=genomeIndexFolder,
					genomeName=self.genomeName,
					refGenomeFolder=refGenomeFolder)

		if  self.aligner == "star":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_star_index)
			print (run_cmd(cmd_run_star_index))


		if  self.aligner == "dart":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_index)
			print (run_cmd(cmd_run_dart_index))


		if  self.aligner == "bowtie2":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_index)
			print (run_cmd(cmd_run_bowtie2_index))

		if self.aligner == "hisat2":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_index)
			print (run_cmd(cmd_run_hisat2_index))

		if self.aligner == "segemehl":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_index)
			print (run_cmd(cmd_run_segemehl_index))



#################################################################################################################################################
class alignReads(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	readType = luigi.Parameter()
	sampleName = luigi.Parameter()
	aligner = luigi.Parameter()
	genomeName = luigi.Parameter()
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	domain = luigi.Parameter()

	def requires(self):


		return [indexGenome(projectName=self.projectName,
							readType=self.readType,
							aligner=self.aligner,
							genomeName=self.genomeName,
							domain=self.domain,
							adapter=self.adapter),

				bbduk(projectName=self.projectName,
					  sampleName=self.sampleName,
					  readType=self.readType,
					  adapter=self.adapter)]

	def output(self):

		#mapFolder = os.path.join(GlobalParameter().basefolder,self.projectName,"TransQuant", self.domain + "_"
								 #+ self.aligner + "_map_" + self.readType, self.sampleName + "/")

		mapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.aligner + "_map_" + self.readType + "/")

		qualimapFolder = os.path.join(GlobalParameter().basefolder,self.projectName,"QCReports", self.aligner +
									  "_qualimap_" + self.readType, self.sampleName + "/")


		#STAR OUTPUT
		if all([self.readType == "paired", self.aligner == "star",self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "star",self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}


		#HISAT2 OUTPUT
		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}


		#BOWTIE2 OUTPUT
		if all([self.readType == "paired", self.aligner == "bowtie2", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "bowtie22", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}


		#DART OUTPUT
		if all([self.readType == "paired", self.aligner == "dart", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "dart", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "dart", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}



		# SEGEMEHL OUTPUT
		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}




	def run(self):

		refGenomeFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "genome", self.genomeName + "/")

		cleanedReadFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "CleanedReads",
										 "Cleaned_" + self.readType + "-end_Reads/")

		genomeIndexFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "genome",
										 self.genomeName + "_genome_index", self.aligner + "_index" + "/")

		mapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.aligner + "_map_" + self.readType + "/")

		#mapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.domain + "_"
								 #+ self.aligner + "_map_" + self.readType, self.sampleName + "/")

		qualimapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "QCReports", self.domain + "_" +
									  self.aligner + "_qualimap_" + self.readType, self.sampleName + "/")

		##########################################################################################################
		# 1 SEGEMEHL Aligner                                                                                          #
		##########################################################################################################

		cmd_run_segemehl_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "segemehl.x -t 2 " \
								   "-i {genomeIndexFolder}index.idx " \
									"-d {refGenomeFolder}{genomeName}.fna " \
								   "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
								   "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
								   "-S " \
								   "|samtools view -bS - | samtools sort " \
								   "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					genomeName=self.genomeName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t 2 " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genomeName}.fna " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					genomeName=self.genomeName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t 2 " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genomeName}.fna " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					genomeName=self.genomeName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t 2 " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genomeName}.fna " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					genomeName=self.genomeName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)




		################################################################################################################
		# 1. HISAT2 Aligner
		################################################################################################################
		cmd_run_hisat2_pe_prok = "[-d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder};hisat2 --dta -x {genomeIndexFolder}{genomeName} " \
								"-p 2 " \
								 "--max-intronlen 20 " \
								 "--no-spliced-alignment " \
								"-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								"-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								"| samtools view -bS - | samtools sort " \
								"-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					sampleName=self.sampleName,
					genomeName=self.genomeName)

		cmd_run_hisat2_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genomeName} " \
								"-p 2 " \
								"--max-intronlen 20 " \
								"--no-spliced-alignment " \
								"-U {cleanedReadFolder}{sampleName}.fastq " \
								"| samtools view -bS - | samtools sort " \
								"-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genomeName=self.genomeName,
					sampleName=self.sampleName)


		cmd_run_hisat2_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genomeName} " \
								 "-p 2 " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								 "| samtools view -bS - | samtools sort " \
								 "-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genomeName=self.genomeName,
					sampleName=self.sampleName)

		cmd_run_hisat2_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genomeName} " \
								 "-p 2 " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "| samtools view -bS - | samtools sort -m 5G " \
								 "-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genomeName=self.genomeName,
					sampleName=self.sampleName)



		################################################################################################################
		# 2. STAR Aligner
		################################################################################################################
		cmd_run_star_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "STAR --runMode alignReads " \
							  "--alignIntronMax 1 " \
							  "--outSAMtype BAM " \
							  "SortedByCoordinate " \
							  "--limitBAMsortRAM 16000000000 " \
							  "--alignSJDBoverhangMin 999 " \
							  "--runThreadN 2 " \
							  "--genomeDir {genomeIndexFolder} " \
							  "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							  "{cleanedReadFolder}{sampleName}_R2.fastq " \
							  "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_star_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							   "cd {mapFolder}; " \
							   "STAR --runMode alignReads " \
							   "--alignIntronMax 1 " \
							   "--outSAMtype BAM " \
							   "SortedByCoordinate " \
							   "--limitBAMsortRAM 16000000000 " \
							   "--alignSJDBoverhangMin 999 " \
							   "--runThreadN 2 " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)



		cmd_run_star_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							   "cd {mapFolder}; " \
							   "STAR --runMode alignReads " \
							   "--outSAMstrandField intronMotif " \
							   "--outSAMtype BAM " \
							   "SortedByCoordinate " \
							   "--limitBAMsortRAM 16000000000 " \
							   "--runThreadN 2 " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							   "{cleanedReadFolder}{sampleName}_R2.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_star_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							   "cd {mapFolder}; " \
							   "STAR --runMode alignReads " \
							   "--outSAMstrandField intronMotif " \
							   "--outSAMtype BAM " \
							   "SortedByCoordinate " \
							   "--limitBAMsortRAM 16000000000 " \
							   "--runThreadN 2 " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)


		cmd_star_bam_rename = "mv {mapFolder}{sampleName}_Aligned.sortedByCoord.out.bam " \
							  "{mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder, sampleName=self.sampleName)


		##########################################################################################################
		#3 DART Aligner                                                                                          #
		##########################################################################################################

		cmd_run_dart_map_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "dart -intron 1 -t 2 " \
							  "-i {genomeIndexFolder}index " \
							  "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
							  "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
							  "-j {mapFolder}{sampleName}_junctions.tab " \
							  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)


		cmd_run_dart_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "dart -intron 1 -t 2 " \
							  "-i {genomeIndexFolder}index " \
							  "-f {cleanedReadFolder}{sampleName}.fastq " \
							  "-j {mapFolder}{sampleName}_junctions.tab " \
							  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_dart_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "dart -t 2 " \
								   "-i {genomeIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
								   "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_dart_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "dart -t 2 " \
								   "-i {genomeIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		#Sort Unsorted BAM file and rempove the unsorted bam file

		cmd_run_dart_sort_bam = "cd {mapFolder}; samtools sort -@8 {sampleName}_unsorted.bam -o {sampleName}.bam".format(
			mapFolder=mapFolder, sampleName=self.sampleName)

		cmd_run_dart_remove_unsorted_bam = "cd {mapFolder}; rm {sampleName}_unsorted.bam".format(
			mapFolder=mapFolder, sampleName=self.sampleName)

		################################################################################################################
		#4 Bowtie2 Aligner
		################################################################################################################

		cmd_run_bowtie2_map_pe = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "bowtie2 -x {genomeIndexFolder}index " \
								 "-p 2 " \
								 "-a " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq |samtools view -bS - | samtools " \
								 "sort " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)


		cmd_run_bowtie2_map_se = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "bowtie2 -x {genomeIndexFolder}index " \
								 "-p 2 " \
								 "-a " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "|samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)



#######################################################################################################################
#Mapping Quality Assessment using Qualimap
########################################################################################################################
		cmd_run_qualimap =  "[ -d {qualimapFolder} ] || mkdir -p {qualimapFolder}; " \
							"cd {qualimapFolder}; " \
							"qualimap bamqc " \
							"-bam {mapFolder}{sampleName}.bam " \
							"--java-mem-size=8G " \
							"-outdir {qualimapFolder} " \
							"-outfile {sampleName}_QualiMap " \
							"-outformat PDF:HTML" \
			.format(qualimapFolder=qualimapFolder,
					mapFolder=mapFolder,
					sampleName=self.sampleName)

########################################################################################################################
# Call aligner commands
########################################################################################################################
#Run Bowtie2: Only for prokaryotes
		if (self.readType == "paired") and (self.aligner == "bowtie2") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_pe)
			print (run_cmd(cmd_run_bowtie2_map_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if (self.readType == "single") and (self.aligner == "bowtie2") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_se)
			print (run_cmd(cmd_run_bowtie2_map_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

# Run Segmehl
		if (self.readType == "paired") and (self.aligner == "segemehl") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_euk)
			print (run_cmd(cmd_run_segemehl_map_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.readType == "single") and (self.aligner == "segemehl") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_euk)
			print (run_cmd(cmd_run_segemehl_map_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.readType == "paired") and (self.aligner == "segemehl") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_prok)
			print (run_cmd(cmd_run_segemehl_map_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.readType == "single") and (self.aligner == "segemehl") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_prok)
			print (run_cmd(cmd_run_segemehl_map_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


#Run DART

		if (self.readType == "paired") and (self.aligner == "dart") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_euk)
			print (run_cmd(cmd_run_dart_map_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print (run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print (run_cmd(cmd_run_dart_remove_unsorted_bam))


		if (self.readType == "single") and (self.aligner == "dart") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_se_euk)
			print (run_cmd(cmd_run_dart_map_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print (run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print (run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.readType == "paired") and (self.aligner == "dart") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_prok)
			print (run_cmd(cmd_run_dart_map_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print (run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print (run_cmd(cmd_run_dart_remove_unsorted_bam))


		if (self.readType == "single") and (self.aligner == "dart") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_se_prok)
			print (run_cmd(cmd_run_dart_map_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print (run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print (run_cmd(cmd_run_dart_remove_unsorted_bam))


#Run HISAT2
####
		if (self.readType == "paired") and (self.aligner == "hisat2") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_euk)
			print (run_cmd(cmd_run_hisat2_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if (self.readType == "single") and (self.aligner == "hisat2") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_euk)
			print (run_cmd(cmd_run_hisat2_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.readType == "paired") and (self.aligner == "hisat2") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_prok)
			print (run_cmd(cmd_run_hisat2_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if (self.readType == "single") and (self.aligner == "hisat2") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_prok)
			print (run_cmd(cmd_run_hisat2_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

#Run STAR

		if (self.readType == "paired") and (self.aligner == "star") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_prok)
			print (run_cmd(cmd_run_star_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if (self.readType == "single") and (self.aligner == "star") and (self.domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_prok)
			print (run_cmd(cmd_run_star_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.readType == "paired") and (self.aligner == "star") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_euk)
			print (run_cmd(cmd_run_star_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.readType == "single") and (self.aligner == "star") and (self.domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_euk)
			print (run_cmd(cmd_run_star_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


#####################################################################################################
class mapReads(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	#sampleName = luigi.Parameter()
	readType = luigi.Parameter()
	aligner = luigi.Parameter()
	genomeName = luigi.Parameter()
	domain = luigi.Parameter()
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")

	def requires(self):

		return [alignReads(readType=self.readType,
						   aligner=self.aligner,
						   genomeName=self.genomeName,
						   projectName=self.projectName,
						   domain = self.domain,
						   adapter=self.adapter,
						   sampleName=i)
				for i in [line.strip()
						  for line in
								open (os.path.join(GlobalParameter().basefolder, self.projectName, "samples.txt"))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))


########################################################################################################################
# FeatureCounts to generate count file
########################################################################################################################
class featureCounts(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
											 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
											 for domain eukaryote """)
	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in 
											folder /raw_data/genome/ string(=[my_genome])""")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	domain = luigi.Parameter(description='''domain of the organism. 
											(string [=prokaryote] OR [=eukaryote])''')

	featureType = luigi.Parameter(default="exon",description='''Features used for read counting will be 
					  extracted from annotation using the provided value. Feature is defined by the 3rd column of 
					  your GTF file. To get all the possible features you can use 'cut -f 3 <GTF_file> | sort | uniq'
					  and select the appropriate feature. string(=[exon]) ''')

	attributeType = luigi.Parameter(default="gene_id",description='''Specify attribute type in GTF annotation. 
											 string(=[gene_id])''')

	strandType = luigi.Parameter(default="0",description='''perform strand-specific read counting. int([=0]unstranded) 
															OR [=1] stranded] OR [=2] reversely-stranded. default[
															=0]''')
	threads = luigi.Parameter(default="2",description='''number of threads. (int [=4])''')

	def requires(self):

		return [mapReads(projectName=self.projectName,
						 readType=self.readType,
						 aligner=self.aligner,
						 genomeName=self.genomeName,
						 domain=self.domain,
						 adapter=self.adapter)]

	def output(self):
		QuantFolder = os.path.join(GlobalParameter().basefolder, self.projectName ,"TransQuant" , self.aligner + "_quant_" + self.readType + "/")


		if (self.readType == "paired") and (self.aligner == "star") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "star") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "paired") and (self.aligner == "star") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "star") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}



		if (self.readType == "paired") and (self.aligner == "dart") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "dart") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "paired") and (self.aligner == "dart") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "dart") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}




		if (self.readType == "paired") and (self.aligner == "hisat2") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "hisat2") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "paired") and (self.aligner == "hisat2") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "hisat2") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}




		if (self.readType == "paired") and (self.aligner == "segemehl") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "segemehl") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}
		if (self.readType == "paired") and (self.aligner == "segemehl") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "segemehl") and (self.domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}



		if (self.readType == "paired") and (self.aligner == "bowtie2") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.readType == "single") and (self.aligner == "bowtie2") and (self.domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}




	def run(self):

		refGenomeFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "genome", self.genomeName + "/")


		genomeIndexFolder = os.path.join(GlobalParameter().basefolder, "raw_data", "genome",
										 self.genomeName + "_genome_index", self.aligner + "_index" + "/")

		mapFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant", self.aligner + "_map_" + self.readType + "/")

		QuantFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
								   self.aligner + "_quant_" + self.readType + "/")



		cmd_run_featureCount_pe = "[ -d {QuantFolder} ] || mkdir -p {QuantFolder}; " \
										"cd {QuantFolder}; " \
										"featureCounts -a {refGenomeFolder}{genomeName}.gtf " \
										"-t {featureType} " \
										"-g {attributeType} " \
										"-s {strandType} " \
										"-T {threads} " \
										"-p " \
										"-o counts.txt " \
										"{mapFolder}*.bam " \
			.format(QuantFolder=QuantFolder,
					refGenomeFolder=refGenomeFolder,
					featureType=self.featureType,
					strandType=self.strandType,
					attributeType=self.attributeType,
					threads=self.threads,
					genomeIndexFolder=genomeIndexFolder,
					mapFolder=mapFolder,
					genomeName=self.genomeName)

		cmd_run_featureCount_se = "[ -d {QuantFolder} ] || mkdir -p {QuantFolder}; " \
									  "cd {QuantFolder}; " \
									  "featureCounts -a {refGenomeFolder}{genomeName}.gtf " \
									  "-t {featureType} " \
									  "-g {attributeType} " \
									  "-s {strandType} " \
									  "-T {threads} " \
									  "-o counts.txt " \
									  "{mapFolder}*.bam " \
			.format(QuantFolder=QuantFolder,
					refGenomeFolder=refGenomeFolder,
					featureType=self.featureType,
					attributeType=self.attributeType,
					strandType=self.strandType,
					threads=self.threads,
					genomeIndexFolder=genomeIndexFolder,
					mapFolder=mapFolder,
					genomeName=self.genomeName)



		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.readType == "single", self.aligner == "star", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "star", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))



		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))



		if all([self.readType == "single", self.aligner == "dart", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.readType == "single", self.aligner == "dart", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))


		if all([self.readType == "single", self.aligner == "bowtie2", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.readType == "paired", self.aligner == "bowtie2", self.domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))



##################################################################
class generateCounts(luigi.Task, TimeTask):

	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
									name should not contain any spaces or special characters (string [=MyProject])""")
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
												 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2] OR = [segemehl])
												 NOTE: star and segemehl demands high memory. bowtie2 should not be used 
												 for domain eukaryote """)
	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in 
												folder /raw_data/genome/ string(=[my_genome])""")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	domain = luigi.Parameter(description='''domain of the organism. 
												(string [=prokaryote] OR [=eukaryote])''')

	featureType = luigi.Parameter(default="exon", description='''Features used for read counting will be 
						  extracted from annotation using the provided value. Feature is defined by the 3rd column of 
						  your GTF file. To get all the possible features you can use 'cut -f 3 <GTF_file> | sort | uniq'
						  and select the appropriate feature. string(=[exon]) ''')

	attributeType = luigi.Parameter(default="gene_id", description='''Specify attribute type in GTF annotation. 
												 string(=[gene_id])''')

	strandType = luigi.Parameter(default="0", description='''perform strand-specific read counting. int([=0]unstranded) 
																OR [=1] stranded] OR [=2] reversely-stranded. default[
																=0]''')
	threads = luigi.Parameter(default="2",description='''number of threads. (int [=4])''')





	def requires(self):
		return [featureCounts(readType=self.readType,
								aligner=self.aligner,
								genomeName=self.genomeName,
								projectName=self.projectName,
								domain=self.domain,
								adapter=self.adapter,
								featureType=self.featureType,
								attributeType=self.attributeType,
								strandType=self.strandType,
								threads=self.threads)]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))


##################################################################
# Differential Gene Expression Analysis
#########################################################################################################################
class genomeBasedDEA(luigi.Task, TimeTask):
	projectName = luigi.Parameter(description="""name of the project used for storing the analysis results.
								name should not contain any spaces or special characters (string [=MyProject])""")
	domain = luigi.Parameter(description="""domain of the organism . (string [=prokaryote] OR [=eukaryote])""")
	readType = luigi.Parameter(description="""sequencing read type. (string [=single] OR [=paired])""")
	aligner = luigi.Parameter(description="""name of the aligner to be used to map clean reads to indexed genome . (
											 string [=star] OR [=hisat2] OR [=dart] OR [=bowtie2])
											 NOTE: star demands high memory. bowtie2 should not be used
											 for domain eukaryote """)
	genomeName = luigi.Parameter(description="""name of the genome file [with out .fna extension] present in
											folder /raw_data/genome/ string(=[my_genome])""")
	adapter = luigi.Parameter(default="./utility/adapters.fasta.gz")
	deaMethod = luigi.Parameter(description="""Differential Expression Analysis Method to be used (string [
											   =deseq2] OR =[edger])""")
	resultTag = luigi.Parameter()
	reportName = luigi.Parameter(default="DEA_Report")


	factorInt = luigi.Parameter(default="group", description="factor of intrest column of the target file (string ["
															 "=group]). ")
	refCond = luigi.Parameter(default="WT", description="reference biological condition.  (string [=WT]")
	targetFile = luigi.Parameter(default="target.txt", description="path to the design/target file. (string ["
																   "=target.txt]")
	alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
	pAdjustMethod = luigi.Parameter(default="BH", description="p-value adjustment method.  (string ["
															  "=BH] OR [=BY]")

	fitType = luigi.Parameter(default="mean", description="mean-variance relationship.  (string ["
														  "=parametric] OR [=local] OR [=mean]")

	sizeFactor = luigi.Parameter(default="median", description="method to estimate the size factors.  (string ["
															   "=median] OR =[shorth]")

	def requires(self):

		return [generateCounts(projectName=self.projectName,
						   domain=self.domain,
						   readType=self.readType,
						   aligner=self.aligner,
						   genomeName=self.genomeName,
						   adapter=self.adapter)]

	def output(self):

		#edgeResultFolter = GlobalParameter().basefolder + "/" + self.projectName + "_DEAnalysis/" + self.assembler + "_" +  self.sampleGroupFile + "_" + self.readType + "/" + "edgeR/"
		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis", self.deaMethod +
									"_" + self.aligner  + "_" +  self.resultTag + "_" + self.readType + "/")
		#STAR Aligner Output
		if all([self.readType == "paired", self.domain == "prokaryote", self.aligner == "star"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "prokaryote", self.aligner == "star"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.aligner == "star"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.aligner == "star"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}


		if all([self.readType == "paired", self.domain == "prokaryote", self.aligner == "hisat2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "prokaryote", self.aligner == "hisat2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.aligner == "hisat2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.aligner == "hisat2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}



		if all([self.readType == "paired", self.domain == "prokaryote", self.aligner == "dart"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "prokaryote", self.aligner == "dart"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.aligner == "dart"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.aligner == "dart"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}



		if all([self.readType == "paired", self.domain == "prokaryote", self.aligner == "segemehl"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "prokaryote", self.aligner == "segemehl"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "paired", self.domain == "eukaryote", self.aligner == "segemehl"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "eukaryote", self.aligner == "segemehl"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}



		if all([self.readType == "paired", self.domain == "prokaryote", self.aligner == "bowtie2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}

		if all([self.readType == "single", self.domain == "prokaryote", self.aligner == "bowtie2"]):
			return {'out1': luigi.LocalTarget(resultFolder + self.reportName + "/" + "index.html")}






	def run(self):
		resultFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "DEAnalysis", self.deaMethod + "_" +
		self.aligner + "_" + self.resultTag + "_" + self.readType + "/")

		QuantFolder = os.path.join(GlobalParameter().basefolder, self.projectName, "TransQuant",
								   self.aligner + "_quant_" + self.readType + "/")

		basefol = GlobalParameter().basefolder + "/"

		rmd_DESeq2File = GlobalParameter().utilityFolder + "PlotDESEQ2.Rmd"
		rmd_edgeRFile = GlobalParameter().utilityFolder + "PlotEDGER.Rmd"


		cmd_run_featuteCount_DESeq2 = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
										"cd {resultFolder};" \
										"featureCount_DESeq2.r " \
										"-t {basefol}{targetFile} " \
										"-q {QuantFolder}counts.txt " \
										"-v {factorInt} " \
										"-c {refCond} " \
										"-f {fitType} " \
										"-a {alpha} " \
										"-p {pAdjustMethod} " \
										"-l {sizeFactor} " \
										"-T {rmd_DESeq2File}" \
											.format(resultFolder=resultFolder,
													basefol=basefol,
													targetFile=self.targetFile,
													QuantFolder=QuantFolder,
													factorInt=self.factorInt,
													refCond=self.refCond,
													fitType=self.fitType,
													alpha=self.alpha,
													pAdjustMethod=self.pAdjustMethod,
													sizeFactor=self.sizeFactor,
													rmd_DESeq2File=rmd_DESeq2File)

		cmd_run_featuteCount_edgeR = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
									  "cd {resultFolder};" \
									  "featureCount_edgeR.r " \
									  "-t {basefol}{targetFile} " \
									  "-q {QuantFolder}counts.txt " \
									  "-v {factorInt} " \
									  "-c {refCond} " \
									  "-f {fitType} " \
									  "-a {alpha} " \
									  "-p {pAdjustMethod} " \
									  "-l {sizeFactor} " \
									  "-T {rmd_edgeRFile}" \
			.format(resultFolder=resultFolder,
					basefol=basefol,
					targetFile=self.targetFile,
					QuantFolder=QuantFolder,
					factorInt=self.factorInt,
					refCond=self.refCond,
					fitType=self.fitType,
					alpha=self.alpha,
					pAdjustMethod=self.pAdjustMethod,
					sizeFactor=self.sizeFactor,
					rmd_edgeRFile=rmd_edgeRFile)

		if all([self.readType == "paired", self.aligner == "star", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))


		if all([self.readType == "single", self.aligner == "star", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))



		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))



		if all([self.readType == "paired", self.aligner == "dart", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "dart", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "dart", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))


		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "eukaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))


		if all([self.readType == "paired", self.aligner == "bowtie2", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.readType == "single", self.aligner == "bowtie2", self.domain == "prokaryote", self.deaMethod=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

####################################################################################################################
#Run EDGER
######################################################################################################################
		if all([self.readType == "paired", self.aligner == "star", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))


		if all([self.readType == "single", self.aligner == "star", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "paired", self.aligner == "star", self.domain == "eukaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "star", self.domain == "eukaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))



		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "paired", self.aligner == "hisat2", self.domain == "eukaryote",  self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "hisat2", self.domain == "eukaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))



		if all([self.readType == "paired", self.aligner == "dart", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "dart", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "paired", self.aligner == "dart", self.domain == "eukaryote",  self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "dart", self.domain == "eukaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))


		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "paired", self.aligner == "segemehl", self.domain == "eukaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "segemehl", self.domain == "eukaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))


		if all([self.readType == "paired", self.aligner == "bowtie2", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.readType == "single", self.aligner == "bowtie2", self.domain == "prokaryote", self.deaMethod=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))


############################################################################################################################
if __name__ == '__main__':
	luigi.run()
