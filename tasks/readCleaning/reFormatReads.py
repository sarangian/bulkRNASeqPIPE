import luigi
import time
import os
import subprocess
from tasks.readCleaning.cleanedReadQC import *

class GlobalParameter(luigi.Config):

	project_name=luigi.Parameter()
	rnaseq_dir=luigi.Parameter()
	read_suffix=luigi.Parameter()
	adapter=luigi.Parameter()
	threads=luigi.Parameter()
	maxMemory=luigi.Parameter()
	read_library_type=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	adapter=luigi.Parameter()
	

def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output

def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)

createFolder("task_logs")

class reformat(luigi.Task):
	project_name=GlobalParameter().project_name
	read_library_type=GlobalParameter().read_library_type
	rnaseq_dir=GlobalParameter().rnaseq_dir
	read_suffix=GlobalParameter().read_suffix

	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory
	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")


	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory
	
	projectName = luigi.Parameter(default="VerifiedReads")
	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")


	def output(self):
				
		if self.read_library_type == "pe":
			verified_read_folder = os.path.join(os.getcwd(), self.project_name, "ReadQC","Verified_PE_Reads" + "/")
			return {'out1': luigi.LocalTarget(verified_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(verified_read_folder + self.sampleName + "_R2.fastq")}

		if self.read_library_type == "se":
			verified_read_folder = os.path.join(os.getcwd(), self.project_name, "Verified_SE_Reads" + "/")
			return {'out1': luigi.LocalTarget(verified_read_folder + self.sampleName + ".fastq")}

		
		
	def run(self):
		if self.read_library_type == "pe":
			verified_read_folder = os.path.join(os.getcwd(), self.project_name, "ReadQC","Verified_PE_Reads" + "/")
		if self.read_library_type == "se":
			verified_read_folder = os.path.join(os.getcwd(), self.project_name, "ReadQC","Verified_SE_Reads" + "/")

				
		read_verification_log_folder = os.path.join(os.getcwd(), "log","VerifiedReads" + "/")
		
		cmd_clean_pe = "[ -d  {verified_read_folder} ] || mkdir -p {verified_read_folder}; mkdir -p {verified_read_folder}; " \
					   "[ -d  {read_verification_log_folder} ] || mkdir -p {read_verification_log_folder}; " \
					   "reformat.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "tossbrokenreads=t " \
					   "verifypaired=t " \
					   "in1={rnaseq_dir}{sampleName}_R1.{read_suffix} " \
					   "in2={rnaseq_dir}{sampleName}_R2.{read_suffix} " \
					   "out={verified_read_folder}{sampleName}_R1.fastq " \
					   "out2={verified_read_folder}{sampleName}_R2.fastq " \
					   " 2>&1 | tee {read_verification_log_folder}{sampleName}_pe_reformat_run.log "\
			.format(Xmx=GlobalParameter().maxMemory,
					cpu=GlobalParameter().threads,
					rnaseq_dir=os.path.join(GlobalParameter().rnaseq_dir),
					read_suffix=GlobalParameter().read_suffix,
					sampleName=self.sampleName,
					verified_read_folder=verified_read_folder,
					read_verification_log_folder=read_verification_log_folder)

			
		##################
		cmd_clean_se = "[ -d  {verified_read_folder} ] || mkdir -p {verified_read_folder}; " \
					   "[ -d  {read_verification_log_folder} ] || mkdir -p {read_verification_log_folder}; " \
					   "reformat.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "tossbrokenreads=t " \
					   "in1={rnaseq_dir}{sampleName}.{read_suffix} " \
					   "out={verified_read_folder}{sampleName}.fastq " \
					   " 2>&1 | tee {read_verification_log_folder}{sampleName}_se_reformat_run.log " \
			.format(Xmx=GlobalParameter().maxMemory,
					cpu=GlobalParameter().threads,
					rnaseq_dir=GlobalParameter().rnaseq_dir,
					read_suffix=GlobalParameter().read_suffix,
					sampleName=self.sampleName,
					verified_read_folder=verified_read_folder,
					read_verification_log_folder=read_verification_log_folder)

	
		if self.read_library_type == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pe)
			print(run_cmd(cmd_clean_pe))

		if self.read_library_type == "se":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_se)
			print(run_cmd(cmd_clean_se))


class reformatReads(luigi.Task):
	
	read_library_type = GlobalParameter().read_library_type

	def requires(self):

		if self.read_library_type == "pe":
			return [

					[reformat(sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]
			       ]

		
		if self.read_library_type == "se":
			return [[reformat(sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]
					]

	
	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.validate.reads.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('read validation finished at {t}'.format(t=timestamp))