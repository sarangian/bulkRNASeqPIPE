import luigi
import time
import os
import subprocess

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

class readqc(luigi.Task):

	read_library_type=GlobalParameter().read_library_type
	rnaseq_dir=GlobalParameter().rnaseq_dir
	read_suffix=GlobalParameter().read_suffix
	adapter=GlobalParameter().adapter


	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory
	projectName = GlobalParameter().project_name
	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
	



	def output(self):
		readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PreQC_" +self.read_library_type+"_reads" "/")
		
		if self.read_library_type == "pe":
			return {'out1': luigi.LocalTarget(readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(readQC_folder + self.sampleName + "_R2_fastqc.html")}

		if self.read_library_type == "se":
			return {'out1': luigi.LocalTarget(readQC_folder + self.sampleName + "_fastqc.html")}

		
	def run(self):
		readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PreQC_" + self.read_library_type+"_reads" "/")
		
		
		read_QC_log_folder = os.path.join(os.getcwd(), "log", self.read_library_type +"_ReadQC" + "/")



		cmd_raw_pe_qc = "[ -d  {readQC_folder} ] || mkdir -p {readQC_folder}; mkdir -p {read_QC_log_folder}; " \
					   "/usr/bin/time -v fastqc " \
						"-t {cpu} " \
						"{rnaseq_dir}{sampleName}_R1.{read_suffix} " \
						"{rnaseq_dir}{sampleName}_R2.{read_suffix} " \
						"-o {readQC_folder} 2>&1 | tee  {read_QC_log_folder}{sampleName}_pe_fastqc.log".format(
													   sampleName=self.sampleName,
													   read_suffix=GlobalParameter().read_suffix,
													   readQC_folder=readQC_folder,
													   cpu=GlobalParameter().threads,
													   rnaseq_dir=GlobalParameter().rnaseq_dir,
													   read_QC_log_folder=read_QC_log_folder)

		
		cmd_raw_se_qc = "[ -d  {readQC_folder} ] || mkdir -p {readQC_folder};   mkdir -p {read_QC_log_folder}; " \
						"fastqc " \
						"--threads {cpu} " \
						"{rnaseq_dir}{sampleName}.{read_suffix} " \
						"-o {readQC_folder} " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_se_fastqc.log".format(
													   sampleName=self.sampleName,
													   read_suffix=GlobalParameter().read_suffix,
													   readQC_folder=readQC_folder,
													   cpu=GlobalParameter().threads,
													   read_QC_log_folder=read_QC_log_folder,
													   rnaseq_dir=GlobalParameter().rnaseq_dir)

		
		if self.read_library_type == "se":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_se_qc)
			print (run_cmd(cmd_raw_se_qc))

		if self.read_library_type == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pe_qc)
			print (run_cmd(cmd_raw_pe_qc))

		
class rawReadsQC(luigi.Task):
	read_library_type = GlobalParameter().read_library_type

	def requires(self):

		if self.read_library_type == "pe":
			return [readqc(
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]

		
		if self.read_library_type == "se":
			return [readqc(
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]


		

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.read.qc.analysis.complete.{t}'.format(
			t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('Read QC Assessment finished at {t}'.format(t=timestamp))










		





