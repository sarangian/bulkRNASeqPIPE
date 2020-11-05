import luigi
import os
import time
import subprocess
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.reFormatReads import reformat
from tasks.rnaSeq.alignmentFree.indexTransctriptome import indexTranscript

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
class GlobalParameter(luigi.Config):

	project_name=luigi.Parameter()
	rnaseq_dir=luigi.Parameter()
	genome_suffix=luigi.Parameter()
	genome_name=luigi.Parameter()
	genome_dir=luigi.Parameter()
	domain=luigi.Parameter()

	read_library_type=luigi.Parameter()
	
	transcriptome_dir=luigi.Parameter()
	transcriptome_name=luigi.Parameter()
	transcriptome_suffix=luigi.Parameter()

	annotation_dir=luigi.Parameter()
	annotation_suffix=luigi.Parameter()
	annotation_name=luigi.Parameter()

	adapter=luigi.Parameter()
	threads=luigi.Parameter()
	maxMemory=luigi.Parameter()
	
	feature_type=luigi.Parameter()
	adapter=luigi.Parameter()

class transQuant(luigi.Task):

	project_name=GlobalParameter().project_name
	genome_dir=GlobalParameter().genome_dir
	genome_name=GlobalParameter().genome_name
	genome_suffix=GlobalParameter().genome_suffix

	organism_domain=GlobalParameter().domain

	transcriptome_dir=GlobalParameter().transcriptome_dir
	transcriptome_name=GlobalParameter().transcriptome_name
	transcriptome_suffix=GlobalParameter().transcriptome_suffix

	rnaseq_dir=GlobalParameter().rnaseq_dir
	read_library_type=GlobalParameter().read_library_type

	annotation_dir=GlobalParameter().annotation_dir
	annotation_name=GlobalParameter().annotation_name
	annotation_file_type=GlobalParameter().annotation_suffix


	threads = GlobalParameter().threads
	adapter = GlobalParameter().adapter

	
	#local parameters
	sampleName = luigi.Parameter()
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	quant_method = luigi.ChoiceParameter(choices=["salmon", "kallisto"], var_type=str)
	

	def requires(self):
		if self.read_library_type=="pe":
			return [indexTranscript(quant_method=self.quant_method,
					pre_process_reads=self.pre_process_reads,
					sampleName=i)
				for i in [line.strip()
						for line in
						open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]
		if self.read_library_type=="se":
			return [indexTranscript(quant_method=self.quant_method,
					pre_process_reads=self.pre_process_reads,
					sampleName=i)
				for i in [line.strip()
						for line in
						open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]

	def output(self):

		TranscriptQuantFolder = os.path.join(os.getcwd(),self.project_name,"AlignmentFreeDEA","ReadQuant",self.quant_method+"_quant_"+self.read_library_type,
					self.quant_method + "_" + "quant" + "/")

		if (self.read_library_type == "pe") and (self.quant_method == "salmon"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "cmd_info.json")
					}
		if (self.read_library_type == "se") and (self.quant_method == "salmon"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "cmd_info.json")
					}

		if (self.read_library_type == "pe") and (self.quant_method == "kallisto"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "abundance.tsv"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "run_info.json")
					}

		if (self.read_library_type == "se") and (self.quant_method == "kallisto"):
			return {'out1': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "abundance.tsv"),
					'out2': luigi.LocalTarget(TranscriptQuantFolder + self.sampleName + "/" + "run_info.json")
					}

	def run(self):

		transcriptome_dir = os.path.join(os.getcwd(),GlobalParameter().transcriptome_dir + "/")
		TranscriptIndexFolder = os.path.join(os.getcwd(), self.project_name, "AlignmentFreeDEA","transcript_index", self.transcriptome_name + "_" + self.quant_method +"_index" + "/")
		TranscriptMapFolder = os.path.join(os.getcwd(), self.project_name,"AlignmentFreeDEA","transcript_map", self.transcriptome_name + "_"+self.quant_method + "_map_" + "/")


		TranscriptQuantFolder = os.path.join(os.getcwd(),self.project_name,"AlignmentFreeDEA","ReadQuant",self.quant_method+"_quant_"+self.read_library_type,
					self.quant_method + "_" + "quant" + "/")

		TranscriptQuantSampleFolder = os.path.join(os.getcwd(),
					self.project_name,"AlignmentFreeDEA",
					"ReadQuant",
					 self.quant_method + "_quant_" + self.read_library_type,
					 self.quant_method + "_" + "quant",
					 self.sampleName + "/")

		
		if self.pre_process_reads == "yes" and self.read_library_type=="pe":
			cleanedReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Cleaned_PE_Reads" + "/")
		if self.pre_process_reads=="yes" and self.read_library_type=="se":
			cleanedReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Cleaned_SE_Reads" + "/")

		if self.pre_process_reads == "no" and self.read_library_type=="pe":
			cleanedReadFolder = os.path.join(os.getcwd(),self.project_name,"ReadQC","Verified_PE_Reads" + "/")
		if self.pre_process_reads=="no" and self.read_library_type=="se":
			cleanedReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Verified_SE_Reads" + "/")


		cmd_run_salmon_quant_pe = "[ -d {TranscriptQuantFolder} ] || mkdir -p {TranscriptQuantFolder}; " \
						  "[ -d {TranscriptMapFolder} ] || mkdir -p {TranscriptMapFolder}; " \
								  "cd {TranscriptMapFolder}; " \
								  "salmon quant --no-version-check " \
								  "-p {threads} " \
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
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					TranscriptIndexFolder=TranscriptIndexFolder)

		cmd_run_salmon_quant_se = "[ -d {TranscriptQuantFolder} ] || mkdir -p {TranscriptQuantFolder}; " \
								  "mkdir -p {TranscriptMapFolder}; " \
								  "cd {TranscriptQuantFolder};" \
								  "salmon quant --no-version-check " \
								  "-p {threads} " \
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
					threads=self.threads,
					cleanedReadFolder=cleanedReadFolder,
					TranscriptIndexFolder=TranscriptIndexFolder)


		cmd_run_kallisto_quant_pe = "[ -d {TranscriptQuantFolder} ] || mkdir -p {TranscriptQuantFolder}; " \
									"mkdir -p {TranscriptMapFolder}; mkdir -p {TranscriptQuantFolder}; " \
									"cd {TranscriptQuantFolder}; " \
									"kallisto quant " \
									"--threads={threads} " \
									"--index={TranscriptIndexFolder}kallisto.idx " \
									"--output-dir={TranscriptQuantFolder}{sampleName} " \
									"{cleanedReadFolder}{sampleName}_R1.fastq " \
									"{cleanedReadFolder}{sampleName}_R2.fastq " \
									"--pseudobam " \
									"--genomebam " \
									"-g {annotation_dir}{annotation_name}.gtf" \
			.format(TranscriptQuantFolder=TranscriptQuantFolder,
					TranscriptIndexFolder=TranscriptIndexFolder,
					annotation_dir=GlobalParameter().annotation_dir,
					annotation_name=GlobalParameter().annotation_name,
					threads=self.threads,
					TranscriptMapFolder=TranscriptMapFolder,
					sampleName=self.sampleName,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder)

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
									"--threads={threads} " \
									"--index={TranscriptIndexFolder}kallisto.idx " \
									"--output-dir={TranscriptQuantFolder}{sampleName} " \
									"--single " \
									"{cleanedReadFolder}{sampleName}.fastq " \
									"$(<{cleanedReadFolder}{sampleName}.txt) " \
									"--pseudobam " \
									"--genomebam " \
									"-g {annotation_dir}{annotation_name}.gtf" \
			.format(TranscriptQuantFolder=TranscriptQuantFolder,
					TranscriptQuantSampleFolder=TranscriptQuantSampleFolder,
					TranscriptIndexFolder=TranscriptIndexFolder,
					annotation_dir=GlobalParameter().annotation_dir,
					annotation_name=GlobalParameter().annotation_name,
					TranscriptMapFolder=TranscriptMapFolder,
					sampleName=self.sampleName,
					threads=self.threads,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder)

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

		'''cmd_tx2gene_from_gtf = "cd {transcriptFolder}; tx2gene.R " \
							   "-a gtf " \
							   "-p {transcriptFolder}{transcriptName}.gtf " \
							   "-o tx2gene.csv" \
							   .format(transcriptFolder=os.path.join(GlobalParameter().basefolder, "raw_data", "transcriptome",self.transcriptName + "/"),
									   transcriptName=self.transcriptName)'''

				
		if (self.read_library_type == "pe") and (self.quant_method == "salmon"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))
			#print("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_from_gtf)
			#print (run_cmd(cmd_tx2gene_from_gtf))


		if (self.read_library_type == "se") and (self.quant_method == "salmon"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))
			#print("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_from_gtf)
			#print (run_cmd(cmd_tx2gene_from_gtf))


		if (self.read_library_type == "pe") and (self.quant_method == "kallisto"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_quant_pe)
			print (run_cmd(cmd_run_kallisto_quant_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_move_kallisto_bam_pe)
			print (run_cmd(cmd_move_kallisto_bam_pe))
			#print("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_from_gtf)
			#print (run_cmd(cmd_tx2gene_from_gtf))

		if (self.read_library_type == "se") and (self.quant_method == "kallisto"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_mean_sd_se)
			print (run_cmd(cmd_run_mean_sd_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_quant_se)
			print (run_cmd(cmd_run_kallisto_quant_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_move_kallisto_bam_se)
			print (run_cmd(cmd_move_kallisto_bam_se))
			#print("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_from_gtf)
			#print (run_cmd(cmd_tx2gene_from_gtf))


##########################################################################################################################
class alignmentFreeQuant(luigi.Task):

	project_name=GlobalParameter().project_name
	organism_domain = GlobalParameter().domain
	genome_dir=GlobalParameter().genome_dir
	genome_name=GlobalParameter().genome_name
	genome_suffix=GlobalParameter().genome_suffix

	transcriptome_dir=GlobalParameter().transcriptome_dir
	transcriptome_name=GlobalParameter().transcriptome_name
	transcriptome_suffix=GlobalParameter().transcriptome_suffix

	rnaseq_dir=GlobalParameter().rnaseq_dir
	read_library_type=GlobalParameter().read_library_type

	annotation_dir=GlobalParameter().annotation_dir
	annotation_name=GlobalParameter().annotation_name
	#annotation_file_type=GlobalParameter().annotation_suffix


	threads = GlobalParameter().threads
	adapter = GlobalParameter().adapter

	
	#local parameters
	#sampleName = luigi.Parameter()
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	quant_method = luigi.ChoiceParameter(choices=["salmon", "kallisto"], var_type=str)


	def requires(self):
		if self.read_library_type=="pe":
			return [transQuant(quant_method=self.quant_method,
						   pre_process_reads=self.pre_process_reads,
						   sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]

		if self.read_library_type=="se":
			return [transQuant(quant_method=self.quant_method,
						   pre_process_reads=self.pre_process_reads,
						   sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "congig", "se_samples.lst")))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(), "task_logs", 'task.generate.transcript.count.complete.{t}'.format(
			t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('Generate Transcript Count finished at {t}'.format(t=timestamp))
