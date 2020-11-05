import luigi
import os
import time
import subprocess
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.reFormatReads import reformat
from tasks.rnaSeq.annotation.make_tx_to_gene import makeTx2Gene

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

class GlobalParameter(luigi.Config):

	project_name=luigi.Parameter()
	genome_suffix=luigi.Parameter()
	genome_name=luigi.Parameter()
	genome_dir=luigi.Parameter()
	domain=luigi.Parameter()
	
	transcriptome_dir=luigi.Parameter()
	transcriptome_name=luigi.Parameter()
	transcriptome_suffix=luigi.Parameter()

	annotation_dir=luigi.Parameter()
	annotation_suffix=luigi.Parameter()
	annotation_name=luigi.Parameter()

	adapter=luigi.Parameter()
	threads=luigi.Parameter()
	maxMemory=luigi.Parameter()


	read_library_type=luigi.Parameter()

	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	feature_type=luigi.Parameter()
	adapter=luigi.Parameter()

class indexTranscript(luigi.Task):

	project_name=GlobalParameter().project_name
	genome_name = GlobalParameter().genome_name
	transcriptome_name=GlobalParameter().transcriptome_name
	transcriptome_suffix=GlobalParameter().transcriptome_suffix
	genome_suffix=GlobalParameter().genome_suffix
	read_library_type =GlobalParameter().read_library_type
	adapter=GlobalParameter().adapter

	#annotation_file_type=GlobalParameter().annotation_suffix

	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	quant_method = luigi.ChoiceParameter(choices=["salmon", "kallisto"], var_type=str)
	Salmon_Index_Parameter = luigi.Parameter(default="-k 31")
	Kallisto_Index_Parameter = luigi.Parameter(default="-k 31")
	sampleName = luigi.Parameter()


	def requires(self):

		if (self.pre_process_reads == "yes"):
			return [bbduk(sampleName=self.sampleName),
				makeTx2Gene()]
		if (self.pre_process_reads == "no"):
			return [reformat(sampleName=self.sampleName),
				   makeTx2Gene()]
		

	def output(self):

		transcriptome_dir = os.path.join(os.getcwd(),GlobalParameter().transcriptome_dir + "/")
		transcriptome_index_dir = os.path.join(os.getcwd(), self.project_name, "AlignmentFreeDEA","transcript_index", self.transcriptome_name + "_" + self.quant_method +"_index" + "/")

		if (self.quant_method == "salmon") and (self.read_library_type == "pe"):
			return {'out1': luigi.LocalTarget(transcriptome_index_dir + "seq.bin"),
					'out2': luigi.LocalTarget(transcriptome_index_dir + "versionInfo.json")}

		if (self.quant_method == "salmon") and (self.read_library_type == "se"):
			return {'out1': luigi.LocalTarget(transcriptome_index_dir + "seq.bin"),
					'out2': luigi.LocalTarget(transcriptome_index_dir + "versionInfo.json")}

		if (self.quant_method == "kallisto") and (self.read_library_type == "pe"):
			return {'out': luigi.LocalTarget(transcriptome_index_dir + "/" + "kallisto.idx")}

		if (self.quant_method == "kallisto") and (self.read_library_type == "se"):
			return {'out': luigi.LocalTarget(transcriptome_index_dir + "/" + "kallisto.idx")}

	def run(self):
		transcriptome_dir = os.path.join(os.getcwd(),GlobalParameter().transcriptome_dir + "/")
		genome_dir = os.path.join(os.getcwd(), GlobalParameter().genome_dir + "/")
		transcriptome_index_dir = os.path.join(os.getcwd(), self.project_name,"AlignmentFreeDEA", "transcript_index", self.transcriptome_name + "_" + self.quant_method +"_index" + "/")


		cmd_run_salmon_index_pe = "[ -d  {transcriptome_index_dir} ] || mkdir -p {transcriptome_index_dir}; " \
								  "salmon index {Salmon_Index_Parameter} -t {transcriptome_dir}{transcriptome_name}.{transcriptome_suffix} " \
								  "-i {transcriptome_index_dir}" \
			.format(transcriptome_index_dir=transcriptome_index_dir,
					transcriptome_name=self.transcriptome_name,
					transcriptome_suffix=self.transcriptome_suffix,
					transcriptome_dir=transcriptome_dir,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_salmon_index_se = "[ -d  {transcriptome_index_dir} ] || mkdir -p {transcriptome_index_dir}; " \
								  "salmon index {Salmon_Index_Parameter} -t {transcriptome_dir}{transcriptome_name}.{transcriptome_suffix} " \
								  "-i {transcriptome_index_dir}" \
			.format(transcriptome_index_dir=transcriptome_index_dir,
					transcriptome_name=self.transcriptome_name,
					transcriptome_suffix=self.transcriptome_suffix,
					transcriptome_dir=transcriptome_dir,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_kallisto_index_pe = "[ -d  {transcriptome_index_dir} ] || mkdir -p {transcriptome_index_dir}; " \
									"cd {transcriptome_index_dir}; " \
									"kallisto index {Kallisto_Index_Parameter} " \
									"--index=kallisto.idx {transcriptome_dir}{transcriptome_name}.{transcriptome_suffix} " \
			.format(transcriptome_index_dir=transcriptome_index_dir,
					transcriptome_name=self.transcriptome_name,
					transcriptome_suffix=self.transcriptome_suffix,
					transcriptome_dir=transcriptome_dir,
					Kallisto_Index_Parameter=self.Kallisto_Index_Parameter)

		cmd_run_kallisto_index_se = "[ -d  {transcriptome_index_dir} ] || mkdir -p {transcriptome_index_dir}; " \
									"cd {transcriptome_index_dir}; " \
									"kallisto index {Kallisto_Index_Parameter} " \
									"--index=kallisto.idx {transcriptome_dir}{transcriptome_name}.{transcriptome_suffix} " \
			.format(transcriptome_index_dir=transcriptome_index_dir,
					transcriptome_name=self.transcriptome_name,
					transcriptome_suffix=self.transcriptome_suffix,
					transcriptome_dir=transcriptome_dir,
					Kallisto_Index_Parameter=self.Kallisto_Index_Parameter)

		if (self.quant_method == "salmon") and (self.read_library_type == "pe"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_pe)
			print (run_cmd(cmd_run_salmon_index_pe))


		if (self.quant_method == "salmon") and (self.read_library_type == "se"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_se)
			print (run_cmd(cmd_run_salmon_index_se))



		if (self.quant_method == "kallisto") and (self.read_library_type == "pe"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_index_pe)
			print (run_cmd(cmd_run_kallisto_index_pe))


		if (self.quant_method == "kallisto") and (self.read_library_type == "se"):
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_kallisto_index_se)
			print (run_cmd(cmd_run_kallisto_index_se))
