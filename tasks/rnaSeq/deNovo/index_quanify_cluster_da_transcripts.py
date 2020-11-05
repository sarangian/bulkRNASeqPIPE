import luigi
import os
import time
import subprocess
import pandas as pd
from tasks.rnaSeq.deNovo.transcriptome_assembly import dnTransAssemble

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

def prepare_corset_input(input_file):
    df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
    df1.columns = ["SampleName", "Condition"]
    df2 = df1.set_index("SampleName", drop = False)
    df_groups = df2.groupby('Condition')

    condition_group = ",".join(df1['Condition'])
    sample_group = ",".join(df1['SampleName'])
    CorSet_Input = '-g {} -n {}'.format(condition_group, sample_group)
    return CorSet_Input


class GlobalParameter(luigi.Config):
	genome_suffix=luigi.Parameter()
	read_library_type=luigi.Parameter()
	domain=luigi.Parameter()
	genome_name=luigi.Parameter()
	transcriptome_dir=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	adapter=luigi.Parameter()
	project_name=luigi.Parameter()

class indexDnTranscript(luigi.Task):

	project_name=GlobalParameter().project_name
	genome_name = GlobalParameter().genome_name
	read_library_type =GlobalParameter().read_library_type
	organism_domain = GlobalParameter().domain
	adapter=GlobalParameter().adapter
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	#annotation_file_type = luigi.ChoiceParameter(choices=["GFF", "GTF", "NA"], var_type=str)
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	Salmon_Index_Parameter = luigi.Parameter(default="--type quasi -k 31")	
	#sampleName = luigi.Parameter()


	def requires(self):
		return [dnTransAssemble(
					pre_process_reads=self.pre_process_reads,
					rnaseq_assembler=self.rnaseq_assembler)]
		

	def output(self):

		TranscriptIndexFolder = os.path.join(os.getcwd(), self.project_name,
											 "deNovoDEA",
											 "transcript_index",
											 self.genome_name +  "_" + self.rnaseq_assembler + "/")

		if self.rnaseq_assembler=="rockhopper":
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
				'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}

		if self.rnaseq_assembler=="spades":
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
				'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}
		if self.rnaseq_assembler=="trinity":
			return {'out1': luigi.LocalTarget(TranscriptIndexFolder + "hash.bin"),
				'out2': luigi.LocalTarget(TranscriptIndexFolder + "versionInfo.json")}


		
	def run(self):

		if self.rnaseq_assembler == "trinity":
			assembled_transcript = os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", "trinity","Trinity.fasta")
		
		if self.rnaseq_assembler == "rockhopper":
			assembled_transcript=os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", "rockhopper","transcripts.fna")
	
		if self.rnaseq_assembler == "spades":
			assembled_transcript = os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", "spades","transcripts.fasta")

		TranscriptIndexFolder = os.path.join(os.getcwd(), self.project_name,
											 "deNovoDEA",
											 "transcript_index",
											 self.genome_name +  "_" + self.rnaseq_assembler + "/")

		cmd_run_salmon_index_pe = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {assembled_transcript} " \
								  "-i {TranscriptIndexFolder} " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					assembled_transcript=assembled_transcript,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)

		cmd_run_salmon_index_se = "[ -d  {TranscriptIndexFolder} ] || mkdir -p {TranscriptIndexFolder}; " \
								  "salmon index -t {assembled_transcript} " \
								  "-i {TranscriptIndexFolder} " \
			.format(TranscriptIndexFolder=TranscriptIndexFolder,
					assembled_transcript=assembled_transcript,
					Salmon_Index_Parameter=self.Salmon_Index_Parameter)


		if self.read_library_type == "pe":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_pe)
			print (run_cmd(cmd_run_salmon_index_pe))


		if self.read_library_type == "se":
			print ("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_index_se)
			print (run_cmd(cmd_run_salmon_index_se))



##Run Salmom Quantification
class denovoQuant(luigi.Task):
	sampleName = luigi.Parameter()
	project_name=GlobalParameter().project_name
	organism_domain = GlobalParameter().domain
	read_library_type = GlobalParameter().read_library_type
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	adapter=GlobalParameter().adapter
	genome_name=GlobalParameter().genome_name
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	threads= GlobalParameter().threads
	memory=GlobalParameter().maxMemory


	def requires(self):
		return [indexDnTranscript(pre_process_reads=self.pre_process_reads,
							      rnaseq_assembler=self.rnaseq_assembler)]

	def output(self):
		salmon_quant_folder = os.path.join(os.getcwd(), self.project_name,"deNovoDEA", "ReadQuant",
									 self.rnaseq_assembler, "SalmonQuant" + "/")

		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			return {'out1': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "quant.sf"),
					'out2': luigi.LocalTarget(salmon_quant_folder + self.sampleName + "/" + "cmd_info.json")
					}

	def run(self):

		salmon_quant_folder = os.path.join(os.getcwd(), self.project_name,"deNovoDEA", "ReadQuant",
									 self.rnaseq_assembler, "SalmonQuant" + "/")

		salmon_index_folder = os.path.join(os.getcwd(), self.project_name,
											 "deNovoDEA",
											 "transcript_index",
											 self.genome_name +  "_" + self.rnaseq_assembler + "/")

		
		salmon_map_folder = os.path.join(os.getcwd(), self.project_name,
											 "deNovoDEA",
											 "transcript_map",
											 self.genome_name +  "_" + self.rnaseq_assembler + "/")


		if self.pre_process_reads == "yes" and self.read_library_type == "pe":
			clean_read_folder =os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Cleaned_PE_Reads" + "/")

		if self.pre_process_reads == "yes" and self.read_library_type == "se":
			clean_read_folder = os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Cleaned_PE_Reads" + "/")

		if self.pre_process_reads == "no" and self.read_library_type == "pe":
			clean_read_folder =os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Verified_PE_Reads" + "/")

		if self.pre_process_reads == "no" and self.read_library_type == "se":
			clean_read_folder =os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Verified_SE_Reads" + "/")



		cmd_run_salmon_quant_pe = "[ -d {salmon_quant_folder} ] || mkdir -p {salmon_quant_folder}; " \
								  "[ -d {salmon_map_folder} ] || mkdir -p {salmon_map_folder}; " \
								  "cd {salmon_map_folder}; " \
								  "salmon quant --no-version-check -p {threads} " \
								  "-i {salmon_index_folder} " \
								  "-l A " \
								  "-1 {clean_read_folder}{sampleName}_R1.fastq " \
								  "-2 {clean_read_folder}{sampleName}_R2.fastq " \
								  "--dumpEq " \
								  "--output {salmon_quant_folder}{sampleName} " \
								  "--validateMappings --hardFilter " \
								  "--writeMappings | samtools view -bS - | samtools sort -m {memory}G " \
								  "-o {salmon_map_folder}{sampleName}.bam" \
			.format(salmon_quant_folder=salmon_quant_folder,
					salmon_map_folder=salmon_map_folder,
					sampleName=self.sampleName,
					clean_read_folder=clean_read_folder,
					threads=self.threads,
					salmon_index_folder=salmon_index_folder,
					memory=self.memory)

		cmd_run_salmon_quant_se = "[ -d {salmon_quant_folder} ] || mkdir -p {salmon_quant_folder}; " \
								  "mkdir -p {salmon_map_folder}; " \
								  "cd {salmon_map_folder};" \
								  "salmon quant --no-version-check -p {threads} " \
								  "--seqBias " \
								  "--gcBias " \
								  "-i {salmon_index_folder} " \
								  "-l A " \
								  "-r {clean_read_folder}{sampleName}.fastq " \
								  "--dumpEq " \
								  "--output {salmon_quant_folder}{sampleName} " \
								  "--validateMappings --hardFilter " \
								  "--writeMappings | samtools view -bS - | samtools sort -m {memory}G " \
								  "-o {salmon_map_folder}{sampleName}.bam" \
			.format(salmon_quant_folder=salmon_quant_folder,
					salmon_map_folder=salmon_map_folder,
					sampleName=self.sampleName,
					clean_read_folder=clean_read_folder,
					salmon_index_folder=salmon_index_folder,memory=self.memory,
					threads=self.threads)

		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))



		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))



		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))


		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))



		if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_pe)
			print (run_cmd(cmd_run_salmon_quant_pe))



		if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_salmon_quant_se)
			print (run_cmd(cmd_run_salmon_quant_se))


##########################################################################################################################

class quantifyDnTranscript(luigi.Task):

	project_name=GlobalParameter().project_name
	organism_domain = GlobalParameter().domain
	read_library_type = GlobalParameter().read_library_type
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	adapter=GlobalParameter().adapter
	genome_name=GlobalParameter().genome_name
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	threads= GlobalParameter().threads



	def requires(self):
		if self.read_library_type =="pe":
			return [denovoQuant(rnaseq_assembler=self.rnaseq_assembler,pre_process_reads=self.pre_process_reads,
							
							sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]

		if self.read_library_type =="se":
			return [denovoQuant(rnaseq_assembler=self.rnaseq_assembler,
							pre_process_reads=self.pre_process_reads,
							sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget('workflow.complete.{t}'.format(t=timestamp))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('workflow finished at {t}'.format(t=timestamp))

#####################################################################################################################################################3
# Run Corsetanalysis
class clusterDnTranscript(luigi.Task):
	project_name=GlobalParameter().project_name
	rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
	organism_domain = GlobalParameter().domain
	read_library_type = GlobalParameter().read_library_type
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	adapter=GlobalParameter().adapter

	def requires(self):

		return [
				[quantifyDnTranscript(
									 pre_process_reads=self.pre_process_reads,
									 rnaseq_assembler=self.rnaseq_assembler)]
				]


	def output(self):
		corset_folder = os.path.join(os.getcwd(), self.project_name,"deNovoDEA", "ReadQuant",
									 self.rnaseq_assembler, "Corset/" + "/")

		return {'out1': luigi.LocalTarget(corset_folder + "clusters.txt"),
				'out2': luigi.LocalTarget(corset_folder + "counts.txt")
				}


	def run(self):
		corset_folder = os.path.join(os.getcwd(), self.project_name,"deNovoDEA", "ReadQuant",
									 self.rnaseq_assembler, "Corset" + "/")


		salmon_quant_folder = os.path.join(os.getcwd(), self.project_name,"deNovoDEA", "ReadQuant",
									 self.rnaseq_assembler, "SalmonQuant" + "/")

		assembled_transcript_folder = os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", self.rnaseq_assembler + "/")
		super_transcript_folder = os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", self.rnaseq_assembler , "SuperTranscript" + "/")

		# Paired end reads

		# Command to generate Rockhopper Paired-end input
		input_group_file = os.path.join(os.getcwd(),"config", "group.tsv")

		# Command to generate Rockhopper Paired-end input
		cmd_corset_read_input = prepare_corset_input(input_group_file)

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
			.format(super_transcript_folder=super_transcript_folder, corset_folder=corset_folder,
					assembled_transcript_folder=assembled_transcript_folder)

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset)
		print(run_cmd(cmd_run_corset))

		print("****** NOW RUNNING COMMAND ******: " + cmd_run_supertrans)
		print(run_cmd(cmd_run_supertrans))



	class denovoDEA(luigi.Task):
		adapter = GlobalParameter().adapter
		organism_domain = GlobalParameter().domain
		read_library_type = GlobalParameter().read_library_type
		result_tag = luigi.Parameter(default="treated_vs_control")

		project_name = GlobalParameter().project_name
		rnaseq_assembler = luigi.ChoiceParameter(choices=["trinity", "spades", "rockhopper"], var_type=str)
		pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
		dea_method = luigi.ChoiceParameter(choices=["deseq2", "edger"], var_type=str)
		report_name = luigi.Parameter(default="Corset_DESeq2_HTML_Report")
		factor_of_intrest = luigi.Parameter(default="condititions",description="factor of intrest column of the target file (string [=condititions]). ")
		reference_condition = luigi.Parameter(default="control", description="reference biological condition.  (string [=control]")
		#target_file = luigi.Parameter(description="path to the design/target file. (string [=target.tsv]")
		alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
		p_adjust_method = luigi.ChoiceParameter(default="BH", description="p-value adjustment method.",choices=["BH", "BY"],var_type=str)
		fit_type = luigi.ChoiceParameter(default="parametric", description="mean-variance relationship.",  choices=["parametric", "local","mean"],var_type=str)
		size_factor = luigi.ChoiceParameter(default="median", description="method to estimate the size factors.", choices=["median", "shorth"],var_type=str)


		def requires(self):
			return [clusterDnTranscript(rnaseq_assembler=self.rnaseq_assembler,
							   pre_process_reads=self.pre_process_reads)]

		def output(self):

			#edgeResultFolter = GlobalParameter().basefolder + "/" + self.projectName + "_DEAnalysis/" + self.rnaseq_assembler + "_" +  self.sampleGroupFile + "_" + self.readType + "/" + "edgeR/"
			resultFolder = os.path.join(os.getcwd(), self.project_name,
										"deNovoDEA",
										"DEAnalysis",
										self.dea_method + "_" + self.rnaseq_assembler + "_" +  self.result_tag +"/")

			if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "deseq2"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "deseq2"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "deseq2"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "deseq2"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "deseq2"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "deseq2"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}


			if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "edger"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "edger"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "edger"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "edger"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "edger"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "edger"]):
				return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}


		def run(self):
			resultFolder = os.path.join(os.getcwd(), self.project_name,
										"deNovoDEA",
										"DEAnalysis",
										self.dea_method + "_" + self.rnaseq_assembler + "_" +  self.result_tag +"/")
			corset_folder = os.path.join(os.getcwd(), self.project_name,"deNovoDEA", "ReadQuant",
									 self.rnaseq_assembler, "Corset" + "/")

			target_file = os.path.join(os.getcwd(),"config","target.tsv")
			rmd_DESeq2File = os.path.expanduser(os.path.join(('~'), 'scriptome','tasks','utility',"PlotDESEQ2.Rmd"))
			rmd_edgeRFile = os.path.expanduser(os.path.join(('~'), 'scriptome', 'tasks', 'utility', "PlotEDGER.Rmd"))

			cmd_run_corset_deseq = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
							  "cd {resultFolder};" \
							  "corset_DESeq2.r " \
							  "-t {target_file} " \
							  "-q {corset_folder}counts.txt " \
							   "-v {factor_of_intrest} " \
							   "-c {reference_condition} " \
							   "-f {fit_type} " \
							   "-a {alpha} " \
							   "-p {p_adjust_method} " \
							   "-l {size_factor} " \
								   "-T {rmd_DESeq2File}" \
							.format(resultFolder=resultFolder,
									target_file=target_file,
									corset_folder=corset_folder,
									factor_of_intrest=self.factor_of_intrest,
									reference_condition=self.reference_condition,
									fit_type=self.fit_type,
									alpha=self.alpha,
									p_adjust_method=self.p_adjust_method,
									size_factor=self.size_factor,
									rmd_DESeq2File=rmd_DESeq2File)

			cmd_run_corset_edger = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
								   "cd {resultFolder};" \
								   "corset_edgeR.r " \
								   "-t {target_file} " \
								   "-q {corset_folder}counts.txt " \
								   "-v {factor_of_intrest} " \
								   "-c {reference_condition} " \
								   "-f {fit_type} " \
								   "-a {alpha} " \
								   "-p {p_adjust_method} " \
								   "-l {size_factor} " \
								   "-T {rmd_edgeRFile}" \
				.format(resultFolder=resultFolder,
						target_file=target_file,
						corset_folder=corset_folder,
						factor_of_intrest=self.factor_of_intrest,
						reference_condition=self.reference_condition,
						fit_type=self.fit_type,
						alpha=self.alpha,
						p_adjust_method=self.p_adjust_method,
						size_factor=self.size_factor,
						rmd_edgeRFile=rmd_edgeRFile)

			if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "deseq2"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
				print (run_cmd(cmd_run_corset_deseq))

			if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "deseq2"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
				print (run_cmd(cmd_run_corset_deseq))

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "deseq2"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
				print (run_cmd(cmd_run_corset_deseq))

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "deseq2"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
				print (run_cmd(cmd_run_corset_deseq))

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "deseq2"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
				print (run_cmd(cmd_run_corset_deseq))

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "deseq2"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_deseq)
				print (run_cmd(cmd_run_corset_deseq))

	##############################################################################################################
	#Run EDGER
			if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "edger"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
				print (run_cmd(cmd_run_corset_edger))

			if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_assembler == "rockhopper",
					self.dea_method == "edger"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
				print (run_cmd(cmd_run_corset_edger))

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "edger"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
				print (run_cmd(cmd_run_corset_edger))

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "trinity",
					self.dea_method == "edger"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
				print (run_cmd(cmd_run_corset_edger))

			if all([self.read_library_type == "pe", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "edger"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
				print (run_cmd(cmd_run_corset_edger))

			if all([self.read_library_type == "se", self.organism_domain == "eukaryote", self.rnaseq_assembler == "spades",
					self.dea_method == "edger"]):
				print("****** NOW RUNNING COMMAND ******: " + cmd_run_corset_edger)
				print (run_cmd(cmd_run_corset_edger))