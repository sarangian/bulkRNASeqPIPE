import luigi
import os
import time
import subprocess
import pandas as pd
from tasks.rnaSeq.alignmentBased.generate_gene_count_file import alignmentBasedQuant


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
	genome_suffix=luigi.Parameter()
	read_library_type=luigi.Parameter()
	domain=luigi.Parameter()
	genome_name=luigi.Parameter()
	genome_dir=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	feature_type=luigi.Parameter()
	adapter=luigi.Parameter()
	#result_tag=luigi.Parameter()

	annotation_suffix=luigi.Parameter()
	project_name=luigi.Parameter()

class alignmentBasedDEA(luigi.Task):
	adapter = GlobalParameter().adapter
	organism_domain = GlobalParameter().domain
	read_library_type = GlobalParameter().read_library_type
	#result_tag = GlobalParameter().result_tag
	feature_type = GlobalParameter().feature_type
	project_name=GlobalParameter().project_name

	# Local parameters
	
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread", "star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	attribute_type = luigi.ChoiceParameter(choices=["gene_id", "transcript_id"],default="gene_id", description='''Specify attribute type in GTF annotation. 
												 string(=[gene_id])''', var_type=str)
	annotation_file_type = GlobalParameter().annotation_suffix


	strand_ype = luigi.ChoiceParameter(default="0", choices=['0', '1', '2'], 
									   description='''perform strand-specific read counting. int([=0]unstranded) 
										OR [=1] stranded] OR [=2] reversely-stranded. default[=0]''')

	dea_method = luigi.ChoiceParameter(default="deseq2",choices=["deseq2", "edger"], var_type=str)
	report_name = luigi.Parameter(default="DEA_Report")
	factor_of_intrest = luigi.Parameter(default="conditions", description="factor of intrest column of the target file (string [=condititions]). ")
	reference_condition = luigi.Parameter(default="control", description="reference biological condition.  (string [=control]")
	# target_file = luigi.Parameter(description="path to the design/target file. (string [=target.tsv]")
	alpha = luigi.Parameter(default="0.05", description="threshold of statistical significance.  (float [=0.05]")
	p_adjust_method = luigi.ChoiceParameter(default="BH", description="p-value adjustment method.", choices=["BH", "BY"], var_type=str)
	fit_type = luigi.ChoiceParameter(default="parametric", description="mean-variance relationship.", choices=["parametric", "local", "mean"],
									 var_type=str)
	size_factor = luigi.ChoiceParameter(default="median", description="method to estimate the size factors.", choices=["median", "shorth"],
										var_type=str)
	result_tag = luigi.Parameter(default="treated_vs_control")

#Local Parameter


	def requires(self):
		return [alignmentBasedQuant(rnaseq_aligner=self.rnaseq_aligner,
						
						pre_process_reads=self.pre_process_reads)
			   ]

			   #annotation_file_type=self.annotation_file_type,

	def output(self):

		resultFolder = os.path.join(os.getcwd(), self.project_name,"AlignmentBasedDEA",
									"DEAnalysis",
									self.dea_method + "_" + self.rnaseq_aligner + "_" + self.result_tag + "_" + self.read_library_type + "/")

		return {'out1': luigi.LocalTarget(resultFolder + self.report_name + "/" + "index.html")}


	def run(self):

		target_file = os.path.join(os.getcwd(),"config","target.tsv")
		#rmd_DESeq2File = os.path.expanduser(os.path.join(('~'), 'scriptome', 'tasks', 'utility', "PlotDESEQ2.Rmd"))
		#rmd_edgeRFile = os.path.expanduser(os.path.join(('~'), 'scriptome', 'tasks', 'utility', "PlotEDGER.Rmd"))



		resultFolder = os.path.join(os.getcwd(), self.project_name,"AlignmentBasedDEA",
									"DEAnalysis",
									self.dea_method + "_" + self.rnaseq_aligner + "_" + self.result_tag + "_" + self.read_library_type + "/")

		QuantFolder = os.path.join(os.getcwd(),GlobalParameter().project_name,
								   "AlignmentBasedDEA",
								   "ReadQuant",
								   self.rnaseq_aligner + "_" + self.read_library_type + "_" + "quant" + "/")



		#rmd_DESeq2File = GlobalParameter().utilityFolder + "PlotDESEQ2.Rmd"
		#rmd_edgeRFile = GlobalParameter().utilityFolder + "PlotEDGER.Rmd"


		cmd_run_featuteCount_DESeq2 = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
										"cd {resultFolder};" \
										"featureCount_DESeq2.r " \
										"-t {target_file} " \
										"-q {QuantFolder}counts.txt " \
										"-v {factor_of_intrest} " \
										"-c {reference_condition} " \
										"-f {fit_type} " \
										"-a {alpha} " \
										"-p {p_adjust_method} " \
										"-l {size_factor} " \
										"-T $(which PlotDESEQ2.Rmd)" \
											.format(target_file=target_file,
													resultFolder=resultFolder,
													QuantFolder=QuantFolder,
													factor_of_intrest=self.factor_of_intrest,
													reference_condition=self.reference_condition,
													fit_type=self.fit_type,
													alpha=self.alpha,
													p_adjust_method=self.p_adjust_method,
													size_factor=self.size_factor)

		cmd_run_featuteCount_edgeR = "[ -d  {resultFolder} ] || mkdir -p {resultFolder}; " \
									  "cd {resultFolder};" \
									  "featureCount_edgeR.r " \
									  "-t {target_file} " \
									  "-q {QuantFolder}counts.txt " \
									  "-v {factor_of_intrest} " \
									  "-c {reference_condition} " \
									  "-f {fit_type} " \
									  "-a {alpha} " \
									  "-p {p_adjust_method} " \
									  "-l {size_factor} " \
									  "-T $(which PlotEDGER.Rmd)" \
			.format(resultFolder=resultFolder,
					target_file=target_file,
					QuantFolder=QuantFolder,
					factor_of_intrest=self.factor_of_intrest,
					reference_condition=self.reference_condition,
					fit_type=self.fit_type,
					alpha=self.alpha,
					p_adjust_method=self.p_adjust_method,
					size_factor=self.size_factor)

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))


		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))


		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))



		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))


		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))
			
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))


		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote", self.dea_method=="deseq2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_DESeq2)
			print (run_cmd(cmd_run_featuteCount_DESeq2))

####################################################################################################################
#Run EDGER
######################################################################################################################
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))


		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))



		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote",  self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))



		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote",  self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))


		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))
			
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))


		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote", self.dea_method=="edger"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featuteCount_edgeR)
			print (run_cmd(cmd_run_featuteCount_edgeR))
