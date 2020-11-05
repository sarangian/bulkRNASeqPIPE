import luigi
import os
import time
import subprocess
from tasks.rnaSeq.alignmentBased.index_genome import indexGenome
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.reFormatReads import reformat
from tasks.rnaSeq.alignmentBased.align_rnaseq_reads_with_genome import alignReadsToGenome
from tasks.rnaSeq.alignmentBased.align_rnaseq_reads_with_genome import alignReadSetsToGenome

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


class featureCounts(luigi.Task):
	#Global Parameters
	project_name = GlobalParameter().project_name
	read_library_type = GlobalParameter().read_library_type
	threads = GlobalParameter().threads
	genome_name = GlobalParameter().genome_name
	organism_domain = GlobalParameter().domain
	feature_type = GlobalParameter().feature_type

	#Local parameters
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	attribute_type = luigi.Parameter(default="gene_id",description='''Specify attribute type in GTF annotation. 
											 string(=[gene_id])''')

	annotation_dir=GlobalParameter().annotation_dir
	annotation_name=GlobalParameter().annotation_name
	annotation_file_type=GlobalParameter().annotation_suffix
	strandType = luigi.ChoiceParameter(default="0",choices=['0','1','2'],description='''perform strand-specific read counting. int([=0]unstranded) 
															OR [=1] stranded] OR [=2] reversely-stranded. default[
															=0]''' )


	def requires(self):

		return [alignReadSetsToGenome(pre_process_reads=self.pre_process_reads,						 
						 rnaseq_aligner=self.rnaseq_aligner)

			   ]

			   #annotation_file_type=self.annotation_file_type,

	def output(self):
		QuantFolder = os.path.join(os.getcwd(), self.project_name,
								   "AlignmentBasedDEA",
								   "ReadQuant",
								   self.rnaseq_aligner + "_" + self.read_library_type + "_" + "quant" + "/")

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}
		
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "subread") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "subread") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "subread") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "subread") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}



		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}




		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}




		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}



		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "bowtie2") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "bowtie2") and (self.organism_domain == "prokaryote"):
			return {'out': luigi.LocalTarget(QuantFolder + "/" + "counts.txt")}




	def run(self):

		QuantFolder = os.path.join(os.getcwd(), GlobalParameter().project_name,
								   "AlignmentBasedDEA",
								   "ReadQuant",
								   self.rnaseq_aligner + "_" + self.read_library_type + "_" + "quant" + "/")

		genomeIndexFolder = os.path.join(os.getcwd(), GlobalParameter().project_name, "AlignmentBasedDEA","genome_index",self.genome_name + "_" + self.rnaseq_aligner + "_index" + "/")
		mapFolder = os.path.join(os.getcwd(), GlobalParameter().project_name, "AlignmentBasedDEA","genome_map",self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")
		
		if self.annotation_file_type=="gff":

			cmd_gff2gtf =  "gffread -E  {annotation_dir}{annotation_name}.gff -T -o {annotation_dir}{annotation_name}.gtf " \
							.format(annotation_dir=self.annotation_dir,
							annotation_name=self.annotation_name)

			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print (run_cmd(cmd_gff2gtf))



		cmd_run_featureCount_pe = "[ -d {QuantFolder} ] || mkdir -p {QuantFolder}; " \
										"cd {QuantFolder}; " \
										"featureCounts -a {annotation_dir}{annotation_name}.gtf " \
										"-t {feature_type} " \
										"-g {attribute_type} " \
										"-s {strandType} " \
										"-T {threads} " \
										"-p " \
										"-o {QuantFolder}counts.txt " \
										"{mapFolder}*.bam " \
			.format(QuantFolder=QuantFolder,
					feature_type=self.feature_type,
					strandType=self.strandType,
					attribute_type=self.attribute_type,
					threads=self.threads,
					annotation_dir=GlobalParameter().annotation_dir,
					annotation_name=GlobalParameter().annotation_name,
					mapFolder=mapFolder)

		cmd_run_featureCount_se = "[ -d {QuantFolder} ] || mkdir -p {QuantFolder}; " \
									  "cd {QuantFolder}; " \
									  "featureCounts -a {genomeIndexFolder}{genome_name}.gtf " \
									  "-t {feature_type} " \
									  "-g {attribute_type} " \
									  "-s {strandType} " \
									  "-T {threads} " \
									  "-o counts.txt " \
									  "{mapFolder}*.bam " \
			.format(QuantFolder=QuantFolder,
					feature_type=self.feature_type,
					attribute_type=self.attribute_type,
					strandType=self.strandType,
					threads=self.threads,
					genomeIndexFolder=genomeIndexFolder,
					mapFolder=mapFolder,
					genome_name=self.genome_name)



		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))
			
		
		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))




		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))



		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))


		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_se)
			print (run_cmd(cmd_run_featureCount_se))

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_featureCount_pe)
			print (run_cmd(cmd_run_featureCount_pe))



##################################################################
class alignmentBasedQuant(luigi.Task):
	project_name = luigi.Parameter(default="RNASeqAnalysis")
	read_library_type = GlobalParameter().read_library_type
	threads = GlobalParameter().threads
	genome_name = GlobalParameter().genome_name
	adapter = GlobalParameter().adapter
	organism_domain = GlobalParameter().domain
	feature_type = GlobalParameter().feature_type
	annotation_file_type = GlobalParameter().annotation_suffix


	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star","hisat2","dart", "segemehl","bowtie2"],var_type=str)

	attribute_type = luigi.Parameter(default="gene_id", description='''Specify attribute type in GTF annotation. 
												 string(=[gene_id])''')

	strandType = luigi.ChoiceParameter(default="0", choices=['0', '1', '2'], description='''perform strand-specific read counting. int([=0]unstranded) 
																OR [=1] stranded] OR [=2] reversely-stranded. default[
																=0]''')
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)


	def requires(self):
		return [featureCounts(pre_process_reads=self.pre_process_reads,					 
					 rnaseq_aligner=self.rnaseq_aligner,
					 attribute_type=self.attribute_type,
					 strandType=self.strandType)
				]

				#annotation_file_type=self.annotation_file_type,

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(), "task_logs", 'task.generate.count.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('Count File Generation finished at {t}'.format(t=timestamp))

