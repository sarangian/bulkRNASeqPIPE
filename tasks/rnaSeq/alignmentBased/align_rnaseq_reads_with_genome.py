import luigi
import os
import time
import subprocess
from tasks.rnaSeq.alignmentBased.index_genome import indexGenome
from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.reFormatReads import reformat

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
	
	rnaseq_dir=luigi.Parameter()
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

#################################################################################################################################################
class alignReadsToGenome(luigi.Task):

	project_name=GlobalParameter().project_name
	genome_name=GlobalParameter().genome_name
	read_library_type=GlobalParameter().read_library_type

	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star","hisat2","dart", "segemehl","bowtie2"],var_type=str)
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	organism_domain=GlobalParameter().domain
	#annotation_file_type = GlobalParameter().annotation_suffix
	sampleName=luigi.Parameter()
	threads=GlobalParameter().threads
	genome_suffix=GlobalParameter().genome_suffix
	maxMemory=GlobalParameter().maxMemory

	adapter = GlobalParameter().adapter
	#only_align = luigi.BoolParameter(default=False, parsing=luigi.BoolParameter.EXPLICIT_PARSING)



	def requires(self):
		if self.read_library_type == "pe" and self.pre_process_reads == "yes":
			return [indexGenome(rnaseq_aligner=self.rnaseq_aligner),
								
					[bbduk(sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]]


		if self.read_library_type == "pe" and self.pre_process_reads == "no":
			return [indexGenome(rnaseq_aligner=self.rnaseq_aligner),
								
					[reformat(sampleName=i)
					 for i in [line.strip()
							   for line in
							   open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]]

		if self.read_library_type == "se" and self.pre_process_reads == "yes":
			return [indexGenome(rnaseq_aligner=self.rnaseq_aligner),
								
					[bbduk(sampleName=i)
					 for i in [line.strip()
							   for line in
							   open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]]

		if self.read_library_type == "se" and self.pre_process_reads == "no":
			return [indexGenome(rnaseq_aligner=self.rnaseq_aligner),
								
					[reformat(sampleName=i)
					 for i in [line.strip()
							   for line in
							   open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]]



	def output(self):


		#mapFolder = os.path.join(os.getcwd(), self.project_name, "rnaseq_genome_alignment",self.genome_name + "_" + self.rnaseq_aligner + "_map" + "/")
		#if self.only_align == True:
			#mapFolder = os.path.join(os.getcwd(), self.project_name, "rnaseq_genome_alignment",
								 #self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")
		#else:
		mapFolder = os.path.join(os.getcwd(), self.project_name, "AlignmentBasedDEA", "genome_map",
									 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

		#STAR OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star",self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star",self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}


		#HISAT2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}


		#BOWTIE2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie22", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}


		#DART OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}



		# SEGEMEHL OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}
		
		# SUBREAD OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.sampleName + ".bam")}




	def run(self):
		project_name=GlobalParameter().project_name

		if self.pre_process_reads=="no" and self.read_library_type=="pe":
			cleanedReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Verified_PE_Reads" + "/")
		if self.pre_process_reads=="no" and self.read_library_type=="se":
			cleanedReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Verified_SE_Reads" + "/")

		if self.pre_process_reads=="yes" and self.read_library_type=="pe":
			cleanedReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Cleaned_PE_Reads" + "/")
		if self.pre_process_reads=="yes" and self.read_library_type=="se":
			cleanedReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Cleaned_SE_Reads" + "/")


		
		refGenomeFolder = os.path.join(GlobalParameter().genome_dir + "/")

		genomeIndexFolder = os.path.join(os.getcwd(), self.project_name,"AlignmentBasedDEA", "genome_index",self.genome_name + "_" + self.rnaseq_aligner + "_index" + "/")

		mapFolder = os.path.join(os.getcwd(), self.project_name,"AlignmentBasedDEA","genome_map",
									 self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")


		qualimapFolder = os.path.join(os.getcwd(),self.project_name,"AlignmentBasedDEA", "QualiMAP",
									  self.rnaseq_aligner + "_qualimap_" + self.read_library_type, self.sampleName + "/")

		##########################################################################################################
		# 1 SEGEMEHL rnaseq_aligner                                                                                          #
		##########################################################################################################

		cmd_run_segemehl_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "segemehl.x -t {threads} " \
								   "-i {genomeIndexFolder}index.idx " \
									"-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
								   "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
								   "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
								   "-S " \
								   "|samtools view -bS - | samtools sort " \
								   "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t {threads} " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					genome_name=self.genome_name,genome_suffix=self.genome_suffix,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t {threads} " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-p {cleanedReadFolder}{sampleName}_R2.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)

		cmd_run_segemehl_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
									  "cd {mapFolder}; " \
									  "segemehl.x -t {threads} " \
									  "-i {genomeIndexFolder}index.idx " \
									  "-d {refGenomeFolder}{genome_name}.{genome_suffix} " \
									  "-q {cleanedReadFolder}{sampleName}_R1.fastq " \
									  "-S " \
									  "|samtools view -bS - | samtools sort " \
									  "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder,
					refGenomeFolder=refGenomeFolder)




		################################################################################################################
		# 1. HISAT2 rnaseq_aligner
		################################################################################################################
		cmd_run_hisat2_pe_prok = "[-d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder};hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								"-p {threads} " \
								 "--max-intronlen 20 " \
								 "--no-spliced-alignment " \
								"-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								"-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								"| samtools view -bS - | samtools sort " \
								"-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					sampleName=self.sampleName,
					genome_name=self.genome_name)

		cmd_run_hisat2_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								"cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								"-p {threads} " \
								"--max-intronlen 20 " \
								"--no-spliced-alignment " \
								"-U {cleanedReadFolder}{sampleName}.fastq " \
								"| samtools view -bS - | samtools sort " \
								"-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)


		cmd_run_hisat2_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								 "-p {threads} " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								 "| samtools view -bS - | samtools sort " \
								 "-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)

		cmd_run_hisat2_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; hisat2 --dta -x {genomeIndexFolder}{genome_name} " \
								 "-p {threads} " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "| samtools view -bS - | samtools sort -m 5G " \
								 "-o {mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					genomeIndexFolder=genomeIndexFolder,
					cleanedReadFolder=cleanedReadFolder,
					genome_name=self.genome_name,
					sampleName=self.sampleName)



		################################################################################################################
		# 2. STAR rnaseq_aligner
		################################################################################################################
		cmd_run_star_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "STAR --runMode alignReads " \
							  "--alignIntronMax 1 " \
							  "--outSAMtype BAM " \
							  "SortedByCoordinate " \
							  "--limitBAMsortRAM 16000000000 " \
							  "--alignSJDBoverhangMin 999 " \
							  "--runThreadN {threads} " \
							  "--genomeDir {genomeIndexFolder} " \
							  "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							  "{cleanedReadFolder}{sampleName}_R2.fastq " \
							  "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
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
							   "--runThreadN {threads} " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
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
							   "--runThreadN {threads} " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}_R1.fastq " \
							   "{cleanedReadFolder}{sampleName}_R2.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					sampleName=self.sampleName,
					threads=GlobalParameter().threads,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_star_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							   "cd {mapFolder}; " \
							   "STAR --runMode alignReads " \
							   "--outSAMstrandField intronMotif " \
							   "--outSAMtype BAM " \
							   "SortedByCoordinate " \
							   "--limitBAMsortRAM 16000000000 " \
							   "--runThreadN {threads} " \
							   "--genomeDir {genomeIndexFolder} " \
							   "--readFilesIn {cleanedReadFolder}{sampleName}.fastq " \
							   "--outFileNamePrefix {sampleName}_ " \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)


		cmd_star_bam_rename = "mv {mapFolder}{sampleName}_Aligned.sortedByCoord.out.bam " \
							  "{mapFolder}{sampleName}.bam" \
			.format(mapFolder=mapFolder, sampleName=self.sampleName)


		##########################################################################################################
		#3 DART rnaseq_aligner                                                                                          #
		##########################################################################################################

		cmd_run_dart_map_pe_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "dart -intron 1 " \
							  "-t {threads} " \
							  "-i {genomeIndexFolder}index " \
							  "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
							  "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
							  "-j {mapFolder}{sampleName}_junctions.tab " \
							  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)


		cmd_run_dart_map_se_prok = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
							  "cd {mapFolder}; " \
							  "dart -intron 1 " \
							  "-t {threads} " \
							  "-i {genomeIndexFolder}index " \
							  "-f {cleanedReadFolder}{sampleName}.fastq " \
							  "-j {mapFolder}{sampleName}_junctions.tab " \
							  "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_dart_map_pe_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "dart -t {threads} " \
								   "-i {genomeIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}_R1.fastq " \
								   "-f2 {cleanedReadFolder}{sampleName}_R2.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=GlobalParameter().threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		cmd_run_dart_map_se_euk = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								   "cd {mapFolder}; " \
								   "dart -t {threads} " \
								   "-i {genomeIndexFolder}index " \
								   "-f {cleanedReadFolder}{sampleName}.fastq " \
								   "-j {mapFolder}{sampleName}_junctions.tab " \
								   "-bo {mapFolder}{sampleName}_unsorted.bam" \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)

		#Sort Unsorted BAM file and rempove the unsorted bam file

		cmd_run_dart_sort_bam = "cd {mapFolder}; samtools sort -@8 {sampleName}_unsorted.bam -o {sampleName}.bam".format(
			mapFolder=mapFolder, sampleName=self.sampleName)

		cmd_run_dart_remove_unsorted_bam = "cd {mapFolder}; rm {sampleName}_unsorted.bam".format(
			mapFolder=mapFolder, sampleName=self.sampleName)

		################################################################################################################
		#4 Bowtie2 rnaseq_aligner
		################################################################################################################

		cmd_run_bowtie2_map_pe = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "bowtie2 -x {genomeIndexFolder}index " \
								 "-p {threads} " \
								 "-a " \
								 "-1 {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-2 {cleanedReadFolder}{sampleName}_R2.fastq |samtools view -bS - | samtools " \
								 "sort " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)


		cmd_run_bowtie2_map_se = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "bowtie2 -x {genomeIndexFolder}index " \
								 "-p {threads} " \
								 "-a " \
								 "-U {cleanedReadFolder}{sampleName}.fastq " \
								 "|samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)
		
		
		################################################################################################################
		#5 SUBREAD rnaseq_aligner
		################################################################################################################

		cmd_run_subread_map_pe = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "subread-align -i {genomeIndexFolder}index " \
								 "-T {threads} " \
								 "-t 0 " \
								 "-r {cleanedReadFolder}{sampleName}_R1.fastq " \
								 "-R {cleanedReadFolder}{sampleName}_R2.fastq " \
								 "--SAMoutput | samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
					sampleName=self.sampleName,
					cleanedReadFolder=cleanedReadFolder,
					genomeIndexFolder=genomeIndexFolder)


		cmd_run_subread_map_se = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; " \
								 "cd {mapFolder}; " \
								 "subread-align -i {genomeIndexFolder}index " \
								 "-T {threads} " \
								 "-t 0 " \
								 "-r {cleanedReadFolder}{sampleName}.fastq " \
								 "--SAMoutput | samtools view -@2 -bS - | samtools sort -@2 " \
								 "-o {mapFolder}{sampleName}.bam " \
			.format(mapFolder=mapFolder,
					threads=self.threads,
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
							"--java-mem-size={maxMemory}G " \
							"-outdir {qualimapFolder} " \
							"-outfile {sampleName}_QualiMap " \
							"-outformat PDF:HTML" \
			.format(qualimapFolder=qualimapFolder,
					maxMemory=self.maxMemory,
					mapFolder=mapFolder,
					sampleName=self.sampleName)

########################################################################################################################
# Call rnaseq_aligner commands
########################################################################################################################
#Run Bowtie2: Only for prokaryotes
		if all([self.read_library_type == "pe", self.organism_domain == "prokaryote", self.rnaseq_aligner == "bowtie2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_pe)
			print (run_cmd(cmd_run_bowtie2_map_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if all([self.read_library_type == "se", self.organism_domain == "prokaryote", self.rnaseq_aligner == "bowtie2"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_map_se)
			print (run_cmd(cmd_run_bowtie2_map_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))
			
#Run Subraed
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_pe)
			print (run_cmd(cmd_run_subread_map_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))
			
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_pe)
			print (run_cmd(cmd_run_subread_map_pe))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_se)
			print (run_cmd(cmd_run_subread_map_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))
			
		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_map_se)
			print (run_cmd(cmd_run_subread_map_se))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

# Run Segmehl
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_euk)
			print (run_cmd(cmd_run_segemehl_map_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_euk)
			print (run_cmd(cmd_run_segemehl_map_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_pe_prok)
			print (run_cmd(cmd_run_segemehl_map_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "segemehl") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_map_se_prok)
			print (run_cmd(cmd_run_segemehl_map_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

#Run DART

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_euk)
			print (run_cmd(cmd_run_dart_map_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print (run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print (run_cmd(cmd_run_dart_remove_unsorted_bam))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_se_euk)
			print (run_cmd(cmd_run_dart_map_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print (run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print (run_cmd(cmd_run_dart_remove_unsorted_bam))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_map_pe_prok)
			print (run_cmd(cmd_run_dart_map_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_sort_bam)
			print (run_cmd(cmd_run_dart_sort_bam))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_remove_unsorted_bam)
			print (run_cmd(cmd_run_dart_remove_unsorted_bam))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "dart") and (self.organism_domain == "prokaryote"):
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
		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_euk)
			print (run_cmd(cmd_run_hisat2_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_euk)
			print (run_cmd(cmd_run_hisat2_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_pe_prok)
			print (run_cmd(cmd_run_hisat2_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "hisat2") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_se_prok)
			print (run_cmd(cmd_run_hisat2_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

#Run STAR

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_prok)
			print (run_cmd(cmd_run_star_pe_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "prokaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_prok)
			print (run_cmd(cmd_run_star_se_prok))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "pe") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_pe_euk)
			print (run_cmd(cmd_run_star_pe_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))

		if (self.read_library_type == "se") and (self.rnaseq_aligner == "star") and (self.organism_domain == "eukaryote"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_se_euk)
			print (run_cmd(cmd_run_star_se_euk))

			print("****** NOW RUNNING COMMAND ******: " + cmd_star_bam_rename)
			print (run_cmd(cmd_star_bam_rename))

			print("****** NOW RUNNING COMMAND ******: " + cmd_run_qualimap)
			print (run_cmd(cmd_run_qualimap))


#####################################################################################################
class  alignReadSetsToGenome(luigi.Task):

	#Global Parameters
	project_name = GlobalParameter().project_name
	read_library_type = GlobalParameter().read_library_type
	adapter = luigi.Parameter(default="./tasks/utility/adapters.fasta.gz")
	genome_name = GlobalParameter().genome_name
	organism_domain = GlobalParameter().domain
	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory

	#Local Parameters
	#sampleName = luigi.Parameter()
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star", "hisat2", "dart", "segemehl", "bowtie2"], var_type=str)
	#annotation_file_type = GlobalParameter().annotation_name
	#only_align = luigi.BoolParameter(default=True, parsing=luigi.BoolParameter.EXPLICIT_PARSING)

	#####################################################################################################3

	def requires(self):
		if self.read_library_type=="pe":
			return [alignReadsToGenome(pre_process_reads=self.pre_process_reads,
						   rnaseq_aligner=self.rnaseq_aligner,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]					   

		if self.read_library_type=="se":
			return [alignReadsToGenome(pre_process_reads=self.pre_process_reads,
						   rnaseq_aligner=self.rnaseq_aligner,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]
	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(), "task_logs", 'task.generate.count.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('Count File Generation finished at {t}'.format(t=timestamp))

			
		#annotation_file_type=self.annotation_file_type,

	'''	
	def output(self):

		mapFolder = os.path.join(os.getcwd(), self.project_name, "genome_map",self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

			# STAR OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

			# HISAT2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

			# BOWTIE2 OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "bowtie2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

			# DART OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "dart", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

			# SEGEMEHL OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "segemehl", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

			# SUBREAD OUTPUT
		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "pe", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

		if all([self.read_library_type == "se", self.rnaseq_aligner == "subread", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(mapFolder + "/" + self.genome_name + ".bam")}

	def run(self):
	
	
		mapFolder = os.path.join(os.getcwd(), self.project_name, "genome_map",self.genome_name + "_" + self.rnaseq_aligner + "_" + self.read_library_type + "_map" + "/")

		
		cmd_merge_bam = "[ -d {mapFolder} ] || mkdir -p {mapFolder}; cd {mapFolder}; " \
						"samtools merge {genome_name}.bam *.bam " \
				.format(mapFolder=mapFolder, genome_name=self.genome_name)
		print("****** NOW RUNNING COMMAND ******: " + cmd_merge_bam)
		print(run_cmd(cmd_merge_bam))
	'''