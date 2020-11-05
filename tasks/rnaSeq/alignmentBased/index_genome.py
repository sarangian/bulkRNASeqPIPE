import luigi
import os
import subprocess
from Bio import SeqIO
import math
from tasks.rnaSeq.annotation.prokaryotic_annotation import annotateGenome

def run_cmd(cmd):
    p = subprocess.Popen(cmd, bufsize=-1,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
    output = p.communicate()[0]
    return output


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

def genomeSAindexNbases(genome_fasta):
    with open(genome_fasta) as genome:
        totalsize=0
        for rec in SeqIO.parse(genome, 'fasta'):
            totalsize = totalsize + len(rec)
        log2 = math.log(totalsize, 2.0)
        index = (log2/2) - 1
        gsanb = int(index)
        return gsanb

class indexGenome(luigi.Task):
	project_name = GlobalParameter().project_name
	genome_name = GlobalParameter().genome_name
	genome_suffix= GlobalParameter().genome_suffix
	genome_dir=GlobalParameter().genome_dir
	organism_domain = GlobalParameter().domain
	read_library_type= GlobalParameter().read_library_type
	

	annotation_name=GlobalParameter().annotation_name
	annotation_dir=GlobalParameter().annotation_dir
	annotation_suffix=GlobalParameter().annotation_suffix


	rnaseq_aligner = luigi.ChoiceParameter(choices=["subread","star","hisat2","dart", "segemehl","bowtie2"],var_type=str)
	annotation_file_type = GlobalParameter().annotation_suffix

	def requires(self):

		if (self.organism_domain == "prokaryote") and (self.annotation_file_type == "NA"):
			return [annotateGenome()]


	def output(self):
		genomeIndexFolder=os.path.join(os.getcwd(), self.project_name, "AlignmentBasedDEA","genome_index", self.genome_name +"_"+ self.rnaseq_aligner +"_index" + "/")

		if all([self.read_library_type== "pe", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}

		if all([self.read_library_type== "se", self.rnaseq_aligner == "star", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}

		if all([self.read_library_type== "pe", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}

		if all([self.read_library_type== "se", self.rnaseq_aligner == "star", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "Genome"),
					'out2': luigi.LocalTarget(genomeIndexFolder + "SAindex")}

		if all([self.read_library_type== "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genome_name + ".1.ht2")}
		if all([self.read_library_type== "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "prokaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genome_name + ".1.ht2")}
		if all([self.read_library_type== "pe", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genome_name + ".1.ht2")}
		if all([self.read_library_type== "se", self.rnaseq_aligner == "hisat2", self.organism_domain == "eukaryote"]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + self.genome_name + ".1.ht2")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "segemehl")]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "index.idx")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "subread")]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "index.files")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "dart")]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "index.bwt")}

		if all([((self.read_library_type== "pe") or (self.read_library_type== "se")),
				((self.organism_domain == "prokaryote") or (self.organism_domain == "eukaryote")),
				(self.rnaseq_aligner == "bowtie2")]):
			return {'out1': luigi.LocalTarget(genomeIndexFolder + "index.1.bt2")}

	def run(self):
		genome_dir = os.path.join(GlobalParameter().genome_dir)
		genome_name = GlobalParameter().genome_name
		genome_suffix = GlobalParameter().genome_suffix

		annotation_dir=os.path.join(GlobalParameter().annotation_dir)
		#annotation_name=GlobalParameter().annotation_name
		genomeIndexFolder=os.path.join(os.getcwd(), self.project_name, "AlignmentBasedDEA","genome_index",self.genome_name +"_"+ self.rnaseq_aligner +"_index" + "/")

		genomeFastaFile = "{genome_dir}{genome_name}.{genome_suffix}".format(genome_dir=genome_dir,
																	  genome_name=self.genome_name,genome_suffix=self.genome_suffix)

		cmd_gff2gtf = "gffread -E  {annotation_dir}{annotation_name}.gff -T " \
					  "-o {annotation_dir}{annotation_name}.gtf " \
			.format(annotation_dir=annotation_dir,
					annotation_name=GlobalParameter().annotation_name)

		
		gsan = genomeSAindexNbases(genomeFastaFile)

		cmd_run_star_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder}; STAR " \
							 "--runMode genomeGenerate " \
							 "--genomeSAindexNbases {gsan} " \
							 "--genomeFastaFiles {genome_dir}{genome_name}.{genome_suffix} " \
							 "--genomeDir {genomeIndexFolder} " \
			.format(genomeIndexFolder=genomeIndexFolder,
					genome_name=self.genome_name,genome_suffix=self.genome_suffix,
					gsan=gsan,
					genome_dir=genome_dir)

		cmd_run_bowtie2_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
								"bowtie2-build -f {genome_dir}{genome_name}.{genome_suffix} {genomeIndexFolder}index " \
			.format(genomeIndexFolder=genomeIndexFolder,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					genome_dir=genome_dir)

		cmd_run_dart_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
							 "bwt_index {genome_dir}{genome_name}.{genome_suffix} {genomeIndexFolder}index " \
			.format(genomeIndexFolder=genomeIndexFolder,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					genome_dir=genome_dir)

		cmd_run_hisat2_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
							   "hisat2-build {genome_dir}{genome_name}.{genome_suffix} {genomeIndexFolder}{genome_name} " \
			.format(genomeIndexFolder=genomeIndexFolder,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					genome_dir=genome_dir)

		cmd_run_segemehl_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
								 "segemehl.x -x {genomeIndexFolder}index.idx -d {genome_dir}{genome_name}.{genome_suffix}  " \
			.format(genomeIndexFolder=genomeIndexFolder,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					genome_dir=genome_dir)

		cmd_run_subread_index = "[ -d  {genomeIndexFolder} ] || mkdir -p {genomeIndexFolder};" \
								"subread-buildindex -o {genomeIndexFolder}index {genome_dir}{genome_name}.{genome_suffix} " \
			.format(genomeIndexFolder=genomeIndexFolder,genome_suffix=self.genome_suffix,
					genome_name=self.genome_name,
					genome_dir=genome_dir)

		if (self.rnaseq_aligner == "star") and (self.annotation_file_type == "gtf"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_index)
			print(run_cmd(cmd_run_star_index))
			
		if (self.rnaseq_aligner == "star") and (self.annotation_file_type == "gff"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_star_index)
			print(run_cmd(cmd_run_star_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "dart") and (self.annotation_file_type == "gtf"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_index)
			print(run_cmd(cmd_run_dart_index))
			
		if (self.rnaseq_aligner == "dart") and (self.annotation_file_type == "gff"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_dart_index)
			print(run_cmd(cmd_run_dart_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "bowtie2") and (self.annotation_file_type == "gtf"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_index)
			print(run_cmd(cmd_run_bowtie2_index))
			
		if (self.rnaseq_aligner == "bowtie2") and (self.annotation_file_type == "gff"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_bowtie2_index)
			print(run_cmd(cmd_run_bowtie2_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "hisat2") and (self.annotation_file_type == "gtf"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_index)
			print(run_cmd(cmd_run_hisat2_index))
			

		if (self.rnaseq_aligner == "hisat2") and (self.annotation_file_type == "gff"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_hisat2_index)
			print(run_cmd(cmd_run_hisat2_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "segemehl") and (self.annotation_file_type == "gtf"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_index)
			print(run_cmd(cmd_run_segemehl_index))
			
		if (self.rnaseq_aligner == "segemehl") and (self.annotation_file_type == "gff"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_segemehl_index)
			print(run_cmd(cmd_run_segemehl_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))

		if (self.rnaseq_aligner == "subread") and (self.annotation_file_type == "gtf"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_index)
			print(run_cmd(cmd_run_subread_index))
			
		if (self.rnaseq_aligner == "subread") and (self.annotation_file_type == "gff"):
			print("****** NOW RUNNING COMMAND ******: " + cmd_run_subread_index)
			print(run_cmd(cmd_run_subread_index))
			print("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print(run_cmd(cmd_gff2gtf))