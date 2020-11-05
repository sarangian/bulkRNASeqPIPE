import luigi
import os
import time
import subprocess
from tasks.rnaSeq.annotation.prokaryotic_annotation import annotateGenome


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
	domain=luigi.Parameter()
	genome_suffix=luigi.Parameter()
	genome_name=luigi.Parameter()
	genome_dir=luigi.Parameter()
	
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

class makeTx2Gene(luigi.Task):

	#Global Parameters
	project_name=GlobalParameter().project_name
	genome_name=GlobalParameter().genome_name
	genome_suffix = GlobalParameter().genome_suffix
	organism_domain=GlobalParameter().domain

	annotation_name=GlobalParameter().annotation_name
	annotation_file_type = GlobalParameter().annotation_suffix
	annotation_dir=GlobalParameter().annotation_dir

	genome_dir=GlobalParameter().genome_dir
	genome_name=GlobalParameter().genome_name
	genome_suffix=GlobalParameter().genome_suffix

	transcriptome_dir=GlobalParameter().transcriptome_dir
	transcriptome_name=GlobalParameter().transcriptome_name
	transcriptome_suffix=GlobalParameter().transcriptome_suffix



	threads = GlobalParameter().threads
	adapter = GlobalParameter().adapter

	minContigLength=luigi.Parameter(default="200")


	def requires(self):
		if (self.organism_domain == "prokaryote") and (self.annotation_file_type == "NA"):
			return [annotateGenome(minContigLength=self.minContigLength)]

	def output(self):

		transcriptome_dir = os.path.join(GlobalParameter().transcriptome_dir + "/")
		annotation_dir=os.path.join(GlobalParameter().annotation_dir+ "/")

		return {'out1': luigi.LocalTarget(annotation_dir + "tx2gene.csv"),
			   'out2': luigi.LocalTarget(annotation_dir + self.annotation_name + ".gtf"),
			   'out3': luigi.LocalTarget(transcriptome_dir + self.transcriptome_name + ".ffn")
				}

	def run(self):
		transcriptome_dir = os.path.join(GlobalParameter().transcriptome_dir + "/")
		genome_dir = os.path.join(os.getcwd(), GlobalParameter().genome_dir + "/")
		annotation_dir=os.path.join(GlobalParameter().annotation_dir+ "/")

		cmd_extract_transcript_from_genome = "[ -d  {transcriptome_dir} ] || mkdir -p {transcriptome_dir}; " \
							  "gffread -w {transcriptome_dir}{transcriptome_name}.{transcriptome_suffix} " \
							   "-g {genome_dir}{genome_name}.{genome_suffix} {annotation_dir}{annotation_name}.gtf -F " \
							   .format(genome_dir=genome_dir,
									   transcriptome_dir=transcriptome_dir,
									   annotation_dir=annotation_dir,annotation_name=self.annotation_name,
									   transcriptome_name=self.transcriptome_name,
									   transcriptome_suffix=self.transcriptome_suffix,
									   genome_name=self.genome_name,
									   genome_suffix=self.genome_suffix)
		
		cmd_remove_genome_index = "rm {genome_dir}{genome_name}.{genome_suffix}.fai".format(genome_dir=genome_dir, genome_name=self.genome_name, genome_suffix=self.genome_suffix)

		
		cmd_gff2gtf =   "gffread -E  {annotation_dir}{annotation_name}.gff -T -o {annotation_dir}{annotation_name}.gtf " \
								  .format(annotation_dir=annotation_dir,
								  	annotation_name=self.annotation_name)

		cmd_tx2gene_for_eukaryote = "cd {annotation_dir}; tx2gene.R " \
							   "-a gtf " \
							   "-p {annotation_dir}{annotation_name}.gtf " \
							   "-o {annotation_dir}tx2gene.csv" \
							   .format(annotation_dir=annotation_dir,
							   		   annotation_name=self.annotation_name
									  )


		tx2g_cmd1 = ''' awk '$3=="transcript"' | cut -f9 | tr -s ";" " " |awk '{print$4","$2}' ''' 
		tx2g_cmd2 = ''' sed 's/\"//g' | sed -e '1i\TRANSCRIPT,GENE' '''

		
		cmd_tx2gene_for_prokaryote = "[ -d  {annotation_dir} ] || mkdir -p {annotation_dir}; cd {annotation_dir};" \
								  "cat {annotation_dir}{annotation_name}.gtf |{tx2g_cmd1} | {tx2g_cmd2} > tx2gene.csv" \
								  .format(annotation_dir=annotation_dir,
									  	  annotation_name=self.annotation_name,
									      tx2g_cmd1=tx2g_cmd1,
										  tx2g_cmd2=tx2g_cmd2)
		

		if (self.annotation_file_type == "gtf") and (self.organism_domain == "eukaryote"):

						
			print ("****** NOW RUNNING COMMAND ******: " + cmd_extract_transcript_from_genome)
			print (run_cmd(cmd_extract_transcript_from_genome))

			print ("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_for_eukaryote)
			print (run_cmd(cmd_tx2gene_for_eukaryote))
			
		
		if (self.annotation_file_type == "gff") and (self.organism_domain == "eukaryote"):
			
			print ("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print (run_cmd(cmd_gff2gtf))

				
			print ("****** NOW RUNNING COMMAND ******: " + cmd_extract_transcript_from_genome)
			print (run_cmd(cmd_extract_transcript_from_genome))

			print ("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_for_eukaryote)
			print (run_cmd(cmd_tx2gene_for_eukaryote))
			

			print ("****** NOW RUNNING COMMAND ******: " + cmd_remove_genome_index)
			print (run_cmd(cmd_remove_genome_index))
			
		
		if (self.annotation_file_type == "gtf") and (self.organism_domain == "prokaryote"):			
					
			print ("****** NOW RUNNING COMMAND ******: " + cmd_extract_transcript_from_genome)
			print (run_cmd(cmd_extract_transcript_from_genome))

			print ("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_for_prokaryote)
			print (run_cmd(cmd_tx2gene_for_prokaryote))


		if (self.annotation_file_type == "gff") and (self.organism_domain == "prokaryote"):
			
			print ("****** NOW RUNNING COMMAND ******: " + cmd_gff2gtf)
			print (run_cmd(cmd_gff2gtf))
					
			print ("****** NOW RUNNING COMMAND ******: " + cmd_extract_transcript_from_genome)
			print (run_cmd(cmd_extract_transcript_from_genome))

			print ("****** NOW RUNNING COMMAND ******: " + cmd_tx2gene_for_prokaryote)
			print (run_cmd(cmd_tx2gene_for_prokaryote))

			print ("****** NOW RUNNING COMMAND ******: " + cmd_remove_genome_index)
			print (run_cmd(cmd_remove_genome_index))
