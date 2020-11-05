import luigi
import os
import time
import subprocess
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
	
	
	transcriptome_dir=luigi.Parameter()
	transcriptome_name=luigi.Parameter()
	transcriptome_suffix=luigi.Parameter()

	annotation_dir=luigi.Parameter()
	annotation_suffix=luigi.Parameter()
	annotation_name=luigi.Parameter()

	adapter=luigi.Parameter()
	threads=luigi.Parameter()
	maxMemory=luigi.Parameter()

class annotateGenome(luigi.Task):
	#Global Parameters
	genome_dir=GlobalParameter().genome_dir
	genome_name=GlobalParameter().genome_name
	genome_suffix=GlobalParameter().genome_suffix

	transcriptome_dir=GlobalParameter().transcriptome_dir
	transcriptome_name=GlobalParameter().transcriptome_name
	transcriptome_suffix=GlobalParameter().transcriptome_suffix

	annotation_dir=GlobalParameter().annotation_dir
	annotation_name=GlobalParameter().annotation_name
	annotation_file_type=GlobalParameter().annotation_suffix

	threads = GlobalParameter().threads
	adapter = GlobalParameter().adapter

	project_name=GlobalParameter().project_name

	minContigLength=luigi.Parameter(default="200")
	

	#Local Parameters

	def requires(self):
		return []

	def output(self):
		transcriptome_dir = os.path.join(GlobalParameter().transcriptome_dir + "/")
		annotation_dir = os.path.join(GlobalParameter().annotation_dir + "/")
		genome_dir = os.path.join(GlobalParameter().genome_dir+ "/")
		return {'out1': luigi.LocalTarget(transcriptome_dir +"/"  + self.genome_name + ".ffn"),
				'out2': luigi.LocalTarget(annotation_dir +"/"  + self.genome_name + ".gff"),
				'out3': luigi.LocalTarget(annotation_dir + "/" + "tx2gene.csv")
				}

	def run(self):
		transcriptome_dir = os.path.join(GlobalParameter().transcriptome_dir+ "/")
		genome_dir = os.path.join(GlobalParameter().genome_dir+ "/")
		prokka_dir = os.path.join(os.getcwd(),self.project_name,"annotation",self.genome_name + "_PROKKA" +"/")
		annotation_dir=os.path.join(GlobalParameter().annotation_dir+ "/")


		cmd_run_prokka = "prokka {genome_dir}{genome_name}.{genome_suffix} " \
						 "--cpu {threads} " \
						 "--prefix {genome_name} " \
						 "--mincontiglen {minContigLength} " \
						 "--outdir {prokka_dir} --force --rfam" \
			.format(genome_name=self.genome_name,
					prokka_dir=prokka_dir,
					genome_suffix=self.genome_suffix,
					minContigLength=self.minContigLength,
					genome_dir=genome_dir,
					threads=self.threads)

		print ("****** NOW RUNNING COMMAND ******: " + cmd_run_prokka)
		print (run_cmd(cmd_run_prokka))

		print("Generating tx2gene")
		awk_cmd = 'BEGIN{FS="\\t"}{print ""$1"," ""$1"*" $7}'
		cmd_run_tx2gene = "[ -d  {annotation_dir} ] || mkdir -p {annotation_dir}; " \
						   "cd {prokka_dir}; " \
						  "awk '{awk}' {genome_name}.tsv > {annotation_dir}tx2gene.csv" \
			.format(prokka_dir=prokka_dir,
					annotation_dir=annotation_dir,
					genome_name=self.genome_name,
					awk=awk_cmd)

		print ("****** NOW RUNNING COMMAND ******: " + cmd_run_tx2gene)
		print (run_cmd(cmd_run_tx2gene))

		
		print("Generating GTF from GFF")
		sed_cmd2 = ''' sed 's/ID//g' '''
		awk_cmd2 = 'BEGIN{OFS="\\t"}{print $1,"PROKKA","CDS",$2,$3,".",$4,".","gene_id " $5}'
		cmd_run_gff2gtf = "cd {prokka_dir}; " \
						  "grep -v '#' {genome_name}.gff | grep 'ID=' | {sed_cmd2} | " \
						  "cut -f1,4,5,7,9 | awk '{awk_cmd2}' > {genome_name}.gtf " \
			.format(genome_dir=genome_dir,
					prokka_dir=prokka_dir,
					genome_name=self.genome_name,
					sed_cmd2=sed_cmd2,
					awk_cmd2=awk_cmd2)

		print ("****** NOW RUNNING COMMAND ******: " + cmd_run_gff2gtf)
		print (run_cmd(cmd_run_gff2gtf))


		cmd_copy_files = "cd {prokka_dir}; " \
						 "cp {genome_name}.ffn {transcriptome_dir}{transcriptome_name}.ffn; " \
						 "cp {genome_name}.gff {annotation_dir}{annotation_name}.gff; " \
						 "cp {genome_name}.gtf {annotation_dir}{annotation_name}.gtf " \
						 .format(prokka_dir=prokka_dir,
						 		genome_name=self.genome_name,
						 		transcriptome_dir=transcriptome_dir,
						 		transcriptome_name=self.transcriptome_name,
						 		annotation_dir=self.annotation_dir,
						 		annotation_name=self.annotation_name						 		
						 		)
					
		print ("****** NOW RUNNING COMMAND ******: " + cmd_copy_files)
		print (run_cmd(cmd_copy_files))