import luigi
import os
import time
import subprocess
import pandas as pd
from luigi import Parameter

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

def rockhopper_pe_input_bbduk(input_file):
	df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
	df1.columns = ["SampleName", "Condition"]
	df2 = df1.set_index("SampleName", drop=False)
	df_groups = df2.groupby('Condition')

	grp_list = []
	for name, group in df_groups:
		grp_list.append(name)
	grp_str = ','.join(grp_list)
	parsed_string = "-L "

	group_name_list = []
	pe_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Cleaned_PE_Reads" + "/")
	for name,group in df_groups:
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'
		left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in group.SampleName.tolist()]
		right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in group.SampleName.tolist()]
		result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]
		result1 = [x+"%" +y for x, y in result]
		group_string = ','.join(result1)
		group_name_list.append(group_string)
	parse_group_string = " ".join(group_name_list)
	parsed_string = "{} {} {}".format(parsed_string, grp_str, parse_group_string)
	return parsed_string

def rockhopper_se_input_bbduk(input_file):
	df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
	df1.columns = ["SampleName", "Condition"]
	df2 = df1.set_index("SampleName", drop=False)
	df_groups = df2.groupby('Condition')
	grp_list = []
	se_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Cleaned_SE_Reads" + "/")
	for name, group in df_groups:
		grp_list.append(name)
	grp_str = ','.join(grp_list)
	parsed_string = "-L "
	group_name_list = []
	for name,group in df_groups:
		read_name_suffix = '.fastq'
		read_name_list = [se_cleaned_read_folder + x + read_name_suffix for x in group.SampleName.tolist()]
		result1 = [x for x in read_name_list]
		group_string = ','.join(result1)
		group_name_list.append(group_string)
	parse_group_string = " ".join(group_name_list)
	parsed_string = "{} {} {}".format(parsed_string, grp_str, parse_group_string)
	return parsed_string

############################################
def rockhopper_pe_input_reformat(input_file):
	df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
	df1.columns = ["SampleName", "Condition"]
	df2 = df1.set_index("SampleName", drop=False)
	df_groups = df2.groupby('Condition')

	grp_list = []
	for name, group in df_groups:
		grp_list.append(name)
	grp_str = ','.join(grp_list)
	parsed_string = "-L "

	group_name_list = []
	pe_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC","Verified_PE_Reads" + "/")
	for name,group in df_groups:
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'
		left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in group.SampleName.tolist()]
		right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in group.SampleName.tolist()]
		result = [sublist for sublist in zip(left_read_name_list, right_read_name_list)]
		result1 = [x+"%" +y for x, y in result]
		group_string = ','.join(result1)
		group_name_list.append(group_string)
	parse_group_string = " ".join(group_name_list)
	parsed_string = "{} {} {}".format(parsed_string, grp_str, parse_group_string)
	return parsed_string

def rockhopper_se_input_reformat(input_file):
	df1 = pd.read_table(input_file, sep='\t+', engine='python', header=None)
	df1.columns = ["SampleName", "Condition"]
	df2 = df1.set_index("SampleName", drop=False)
	df_groups = df2.groupby('Condition')
	grp_list = []
	se_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC","Verified_SE_Reads" + "/")
	for name, group in df_groups:
		grp_list.append(name)
	grp_str = ','.join(grp_list)
	parsed_string = "-L "
	group_name_list = []
	for name,group in df_groups:
		read_name_suffix = '.fastq'
		read_name_list = [se_cleaned_read_folder + x + read_name_suffix for x in group.SampleName.tolist()]
		result1 = [x for x in read_name_list]
		group_string = ','.join(result1)
		group_name_list.append(group_string)
	parse_group_string = " ".join(group_name_list)
	parsed_string = "{} {} {}".format(parsed_string, grp_str, parse_group_string)
	return parsed_string
############################################
def rockhopper_tab_to_fasta_pe(input_file):
	data = pd.read_table(input_file, sep='\t+', engine='python')
	transcripts = data[data.Length >= int(100)]

	df = transcripts.assign(Transcript_ID=[1 + i for i in range(len(transcripts))])[['Transcript_ID'] + transcripts.columns.tolist()]

	transcript_file=os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", "rockhopper","transcripts.fna")
	with open(transcript_file, 'w') as f:
		for index, row in df.iterrows():
			f.write(">transcript_" + str(row['Transcript_ID']) + "_length_" + str(row['Length']))
			f.write("\n")
			f.write(row['Sequence'])
			f.write("\n")

def rockhopper_tab_to_fasta_se(input_file):
	data = pd.read_table(input_file, sep='\t+', engine='python')
	transcripts = data[data.Length >= int(100)]

	df = transcripts.assign(Transcript_ID=[1 + i for i in range(len(transcripts))])[['Transcript_ID'] + transcripts.columns.tolist()]

	transcript_file=os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", "rockhopper","transcripts.fna")
	with open(transcript_file, 'w') as f:
		for index, row in df.iterrows():
			f.write(">transcript_" + str(row['Transcript_ID']) + "_length_" + str(row['Length']))
			f.write("\n")
			f.write(row['Sequence'])
			f.write("\n")


############################################
class GlobalParameter(luigi.Config):
	genome_suffix=luigi.Parameter()
	read_library_type=luigi.Parameter()
	domain=luigi.Parameter()
	genome_name=luigi.Parameter()
	genome_dir=luigi.Parameter()
	transcriptome_dir=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	feature_type=luigi.Parameter()
	adapter=luigi.Parameter()
	project_name=luigi.Parameter()

class rockhopper(luigi.Task):
	project_name=GlobalParameter().project_name
	adapter = GlobalParameter().adapter
	read_library_type = GlobalParameter().read_library_type
	organism_domain = GlobalParameter().domain
	threads = GlobalParameter().threads
	genome_name = GlobalParameter().genome_name
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)


	def requires(self):
		if self.pre_process_reads=="yes" and self.read_library_type=="pe":
			return [bbduk(
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]
		if self.pre_process_reads=="yes" and self.read_library_type=="se":
			return [bbduk(
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]
		if self.pre_process_reads=="no" and self.read_library_type=="pe":
			return [reformat(
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]
		if self.pre_process_reads=="no" and self.read_library_type=="se":
			return [reformat(
					  sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "config", "se_samples.lst")))]]

	def output(self):

		assembled_transcript_folder = os.path.join(os.getcwd(), self.project_name, "deNovoDEA","denovo_assembly", "rockhopper" + "/")
		return {'out_1': luigi.LocalTarget(assembled_transcript_folder + "transcripts.txt"),
				'out_2': luigi.LocalTarget(assembled_transcript_folder + "transcripts.fna")
				}


	def run(self):
		assembled_transcript_folder = os.path.join(os.getcwd(), self.project_name, "deNovoDEA","denovo_assembly", "rockhopper" + "/")
		input_group_file = os.path.join(os.getcwd(), "config", "group.tsv")

		if self.pre_process_reads == "yes" and self.read_library_type == "pe":
			ReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC", "Cleaned_PE_Reads" + "/")
			cmd_rockhopper_read_input = rockhopper_pe_input_bbduk(input_group_file)

		if self.pre_process_reads == "yes" and self.read_library_type == "se":
			ReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC", "Cleaned_SE_Reads" + "/")
			cmd_rockhopper_read_input = rockhopper_se_input_bbduk(input_group_file)

		if self.pre_process_reads == "no" and self.read_library_type == "pe":
			ReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Verified_PE_Reads" + "/")
			cmd_rockhopper_read_input = rockhopper_pe_input_reformat(input_group_file)
		if self.pre_process_reads == "no" and self.read_library_type == "se":
			ReadFolder = os.path.join(os.getcwd(), self.project_name,"ReadQC","Verified_SE_Reads" + "/")
			cmd_rockhopper_read_input = rockhopper_se_input_reformat(input_group_file)



		#Command to run Rockhopper
		cmd_run_rockhopper = "[ -d  {assembled_transcript_folder} ] || mkdir -p {assembled_transcript_folder}; " \
					   "cd {ReadFolder};" \
					   "java -Xmx12G " \
					   "-cp $(which Rockhopper.jar) Rockhopper  {cmd_rockhopper_read_input} " \
					   "-s false " \
					   "-p {threads} " \
					   "-o {assembled_transcript_folder}" \
			.format(assembled_transcript_folder=assembled_transcript_folder,
					ReadFolder=ReadFolder,
					threads=GlobalParameter().threads,
					cmd_rockhopper_read_input=cmd_rockhopper_read_input)



		print("****** NOW RUNNING COMMAND ******: " + cmd_run_rockhopper)
		print (run_cmd(cmd_run_rockhopper))

		time.sleep(2)
		transcript_tab_file = os.path.join(assembled_transcript_folder, "transcripts.txt").format(
			assembled_transcript_folder=assembled_transcript_folder)
		if self.read_library_type=="pe":
			rockhopper_tab_to_fasta_pe(transcript_tab_file)
		if self.read_library_type=="se":
			rockhopper_tab_to_fasta_se(transcript_tab_file)






