import luigi
import os
import subprocess

from tasks.readCleaning.preProcessReads import bbduk
from tasks.readCleaning.reFormatReads import reformat

class GlobalParameter(luigi.Config):
    threads = luigi.Parameter()
    maxMemory = luigi.Parameter()
    project_name=luigi.Parameter()

def run_cmd(cmd):
    p = subprocess.Popen(cmd, bufsize=-1,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
    output = p.communicate()[0]
    return output


def spades_pe_formater_bbduk(pefile):
    with open(pefile) as fh:
        sample_name_list = fh.read().splitlines()

        left_read_name_suffix = '_R1.fastq'
        left_read_name_prefix = '--pe-1 ' 
        right_read_name_suffix = '_R2.fastq'
        right_read_name_prefix = '--pe-2 '

        pe_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Cleaned_PE_Reads" + "/")
        left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        right_read_name_list =[pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]

        lib_no_list = [x for x in range (1, len(sample_name_list)+1)]

        left_prefix=[left_read_name_prefix + str(y) for y in lib_no_list]
        right_prefix=[right_read_name_prefix+ str(y) for y in lib_no_list]

        left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]

        right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]

        result = [item for sublist in zip(left_prefix, left_read_name_list, right_prefix, right_read_name_list) for item in sublist]

        parse_string = " ".join(result)

        return parse_string


#Reformat"Verified_PE_Reads"
##########################################
def spades_pe_formater_reformat(pefile):
    with open(pefile) as fh:
        sample_name_list = fh.read().splitlines()

        left_read_name_suffix = '_R1.fastq'
        left_read_name_prefix = '--pe-1 ' 
        right_read_name_suffix = '_R2.fastq'
        right_read_name_prefix = '--pe-2 '

        pe_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().project_name,"ReadQC", "Verified_PE_Reads" + "/")
        left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
        right_read_name_list =[pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]

        lib_no_list = [x for x in range (1, len(sample_name_list)+1)]

        left_prefix=[left_read_name_prefix + str(y) for y in lib_no_list]
        right_prefix=[right_read_name_prefix+ str(y) for y in lib_no_list]

        left_read_name_list = [pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]

        right_read_name_list = [pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]

        result = [item for sublist in zip(left_prefix, left_read_name_list, right_prefix, right_read_name_list) for item in sublist]

        parse_string = " ".join(result)

        return parse_string
        

class spades(luigi.Task):
    project_name = GlobalParameter().project_name
    read_library_type = luigi.Parameter(default="pe")
    pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)

    cov_cutoff = luigi.Parameter(default="off",description="coverage cutoff value (a positive float number, or 'auto', or 'off'). Default 'off'")
    kmers = luigi.Parameter(default="auto",description="comma separated list of kmers. must be odd and less than 128. Default 'auto' ")


    def requires(self):
        if self.read_library_type == "pe" and self.pre_process_reads=="yes":
            return [bbduk(
                          sampleName=i)
                    for i in [line.strip()
                              for line in
                              open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]
        if self.read_library_type == "pe" and self.pre_process_reads =="no":
            return [[reformat(
                           sampleName=i)
                     for i in [line.strip()
                               for line in
                               open((os.path.join(os.getcwd(), "config", "pe_samples.lst")))]]]

        
    def output(self):
        spades_assembly_folder = os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", "spades" + "/")
        return {'out': luigi.LocalTarget(spades_assembly_folder + "transcripts.fasta")}
        

    def run(self):
       
        spades_assembly_folder = os.path.join(os.getcwd(), GlobalParameter().project_name, "deNovoDEA","denovo_assembly", "spades" + "/")
        spades_assembly_log_folder = os.path.join(os.getcwd(),self.project_name , "log", "denovo_assembly",  "spades" + "/")


        pe_sample_list = os.path.join(os.getcwd(), "config", "pe_samples.lst")
        

        if self.pre_process_reads=="yes":
            cmd_spades_pe = spades_pe_formater_bbduk(pe_sample_list)
            
        if self.pre_process_reads == "no":
            cmd_spades_pe = spades_pe_formater_reformat(pe_sample_list)
            
        

        spades_pe_cmd = "[ -d  {spades_assembly_folder} ] || mkdir -p {spades_assembly_folder}; " \
                        "mkdir -p {spades_assembly_log_folder}; cd {spades_assembly_folder}; " \
                        "spades.py --rna " \
                        "--threads {threads} " \
                        "--cov-cutoff {cov_cutoff} " \
                        "--memory {maxMemory} " \
                        "-k {kmers} " \
                        "{cmd_spades_pe} " \
                        "-o {spades_assembly_folder} " \
                        "2>&1 | tee {spades_assembly_log_folder}spades_assembly.log " \
            .format(spades_assembly_folder=spades_assembly_folder,
                    spades_assembly_log_folder=spades_assembly_log_folder,
                    kmers=self.kmers,
                    cov_cutoff=self.cov_cutoff,
                    threads=GlobalParameter().threads,
                    maxMemory=GlobalParameter().maxMemory,
                    cmd_spades_pe=cmd_spades_pe)

        
        print("****** NOW RUNNING COMMAND ******: " + spades_pe_cmd)
        run_cmd(spades_pe_cmd)