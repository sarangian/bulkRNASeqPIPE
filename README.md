Installation Instructions
================================
  
To install RNASeq Analysis Toolkit, you must have a minimum of 6 GiB free disk space and minimum of 16 GiB free RAM to test run. 

To provide an easier way to install, we provide a miniconda based installer.
Installation also requires **pre-instaled** ``git``, ``gcc``, ``cpp`` and ``zlib1g-dev``.
  
    git clone https://github.com/sarangian/bulkRNASeqPIPE.git
    cd bulkRNASeqPIPE
    chmod 755 INSTALL.sh
    ./INSTALL.sh

    
**Post Installation Instructions**
After successful installation, close the current terminal. 
In a new terminal. source the bashrc file:  ``source ~/.bashrc``
Activate ``gabtk`` environment using command: ``conda activate`` 

All the third party tools installed using conda are available at $HOME/bulkRNASeqPIPE/ [default location]
or the user specified location during the installation process.

The script to run RNASeq Analysis Pipeline is rnaseq.py is available inside the bulkRNASeqPIPE folder, that you cloned from github.

Input Files
===========

Raw RNASeq Reads
----------------

.. code-block:: none

  The raw-RNASeq reads in FASTQ format must be placed inside a folder with read permission
   Allowed ``extension`` for FASTQ reads: ``fq`` , ``fq.gz`` , ``fastq``, ``fastq.gz``
 
   For **paired-end** RNAseq reads, sample name must be suffixed with _R1. ``extension`` and _R2. ``extension`` for forward and reverse reads respectively

          *Example*

           sample1_R1.fastq 
           sample1_R2.fastq

           sample2_R1.fastq.gz 
           sample2_R2.fastq.gz

           sample3_R1.fq
           sample3_R2.fq

           sample4_R1.fq.gz 
           sample4_R2.fq.gz
          
   where  ``sample1`` , ``sample2`` , ``sample3`` , ``sample4`` are the sample names
           sample_name must be suffixed with _R1.{extension} and _R2.{extension}

   For **single-end** RNAseq reads sample name must be suffixed with ``.{extension}``


          *Example*

           example_one.fastq.gz
           example_two.fq.gz
           example_three.fastq
           example_four.fq


Genome FASTA and annotation files
----------------------------------

**Required Files**
Reference / draft genome of the organism along with gene annotation files (optional for prokaryotes mandatory for eukaryotes)  transcriptome of the organism in FASTA format with ``.ffn , .fa or .fasta`` extension (optional for prokaryotes eukaryotes) must be made available in the same folder containing RNASeq reads.

 1. The reference genome or assembly of the organism in FASTA format with ``.fna , .fa or .fasta`` extention
 2. Gene Annotation files in GFF3 or GTF format 
    Note: 1.  GFF3 or GTF files must be obtained form NCBI
          2.  GFF3/GTF file is optional if domain is prokaryote and differential expression analysis method is Alignment Free (Salmon/Kallisto)
          3.  GFF3/GTF file is mandatory if differential expression analysis method is Alignment Based (STAR / hisat2 / dart / bowtie2 /subread)

 3. The transcriptome of the organism in FASTA format with ``.ffn , .fa or .fasta`` extension (optional for both Prokaryotes and eukaryotes)
    Note: 1.  For prokaryotes, if user provides the GTF/GFF file and the corresponding transcriptome file, then pipeline will not run PROKKA annotation.
          2.  The user provided transcriptome file and corresponding GTF file will be used for transcript quantification using SALMON/KAllisto
	  3.  For prokaryotes, if user provides the reference genome/assembly of the organism in FASTA format, the pipeline will trigger prokka to annotate the reference                       genome and the corresponding transcriptome, gff, gtf and tx2gene files will be made available for downstream analysis.
	  4.  For eukaryotes, user need to provide the gene annotation file in GTF/GFF3 format along with the corresponding reference genome or genome assembly. The pipeline                   will use gffread program to extract transcript form the genome using the GFF3/GTF file, which will be used for Alignment Free Differential Expression Analysis.



 Allowed ``externsion`` for 
     a. genome FASTA file: ``.fna, .fa , .fasta``
     b. transcriptome FASTA file: ``.ffn, .fa , .fasta``
     c. genome annotation file: ``.gff`` or ``.gtf``



Commands
========

1. Prepare Project
-------------------

To design and perform a RNASeq experiment, a Project need to be prepaired using the ``projectConfig.py`` script. Conda environment must be activated before running the script.

Usage:  projectConfig.py -h

    prepareProject.py <arguments>
    -h       --help             Show this help message and exit

    mandatory arguments         Description
    
    -i      --inputDir          Path to Directory containing raw RNASeq reads, annotation file (gff / gtf), genome and (or) 
                                transcriptome FASTA file
                                type: string
                                Example: $HOME/bulkRNASeqPIPE/sample_data/mastigocladus

    -r      --readType          RNASeq read type 
                                type: string
                                allowed values: [single, paired]

    -d      --domain            Organism Domain
                                type: string
                                allowed values: [prokaryote, eukaryote]

    
    Optional arguments
    ------------------
    -p     --projectDir         Name of the Project Directory to be created
                                Must not contain blank spaces and (or) special characters
                                type: string
                                Default: MyProject

    -s     --schedulerPort      Scheduler Port Number for luigi
                                type: int
                                Default: 8082

    -e     --emailAddress       Provide your email address
                                type: string
                                Default: Null

    -c     --threads            Number of threads to be used
                                type: int
                                Default = (total threads -1)

    -m     --maxMemory          Maximum allowed memory in GB. 
                                type: int
                                Default = [(available memory in GB) -1)
    

**Run Example**

    mkdir RNASeq-Analysis
    cd RNASeq-Analysis
   
    [RNASeq-Analysis]$ projectConfig.py  -i /path/to/bulkRNASeqPIPE/sample_data/ \
                      -r paired \
                      -d prokaryote \
                      -p halomicronema \
                      -o halomicronema_symlink
 
    Running the prepareProject.py script with the above 4 required parameters asks for the individual file types present inside the inputData folder.

    User has to choose  [pe:   paired end
                        se:   single end
                        geno: genome fasta file
                        tran: transcriptome fasta file
                        anot: annotation file
                        ex:   exclude file from analysis
                        ]


**Output**

    Successful run of the projectConfig.py script with appropriate parameters will generate 

    1. Luigi Configuration file ``luigi.cfg`` in the parent folder
       
      Edit the luigi.cfg file if required.
      
      Note:
      It is mandatory to provide the path of the adapter file (default location: bulkRNASeqPipe/tasks/utility/adapter.fastq.gz)

    2. a project folder in the name of ``halomicronema`` containing three files

    3. a configuration folder in the name of config containing 3 files   

      a. group.tsv
      b. pe_samples.lst
      b. samples.txt
      c. target.tsv

The ``target.tsv`` file contains the sample names with their associated biological conditions, which will be used for differential expression analysis. Kindly     check the target.tsv file and modify if required
  

Commands to run RNA-Seq Workflow
--------  ----------------------------
 

rnaseq.py <command> - -help

    Command                      Description   
    
    1. rawReadsQC                   Raw Reads Quality Assessment using FASTQC tool 
    2. cleanReads                   Process Raw Reads using BBDUK

    3. alignmentFreeQuant           Quantify transcripts using salmon (or) kallisto
    4. alignmentFreeDEA             Alignment free differential expression analysis
    
    5. mapReadsToGenome             Map RNASeq reads to the genome of the organism
    6. alignmentBasedQuant          Generates gene counts from .bam files using featureCount tool
    7. alignmentBasedDEA            Genome alignment based differential expression analysis

    8. dnTransAssemble              Denovo Assembly of Prokaryotic and Eukaryotic Transcripts
    9. quantifyDAT                  Quantify denovo assembled transcripts using Salmon
    10. clusterDAT                  Clustred Assembled transcripts based on equivalence class
    11. denovoDEA                   denovo transcriptome assembly based differential expression analysis



2. Raw reads quality assessment
--------------------------------

  **Note**
    Before running any of the RNA-Seq Workflow commands, a project must be prepared using ``projectConfig.py`` script.
    The parent forlder must have the luigi.cfg file, in which the globalparameters are defined.
    Running any of the  RNA-Seq Workflow commands without generating the project folder will give rise to ``luigi.parameter.MissingParameterException``


    **Steps**
    1. Run Prepare Projcet with project name mastigocladus as discussed before 
       and inspect the pe_samples.lst or se_samples.lst file generated inside config folder

    2. Run rawReadsQC
       [RNASeq-Analysis]$ rnaseq.py rawReadsQC --local-scheduler

      Successful execution of rawReadsQC will generate a folder ReadQC/PreQC_pe_reads/
      which contains the FASTQC reports of the raw paired-end fastq files




3. Raw samples quality control
------------------------------
Quality control analysis of the raw samples can be done using command ``preProcessSamples``

  **Requirements**
  
    1. Execution of prepareProject.py command 
    2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config``.
                                 

    [RNASeq-Analysis]$ rnaseq-wf.py cleanReads <arguments> --local-scheduler
    
    arguments               type      Description
      

    --bbduk-Xms                 int   Initial Java heap size in GB 
                                      Example: 10
                                      Default: 2

    --bbduk-Xmx                 int   Maximum Java heap size in GB
                                      Example: 40
                                      Default: 20

    --bbduk-kmer                int   Kmer length used for finding contaminants
                                      Examle: 13  
                                      Default: 11 

    --bbduk-k-trim              str   Trimming protocol to remove bases matching reference
                                      kmers from reads. Choose From['f: dont trim','r: trim to right','l: trim to left] 
                                      Choices: {f, r, l} 

    --bbduk-quality-trim        int   Trim read ends to remove bases with quality below trimq.
                                      Performed AFTER looking for kmers.  Values: 
                                          rl  (trim both ends), 
                                          f   (neither end), 
                                          r   (right end only), 
                                          l   (left end only),
                                          w   (sliding window).

                                          Default: f


    --bbduk-trim-quality        float Regions with average quality BELOW this will be trimmed,
                                      if qtrim is set to something other than f.  Can be a 
                                      floating-point number like 7.3  
                                      Default: 6                 

    --bbduk-min-length          int   Minimum read length after trimming
                                      Example: 50
                                      Default:50

    --bbduk-trim-front          int   Number of bases to be trimmed from the front of the read
                                      Example: 5
                                      Default: 0

    --bbduk-trim-tail           int   Number of bases to be trimmed from the end of the read
                                      Example: 5
                                      Default: 0

    --bbduk-min-average-quality int   Minimum average quality of reads.
                                      Reads with average quality (after trimming) below 
                                      this will be discarded
                                      Example: 15
                                      Default: 10

    --bbduk-mingc              float Minimum GC content threshold
                                      Discard reads with GC content below minGC
                                      Example: 0.1 
                                      Default: 0.0

    --bbduk-maxgc              float Maximum GC content  threshold
                                      Discard reads with GC content below minGC
                                      Example: 0.99 
                                      Default: 1.0
    --local-scheduler



**Example Run**

      [RNASeq-Analysis]$ python rnaseq.py cleanReads \
                            --bbduk-min-average-quality 15 \
                            --bbduk-mingc 0.20 \
                            --bbduk-maxgc 0.70 \
                            --bbduk-quality-trim w  \
                            --local-scheduler

      **Output**
      /path/to/ProjectFolder/ReadQC/Cleaned_PE_Reads --contains the processed FastQ-reads


**4. Alignment Free Differential Expression Analysis [salmon / kallisto]**
-----------------------------------------------------------------------

4.a. Quantify transcripts


Quantification of the transcripts can be done using command ``alignmentFreeQuant``


  **Requirements**
 1. Pre execution of prepareProject.py command 
 2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst or se_samples.lst`` inside the ``config`` folder.
 Note: If the corresponding gtf/gff file is not provided, pipeline will run prokka to annotate the draft / assembly and generate the transcript file
  

.. code-block:: none   

    [RNASeq-Analysis]$ rnaseq-wf.py alignmentFreeQuant <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the RNASeq reads or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

    --quantMethod            str       Read quantification method
                                      [salmon / kallisto]
    --local-scheduler




  **Example Run 1**
  **quantifyTranscripts** 
  1.  with out read quality control analysis  ``--pre-process-reads no``
  2.  with read quantification method ``salmon``

    [RNASeq-Analysis]$ rnaseq.py  alignmentFreeQuant  --pre-process-reads  ``no`` \
                                   --quantMethod ``salmon`` \
                                   --local-scheduler



  **Example Run 2**
  **quantifyTranscripts** 

  1.  with out read quality control analysis  ``--pre-process-reads no``
  2.  with read quantification method ``kallisto``

    [RNASeq-Analysis]$ rnaseq.py  alignmentFreeQuant  --pre-process-reads  ``no`` \
                                   --quantMethod ``kallisto`` \
                                   --local-scheduler


4.b. Alignment Free Differential Expression Analysis

  Transcript Quantification using ``salmon`` / ``kallisto`` followed by Differential expression analysis with ``DESeq2`` / ``edgeR``


  **Requirements**
  1. Pre execution of prepareProject.py command 
  2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst or se_samples.lst`` inside the ``config`` folder.


    [RNASeq-Analysis]$ rnaseq.py alignmentFreeDEA <arguments> --local-scheduler

    argument               type      Description
    
    [required arguments]  

    --pre-process-reads     str       Run Quality Control Analysis of the RNASeq reads or Not
                                      [yes / no]

                                      If yes, cleanReads command will be run with default parameters.
                                      If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                      script of bbmap will be run based on paired-end or single-end reads.
    --quant-method           str       Read quantification method
                                      [salmon / kallisto]

    --dea-method             str       Method to be used for differential expression analysis. 
                                      [deseq2 / edger]

    --reference-condition    str       Reference biological condition. 
                                      example: control

    --local-scheduler

    [optional arguments] 

    --attribute-type        str       Atrribute type in GTF annotation
                                     choose from {gene_id, transcript_id}

    --strand-type           int       perform strand-specific read counting.
                                     choose from [0: unstranded,
                                                  1: stranded,
                                                  2: reversely-stranded]

    --report-name           str       Name of the differential expression analysis report
                                     Default: DEA_Report

    --factor-of-intrest     str       Factor of intrest column of the target file
                                     Default: conditions

    --fit-type              str       mean-variance relationship. 
                                     Choices: {local, parametric, mean}

    --size-factor           str       method to estimate the size factors. Choices: {shorth, median} 

    --result-tag            str       Tag need to be apended to result file
                                     Default: treated_vs_control
    
    

**Example Run**
**Alignment Free Differential Expression Analysis**

      [RNASeq-Analysis]$ python rnaseq.py alignmentFreeDEA \
                                --pre-process-reads no \
                                --dea-method deseq2 \
                                --reference-condition control \
                                --local-scheduler


5. Alignment based Differential Expression Analysis
----------------------------------------------------

  **Requirements**
  1. Pre execution of prepareProject.py command 
  2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst or se_samples.lst`` inside the ``config`` folder.
  3. Raw Read Quality Assessemnt and Quality Control Analysis as described in step 2 and 3 respectively.


5.a. Index Genome


      [RNASeq-Analysis]$ rnaseq.py indexGenome <arguments> --local-scheduler
      
      argument               type     Description
      [required arguments]      
      
                                       
      --rnaseq-aligner       str      Name of the aligner to map RNASeq reads to genome FASTA file
                                      [bowtie2, dart, subread, star, segemehl, hisat2]
                                       
         

5.b. Generate Gene Counts

     Quantification of the gene can be done using command ``alignmentBasedQuant``


     [required arguments]  
     
                                 
 

     --rnaseq-aligner         str      Name of the aligner to be used to map clean reads to indexed genome. 
                                       Options [star | hisat2 | dart | bowtie2 | segemehl | subread]
                                       NOTE:   star and segemehl demands high memory. 
				       bowtie2 should not be used for domain eukaryote
    
     [optional arguments] 

     --attributeType	      str      Specify attribute type in GTF annotation.
     				       [gene_id / transcript_id]   
                                       Default: gene_id
				 
     --strandType             int     Perform strand-specific read counting.  
                                      [0: unstranded / 1: stranded / 2: reversely-stranded]
                                      Default: 0



5.c. Alignment Based Differential Expression Analysis


Transcript Quantification using `featureCounts` followed by Differential expression analysis with ``DESeq2`` / ``edgeR``


  **Requirements**
  1. Pre execution of prepareProject.py command 
  2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst or se_samples.lst`` inside the ``config`` folder.
  3. Raw Read Quality Assessemnt and Quality Control Analysis as described in step 2 and 3 respectively.



    rnaseq.py alignmentBasedDEA --help

    [required arguments]    
  
    --rnaseq-aligner      	str   Name of the aligner to be used to map clean reads to indexed genome. 
                                      Options [star | hisat2 | dart | bowtie2 | segemehl | subread]
                                      NOTE:   star and segemehl demands high memory. 
                                      bowtie2 should not be used for domain eukaryote


    --pre-process-reads        str    Run Quality Control Analysis of the RNASeq reads or Not
                                      [yes / no]

                                      If yes, cleanReads command will be run with default parameters.
                                      If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                      script of bbmap will be run based on paired-end or single-end reads.
  
    --dea-method		str   Differential Expression Analysis Method to be used
                                      [deseq2 / edger]
   
    --report-name		str   Name of the html report
  
  	  
    --reference-condition       str   Reference biological condition. Default: control


    [optional parameters]
      
    --result-tag 	     str      Tag to be appended to the result folder
                                      Default: treated_vs_control
	
    --alpha                float      Threshold of statistical significance.
	                              Default: 0.05
  
    --p-adjust-method        str      p-value adjustment method. [BH | BY ]
                                      Default [Benjamini & Hochberg(BH)]
					
    --factor-of-intrest      str      Factor of intrest column of the target file 
                                      [=conditions]
  
    --fit-type      	     str      Mean-variance relationship. 
				      [parametric | local | mean]
                                      Default: parametric
  
    --size-factor            str      Method to estimate the size factors. 
				      [median / short]
                                      Default: median


**Example Run**
**Genome based Differential Expression Analysis**


 **Requirements**
 1. Pre execution of prepareProject.py command
 2. Availability of luigi.cfg file in parent folder and samples.txt inside the project folder

        [RNASeq-Analysis]$ rnaseq.py alignmentBasedDEA \
                               --pre-process-reads yes \
                               --rnaseq-aligner hisat2  \
			       --dea-method deseq2  \
                               --reference-condition control  \
                               --local-scheduler


6. Denovo Transcript Assembly Based Differential Expression Analysis
-----------------------------------------------------------------------

6.a. Denovo Transcript Assembly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Denovo assembly of the RNASeq can be done using command ``dnTransAssemble``

  **Requirements**
  1. Pre execution of prepareProject.py command 
  2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
  3. Only paired-end RNASeq reads are used for Denovo Transcript Assembly
  4. Rockhopper Assembler must be used for prokaryotic RNASeq reads and Trinity (or) spades assemblers must be used for eukaryotic RNASeq reads
         
    [RNASeq-Analysis]$ rnaseq.py dnTransAssemble <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the RNASeq reads or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

    --rnaseq-assembler      str       choose from [spades, rockhopper, trinity]
                                     Note: please choose rockhopper if domain is prokaryote
                                           choose trinity / spades if domain is eukaryote
    --local-scheduler

   Note: 
   To use rockhopper as assembler, domain must be prokaryote in luigi.cfg file 
   To use trinity / spades as assembler, domain must be eukaryote in luigi.cfg file 


6.b. Quantify De-novo Assembled Transcripts

Quantification of Denovo assembled transcript can be done using command ``quantifyDAT``


  **Requirements**
 1. Pre execution of prepareProject.py command 
 2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
 3. Only paired-end RNASeq reads are used for Denovo Transcript Assembly
 4. Rockhopper Assembler must be used for prokaryotic RNASeq reads and Trinity (or) spades assemblers must be used for eukaryotic RNASeq reads


        [RNASeq-Analysis]$ rnaseq.py quantifyDAT <arguments> --local-scheduler

        argument               type      Description

        --pre-process-reads    str       Run Quality Control Analysis of the RNASeq reads or Not
                                         [yes / no]

                                         If yes, cleanReads command will be run with default parameters.
                                         If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                         script of bbmap will be run based on paired-end or single-end reads.

        --rnaseq-assembler      str       choose from [spades, rockhopper, trinity]
                                          Note: please choose rockhopper if domain is prokaryote
                                          choose trinity / spades if domain is eukaryote
        --local-scheduler

   Note: 
   To use rockhopper as assembler, domain must be prokaryote in luigi.cfg file 
   To use trinity / spades as assembler, domain must be eukaryote in luigi.cfg file 


  **Example Run 1**
 **quantify De-novo Assembled Transcripts (prokaryotes)** 

  1.  with out read quality control analysis  ``--pre-process-reads no`` and ``--rnaseq-assembler rockhopper``
     Note: To use rockhopper as assembler, domain must be prokaryote in luigi.cfg file 

     [RNASeq-Analysis]$ rnaseq.py  quantifyDAT \
                                ``--pre-process-reads  no`` \
                                ``--rnaseq-assembler rockhopper`` \
                                ``--local-scheduler``




  **Example Run 2**
  **quantify De-novo Assembled Transcripts (eukaryotes)** 

  1.  with out read quality control analysis  ``--pre-process-reads yes`` and ``--rnaseq-assembler trinity``
      Note: To use trinity as assembler, domain must be eukaryote in luigi.cfg file 

     [RNASeq-Analysis]$ rnaseq.py  quantifyDAT \
                                ``--pre-process-reads  yes`` \
                                ``--rnaseq-assembler spades`` \
                                ``--local-scheduler``



6.c. Cluster De-novo Assembled Transcripts
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Cluster of Denovo assembled transcript can be done using command ``clusterDAT``


  **Requirements**
  1. Pre execution of prepareProject.py command 
  2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
  3. Only paired-end RNASeq reads are used for Denovo Transcript Assembly
  4. Rockhopper Assembler must be used for prokaryotic RNASeq reads and Trinity (or) spades assemblers must be used for eukaryotic RNASeq reads
 


.. code-block:: none   

    [RNASeq-Analysis]$ rnaseq.py clusterDAT <arguments> --local-scheduler

    argument               type      Description

    --pre-process-reads    str       Run Quality Control Analysis of the RNASeq reads or Not
                                     [yes / no]

                                     If yes, cleanReads command will be run with default parameters.
                                     If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                     script of bbmap will be run based on paired-end or single-end reads.

    --rnaseq-assembler      str      choose from [spades, rockhopper, trinity]
                                     Note: please choose rockhopper if domain is prokaryote
                                           choose trinity / spades if domain is eukaryote
    --local-scheduler

   Note: 
   To use rockhopper as assembler, domain must be prokaryote in luigi.cfg file 
   To use trinity / spades as assembler, domain must be eukaryote in luigi.cfg file 


  **Example Run 1**
  **cluster De-novo Assembled Transcripts (prokaryotes)** 

  1.  with out read quality control analysis  ``--pre-process-reads no`` and ``--rnaseq-assembler rockhopper``
      Note: To use rockhopper as assembler, domain must be prokaryote in luigi.cfg file 

    [RNASeq-Analysis]$ rnaseq.py  clusterDAT \
                                ``--pre-process-reads  no`` \
                                ``--rnaseq-assembler rockhopper`` \
                                ``--local-scheduler``

  **Example Run 2**
  **cluster De-novo Assembled Transcripts (eukaryotes)** 

    1.  with out read quality control analysis  ``--pre-process-reads yes`` and ``--rnaseq-assembler trinity``
        Note: To use trinity as assembler, domain must be eukaryote in luigi.cfg file 

        [RNASeq-Analysis]$ rnaseq.py  clusterDAT \
                                ``--pre-process-reads  yes`` \
                                ``--rnaseq-assembler spades`` \
                                ``--local-scheduler``



6.d. Denovo transcript assembly Based Differential Expression Analysis

    Transcript Quantification using ``salmon`` Clustering using ``Corset`` followed by Differential expression analysis with ``DESeq2`` / ``edgeR``


  **Requirements**
  1. Pre execution of prepareProject.py command 
  2. Availability of ``luigi.cfg`` file in ``parent folder`` and ``pe_samples.lst`` inside the ``config`` folder.
  3. Raw Read Quality Assessemnt and Quality Control Analysis as described in step 2 and 3 respectively.
 

    rnaseq.py denovoDEA --help

    [required arguments]    
  
    --rnaseq-assembler       str     Choose from [spades, rockhopper, trinity]
                                     Note: please choose rockhopper if domain is prokaryote
                                           choose trinity / spades if domain is eukaryote

    --pre-process-reads      str      Run Quality Control Analysis of the RNASeq reads or Not
                                      [yes / no]

                                      If yes, cleanReads command will be run with default parameters.
                                      If no, quality control analysis will not be done, instead re-pair.sh or reformat.sh 
                                      script of bbmap will be run based on paired-end or single-end reads.
  
    --dea-method             str      Differential Expression Analysis Method to be used
                                      [deseq2 / edger]
   
     
    --reference-condition    str      Reference biological condition. Default: control
                                      Note: please check ``target.tsv`` file inside ``config`` folder for `` --reference-condition``

    [optional parameters]
      
    --result-tag             str      Tag to be appended to the result folder
                                      Default: treated_vs_control
  
    --alpha                  float    Threshold of statistical significance.
                                      Default: 0.05
  
    --p-adjust-method        str      p-value adjustment method. [BH | BY ]
                                      Default [Benjamini & Hochberg(BH)]
          
    --factor-of-intrest      str      Factor of intrest column of the target file 
                                      [=conditions]
  
    --fit-type               str      Mean-variance relationship. 
                                      [parametric | local | mean]
                                      Default: parametric
  
    --size-factor            str      Method to estimate the size factors. 
                                      [median / short]
                                      Default: median
