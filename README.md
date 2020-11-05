Installation Instructions
================================
  
To install RNASeq Analysis Toolkit, you must have a minimum of 6 GiB free disk space and minimum of 16 GiB free RAM to test run. 

To provide an easier way to install, we provide a miniconda based installer.
Installation also requires **pre-instaled** ``git``, ``gcc``, ``cpp`` and ``zlib1g-dev``.
  
    git clone https://github.com/computational-genomics-lab/bulkRNASeqPIPE.git
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

 Required Files

 1. The reference genome or assembly of the organism in FASTA format with ``.fna , .fa or .fasta`` extention
 2. Gene Annotation files in GFF3/GTF format 
    Note: GFF3/GTF file is optional if domain is prokaryote and differential expression analysis method is Alignment Free (Salmon/Kallisto)
          GFF3/GTF file is mandatory if differential expression analysis method is Alignment Based (STAR, hisat2, dart, bowtie2, subread)

 3. The transcriptome of the organism in FASTA format with ``.ffn , .fa or .fasta`` extension (optional for both Prokaryotes and eukaryotes)
   
For prokaryotes, if user provides the reference genome/assembly of the organism in FASTA format, the pipeline will trigger    prokka to annotate the reference      genome and the corresponding transcriptome, gff, gtf and tx2gene files will be made available for downstream analysis.

For eukaryotes, user need to provide the gene annotation file in GTF/GFF3 format along with the corresponding reference genome   or genome assembly. The pipeline  will use gffread program to extract transcript form the genome using the GFF3/GTF file, which will be used for Alignment Free Differential Expression Analysis.

Reference / draft genome of the organism along with gene annotation files (optional for prokaryotes mandatory for eukaryotes)  transcriptome of the organism in FASTA format with ``.ffn , .fa or .fasta`` extension (optional for prokaryotes eukaryotes) must be made available in the same folder containing RNASeq reads.

 Allowed ``externsion`` for 
     a. genome FASTA file: ``.fna, .fa , .fasta``
     b. transcriptome FASTA file: ``.ffn, .fa , .fasta``
     c. genome annotation file: ``.gff`` or ``.gtf``


