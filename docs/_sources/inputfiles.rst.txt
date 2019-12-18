
.. _inputfiles:

Input Files
===========

Raw RNASeq Reads
----------------

| The raw-RNASeq reads must be placed inside raw_data/RNASeq folder.
| Reads must be gzipped.
| 
| For **paired-end** RNAseq reads, sample name must be suffixed with ``_R1.fastq.gz`` and ``_R2.fastq.gz``. 
|          *Example*
|           hexcentricum_s1_rep1_R1.fastq.gz 
|           hexcentricum_s1_rep1_R2.fastq.gz
|           
|  where  ``hexcentricum_s1_rep1`` is sample name
|         ``_R1.fastq.gz`` and ``_R2.fastq.gz`` are sufixes to sample name
|      
| For **single-end** RNAseq reads sample name must be suffixed with ``.fastq.gz``
|          *Example*
|           yeast_sample_1_replicate_1.fastq.gz
| 
|  where  ``yeast_sample_1_replicate_1`` is sample name
|         ``.fastq.gz`` is sufixes to sample name
|
|


sampleMetadataFile or targetFile
--------------------------------

The user has to provide the sampleMetadataFile --a ``tab delimited`` file which describes the experiment. This file has one row per sample and is composed of at least three columns with headers. The sampleMetadataFile must be present in the **bulkRNASeqPipe** folder.

**Example:**  halomicronema_target.txt

====================  ====================  ==========
lable                 samples               conditions
====================  ====================  ==========
hexcentricum_s1_rep1  hexcentricum_s1_rep1  Control 
hexcentricum_s1_rep2  hexcentricum_s1_rep2  Control
hexcentricum_s2_rep1  hexcentricum_s2_rep1  Treated
hexcentricum_s2_rep2  hexcentricum_s2_rep2  Treated    
====================  ====================  ==========

|   Column 1 : sample labels --must be unique. Header of column 1 is ``lable``
|   Column 2 : name of the samples --must be unique. Header of column 2 must be ``samples``
|   Column 3 : biological conditions associated with the samples. Header of column 3 is ``conditions``
|   Column 3 : reference condition name is ``Control``
| 

|   NOTE: *Each sample must have at-least two replicates*
|

Genome FASTA file
-----------------

The reference genome or assembly of the organism in FASTA format with ``.fna`` extension is required for 

| 1. Genome Guided Transcript Assembly Pipeline
|
| 2. Genome Based (Prokaryote and Eukaryote) Differential Expression Analysis Pipeline
|
| 3. Prokaryotic Transcriptome based Differential Expression Analysis Pipeline (Transcriptome will be predicted from genome using PROKKA)
|
| **NOTE**
| User must create a sub-folder, **named the same** as the **name of the organism** inside **raw_data/genome** folder. Sub-folder name must not contain any gap or special characters. The genome of organism in FASTA format with ``.fna`` extension must made available inside the Sub-folder.
|
| **Example:**
| The genome of *Yeast* in FASTA format (yeast.fna) is available inside the folder **raw_data/genome/yeast**
| The genome of *Arthrinium*  in FASTA format (arthrinium.fna) is available inside the folder **raw_data/genome/arthrinium**
| 

.. code-block:: none

        genome                           
           ├── arthrinium          
           │     ├── arthrinium.fna
           │     └── arthrinium.gtf
           ├── hexcentricum
           │     ├── hexcentricum.fna
           │     └── hexcentricum.gtf
           └── yeast              
                   ├── yeast.fna
                   └── yeast.gtf


Transcriptome FASTA file
------------------------

| The reference transcriptome (cDNA) or assembled transcriptome of the organism in FASTA format with ``.ffn`` extension is mandatory for Eukaryotic Transcriptome Based Differential Expression Analysis Pipeline.
| 
| For Prokaryotes, transcriptome will be predicted from the genome (if provided )using PROKKA and the predicted transcriptome will be used for downstream analysis.
|
| User must create a sub-folder, **named the same** as the **name of the organism** inside **raw_data/transcriptome** folder. Sub-folder name must not contain any gap or special characters. The transcriptome of organism in FASTA format with ``.ffn`` extension must made available inside the Sub-folder.
|

|  **Example:**
|

.. code-block:: none

            transcriptome
                       ├── arthrinium
                       │   ├── arthrinium.ffn
                       │   └── arthrinium.gtf
                       └── yeast
                           ├── yeast.ffn
                           └── yeast.gtf

|
| The transcriptome of *Yeast* in FASTA format (yeast.ffn) is available inside the folder **raw_data/transcriptome/yeast**
| The transcriptome of *Arthrinium*  in FASTA format (arthrinium.fna) is available inside the folder **raw_data/transcriptome/arthrinium**
|

Genome Annotation File
----------------------
| 
| The genome annotation file of the organism in Gene transfer format (GTF) must made available along with the corresponnding genome (.fna) and (or) transcriptome (.ffn) FASTA file of the organism.
|

| **Example**
|
| The annotation of *Yeast* in GTF format (yeast.gtf) is available inside the folder **raw_data/genome/yeast** and **raw_data/transcriptome/yeast**
|
| The annotation of *Arthrinium* in GTF format (arthrinium.gtf) is available inside the folder **raw_data/genome/arthrinium** and **raw_data/transcriptome/arthrinium**
