
.. _install:

Installation
================================

.. toctree::
   :hidden:
  
To install bulkRNASeqPIPE, you must have a minimum of 6 GiB free disk space and minimum of 16 GiB free RAM to test run. 
bulkRNASeqPIPE has a lot of dependencies and therefore installation is the most difficult part of executing the pipeline. 
To provide an easier way to install, we provide a miniconda based installer.
Installation also requires **pre-instaled** ``git``, ``gcc``, ``cpp`` and ``zlib1g-dev``.

.. code-block:: none
    
    git clone https://github.com/sarangian/bulkRNASeqPIPE.git/
    cd bulkRNASeqPIPE
    chmod 755 INSTALL.sh
    /INSTALL.sh 2>&1 | tee install.log

    
Post Installation Instruction:

.. code-block:: none
	
    After successful installation, close the current terminal. 
    In a new terminal. source the bashrc file:  source ~/.bashrc
    Activate bulkRNASeqPIPE environment using command: conda activate 


Check Installation
^^^^^^^^^^^^^^^^^^^
**Prepare a Project for Arthrinium RNASeq Analysis**

.. code-block:: none
    
    python RNASeqPIPE.py prepareProject \
                         --projectName arthrinium \
                         --sampleMetadataFile arthrinium_target.txt \
                         --sampleHeader samples \
                         --groupHeader group \
                         --refCondName control \
                         --local-scheduler

**Denovo Transcriptome Assembly**

.. code-block:: none

    python RNASeqPIPE.py denovoTransAssemble \
                         --projectName arthrinium \
                         --domain eukaryote \
                         --readType paired \
                         --assembler trinity \
                         --local-scheduler 2>&1 | tee arthrinium_deovo_trans_assembly_trinity.log

**Denovo Transcriptome Assembly and differential expression analysis**

.. code-block:: none

    python RNASeqPIPE.py denovoDEA \
                         --projectName arthrinium \
                         --assembler trinity \
                         --deaMethod deseq2 \
                         --resultTag treated_vs_control \
                         --readType paired \
                         --reportName denovo_trinity_DESeq2 \
                         --factorInt group --refCond control \
                         --domain eukaryote \
                         --targetFile arthrinium_target.txt \
                         --local-scheduler 2>&1 | tee arthrinium_deovo_DEA_trinity.log

**Genome Guided Transcriptome Assembly and differential expression analysis**

.. code-block:: none

    python RNASeqPIPE.py genomeGuidedDEA \
                        --aligner hisat2 \
                        --projectName arthrinium \
                        --deaMethod deseq2 \
                        --resultTag treated_vs_control \
                        --readType paired \
                        --reportName GGTrinity-DESeq2 \
                        --factorInt group \
                        --refCond control \
                        --domain eukaryote \
                        --targetFile arthrinium_target.txt \
                        --genomeName arthrinium \
                        --local-scheduler 2>&1 | tee arthrinium_genome_guided_deseq2.log

**Transcript based differential expression analysis**

.. code-block:: none

    python RNASeqPIPE.py transcriptomeBasedDEA \
                        --projectName arthrinium \
                        --domain eukaryote \
                        --predTranscript no \
                        --readType paired \
                        --quantMethod salmon \
                        --genomeName arthrinium \
                        --transcriptName arthrinium \
                        --deaMethod deseq2 \
                        --factorInt group \
                        --refCond control \
                        --resultTag arthrinium_treated_vs_control \
                        --targetFile arthrinium_target.txt \
                        --local-scheduler 2>&1 | tee arthrinium_transcriptome_salmon_based_DEA_DESeq2.log

**Genome based differential expression analysis**

.. code-block:: none

    python RNASeqPIPE.py genomeBasedDEA \
                        --projectName arthrinium \
                        --domain eukaryote \
                        --aligner hisat2 \
                        --readType paired \
                        --genomeName arthrinium \
                        --deaMethod deseq2 \
                        --resultTag arthrinium_treated_vs_control \
                        --factorInt group \
                        --refCond control \
                        --targetFile arthrinium_target.txt \
                        --local-scheduler 2>&1 | tee arthrinium_genome_based_DEA_DESeq2.log
