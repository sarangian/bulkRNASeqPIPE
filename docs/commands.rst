
.. _commands:

Commands
========

.. toctree::
   :hidden:


Major Commnads
--------------


Usage:  python RNASeqPIPE.py <command> <arguments>

.. code-block:: none

    Command                      Description   
    
    prepareProject               Creates a project folder for storing the analysis results  
    rawReadsQC                   Raw Reads Quality Assessment using FASTQC tool  
    preProcessSamples            Process Raw Reads using BBDUK
    denovoTransAssemble          Denovo Assembly of Prokaryotic and Eukaryotic Transcripts
    quantifyAssembledTranscripts Quantify Assembled transcripts using Salmon 
    clusterContigs               Clustred Assembled transcripts based on equivalence class       

