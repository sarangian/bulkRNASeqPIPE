
.. _tutorial:

Tutorial
========

.. toctree::
   :hidden:


Prokaryotic Transcriptome based DE Analysis
-------------------------------------------


Prepare a Project for Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none
    
    python RNASeqPIPE.py prepareProject \
    	   --projectName halomicronema \
    	   --sampleMetadataFile halomicronema_target.txt \
    	   --sampleHeader samples \
    	   --groupHeader conditions \
    	   --refCondName Control \
    	   --local-scheduler


Raw Read Quality Assessment
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    python RNASeqPIPE.py rawReadsQC \
           --projectName halomicronema \
           --readType paired \
           --local-scheduler


Raw Read Quality Control Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    python RNASeqPIPE.py preProcessSamples \
           --projectName halomicronema \
           --readType paired \
           --bbduk-cpu 16 \
           --bbduk-minGC 0.2 \
           --bbduk-maxGC 0.7 \
           --bbduk-minAQ 20 \
           --local-scheduler


Quantify Transcripts
^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    python RNASeqPIPE.py quantifyTranscripts \
           --projectName halomicronema \
           --readType paired \
           --quantMethod kallisto \
           --domain prokaryote \
           --genomeName hexcentricum \
           --transcriptName hexcentricum \
           --predTranscript yes --local-scheduler


Transcriptome Based Differential Expression Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    python RNASeqPIPE.py transcriptomeBasedDEA \
           --projectName halomicronema \
           --domain prokaryote \
           --predTranscript yes \
           --readType paired \
           --quantMethod kallisto \
           --genomeName hexcentricum \
           --transcriptName hexcentricum \
           --deaMethod deseq2 \
           --factorInt conditions \
           --refCond Control \
           --resultTag hexcentricum_treated_vs_control \
           --targetFile halomicronema_target.txt \
           --local-scheduler



Eukaryotic Transcriptome based DE Analysis
-------------------------------------------






