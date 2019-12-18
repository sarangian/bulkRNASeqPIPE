
.. _folderstructure:


bulkRNASeqPipe folder structure
===============================

.. code-block:: none

    bulkRNASeqPipe
    ├── raw_data     [contains three sub-folders i.e. genome, RNASeq and transcriptome]
    │   ├── genome
    │   │   ├── arthrinium
    │   │   ├── hexcentricum
    │   │   └── yeast
    │   ├── RNASeq   [contains the raw RNASeq Reads]
    │   └── transcriptome
    │       ├── arthrinium
    │       └── yeast
    └── utility      [contains R Scripts for Differential Expression Analysis]
        └── deaRscripts
            ├── man
            └── R


bulkRNASeqPipe sub-folder structure
===================================

.. code-block:: none

    bulkRNASeqPipe
    ├── arthrinium_target.txt
    ├── RNASeqPIPE.py
    ├── commands.txt
    ├── halomicronema_target.txt
    ├── INSTALL.sh
    ├── LICENSE
    ├── raw_data  
    │   ├── genome                           
    │   │   ├── arthrinium
    │   │   │   ├── arthrinium.fna.gz
    │   │   │   └── arthrinium.gtf
    │   │   ├── hexcentricum
    │   │   │   ├── hexcentricum.fna
    │   │   │   └── hexcentricum.gtf
    │   │   └── yeast
    │   │       ├── yeast.fna
    │   │       └── yeast.gtf
    │   ├── RNASeq
    │   │   ├── arthrinium_s1_rep1_R1.fastq.gz
    │   │   ├── arthrinium_s1_rep1_R2.fastq.gz
    │   │   ├── arthrinium_s1_rep2_R1.fastq.gz
    │   │   ├── arthrinium_s1_rep2_R2.fastq.gz
    │   │   ├── arthrinium_s2_rep1_R1.fastq.gz
    │   │   ├── arthrinium_s2_rep1_R2.fastq.gz
    │   │   ├── arthrinium_s2_rep2_R1.fastq.gz
    │   │   ├── arthrinium_s2_rep2_R2.fastq.gz
    │   │   ├── hexcentricum_s1_rep1_R1.fastq.gz
    │   │   ├── hexcentricum_s1_rep1_R2.fastq.gz
    │   │   ├── hexcentricum_s1_rep2_R1.fastq.gz
    │   │   ├── hexcentricum_s1_rep2_R2.fastq.gz
    │   │   ├── hexcentricum_s2_rep1_R1.fastq.gz
    │   │   ├── hexcentricum_s2_rep1_R2.fastq.gz
    │   │   ├── hexcentricum_s2_rep2_R1.fastq.gz
    │   │   ├── hexcentricum_s2_rep2_R2.fastq.gz
    │   │   ├── yeast_sample_1_replicate_1.fastq.gz
    │   │   ├── yeast_sample_1_replicate_2.fastq.gz
    │   │   ├── yeast_sample_2_replicate_1.fastq.gz
    │   │   └── yeast_sample_2_replicate_2.fastq.gz
    │   └── transcriptome
    │       ├── arthrinium
    │       │   ├── arthrinium.ffn
    │       │   └── arthrinium.gtf
    │       └── yeast
    │           ├── yeast.ffn
    │           └── yeast.gtf
    ├── README.md
    ├── run_arthrinium.sh
    └── utility
        ├── adapters.fasta.gz
        ├── deaRscripts
        │   ├── check.Target.R
        │   ├── checkTarget.R
        │   ├── corset_DESeq2.r
        │   ├── corset_edgeR.r
        │   ├── dea_tx2gene_from_gff3.R
        │   ├── DESCRIPTION
        │   ├── featureCount_DESeq2.r
        │   ├── featureCount_edgeR.r
        │   ├── kallisto_DESeq2.r
        │   ├── kallisto_edgeR.r
        │   ├── man
        │   │   ├── exportResults.DESeq2.Rd
        │   │   ├── exportResults.edgeR.Rd
        │   │   ├── loadTargetFile.Rd
        │   │   ├── nDiffTotal.Rd
        │   │   ├── run.DESeq2_corset.Rd
        │   │   ├── run.DESeq2_trans.Rd
        │   │   ├── run.edgeR.Rd
        │   │   ├── run.edgeR_trans.Rd
        │   │   ├── summarizeResults.DESeq2.Rd
        │   │   ├── summarizeResults.edgeR.Rd
        │   │   ├── tr2g_EnsDb.Rd
        │   │   ├── tr2g_ensembl.Rd
        │   │   ├── tr2g_fasta.Rd
        │   │   ├── tr2g_gff3.Rd
        │   │   ├── tr2g_GRanges.Rd
        │   │   ├── tr2g_gtf.Rd
        │   │   ├── tr2g_junction.Rd
        │   │   └── tr2g_TxDb.Rd
        │   ├── PlotDESEQ2.Rmd
        │   ├── PlotEDGER.Rmd
        │   ├── R
        │   │   ├── exportResults.DESeq2.R
        │   │   ├── exportResults.edgeR.R
        │   │   ├── load.TargetFile.R
        │   │   ├── nDiffTotal.r
        │   │   ├── run.DESeq2_corset.r
        │   │   ├── run.DESeq2_Salmon.r
        │   │   ├── run.DESeq2_trans.r
        │   │   ├── run.edgeR.r
        │   │   ├── run.edgeR_trans.r
        │   │   ├── summarizeResults.DESeq2.r
        │   │   ├── summarizeResults.edgeR.r
        │   │   ├── tr2g.R
        │   │   └── utils.R
        │   ├── README.md
        │   ├── rnaseqflow.Rproj
        │   ├── salmon_DESeq2.r
        │   ├── salmon_edgeR.r
        │   └── tx2gene.R
        ├── FormatRockhopTrans.py
        ├── __init__.py
        ├── PlotDESEQ2.Rmd
        ├── PlotEDGER.Rmd
        └── utility.py