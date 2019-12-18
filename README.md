https://bulkrnaseqpipe.readthedocs.io/en/latest/inputfiles.html

# bulkRNASeqPIPE
#### Differential Expression Analysis of bulk-RNASeq data


## Dependencies

**Installation of bulkRNASeqPIPE requires pre-instaled git, gcc, cpp and zlib1g-dev.**


## Download and install
```
git clone https://github.com/sarangian/bulkRNASeqPIPE.git
cd bulkRNASeqPIPE
chmod 755 INSTALL.sh
./INSTALL.sh
```
## Post installation instruction

After successful installation, close the current terminal. 
In a new terminal `source ~/.bashrc`
Activate bulkRNASeqPIPE environment using command `conda activate`

## bulkRNASeqPIPE folder structure
```
..
├── arthrinium_target.txt
├── RNASeqPIPE.py
├── commands.txt
├── halomicronema_target.txt
├── INSTALL.sh
├── LICENSE
├── miniconda_x64_installer.sh
├── raw_data
│   ├── genome
│   │   ├── arthrinium
│   │   │   ├── arthrinium.fna.gz
│   │   │   └── arthrinium.gtf
│   │   ├── hexcentricum
│   │   │   ├── hexcentricum.fna
│   │   │   └── hexcentricum.gtf
│   │   └── yeast
│   │       ├── yeast.fna
│   │       └── yeast.gtf
│   ├── RNASeq
│   │   ├── arthrinium_s1_rep1_R1.fastq.gz
│   │   ├── arthrinium_s1_rep1_R2.fastq.gz
│   │   ├── arthrinium_s1_rep2_R1.fastq.gz
│   │   ├── arthrinium_s1_rep2_R2.fastq.gz
│   │   ├── arthrinium_s2_rep1_R1.fastq.gz
│   │   ├── arthrinium_s2_rep1_R2.fastq.gz
│   │   ├── arthrinium_s2_rep2_R1.fastq.gz
│   │   ├── arthrinium_s2_rep2_R2.fastq.gz
│   │   ├── hexcentricum_s1_rep1_R1.fastq.gz
│   │   ├── hexcentricum_s1_rep1_R2.fastq.gz
│   │   ├── hexcentricum_s1_rep2_R1.fastq.gz
│   │   ├── hexcentricum_s1_rep2_R2.fastq.gz
│   │   ├── hexcentricum_s2_rep1_R1.fastq.gz
│   │   ├── hexcentricum_s2_rep1_R2.fastq.gz
│   │   ├── hexcentricum_s2_rep2_R1.fastq.gz
│   │   ├── hexcentricum_s2_rep2_R2.fastq.gz
│   │   ├── yeast_sample_1_replicate_1.fastq.gz
│   │   ├── yeast_sample_1_replicate_2.fastq.gz
│   │   ├── yeast_sample_2_replicate_1.fastq.gz
│   │   └── yeast_sample_2_replicate_2.fastq.gz
│   └── transcriptome
│       ├── arthrinium
│       │   ├── arthrinium.ffn
│       │   └── arthrinium.gtf
│       └── yeast
│           ├── yeast.ffn
│           └── yeast.gtf
├── README.md
├── run_arthrinium.sh
└── utility
    ├── adapters.fasta.gz
    ├── deaRscripts
    │   ├── check.Target.R
    │   ├── checkTarget.R
    │   ├── corset_DESeq2.r
    │   ├── corset_edgeR.r
    │   ├── dea_tx2gene_from_gff3.R
    │   ├── DESCRIPTION
    │   ├── featureCount_DESeq2.r
    │   ├── featureCount_edgeR.r
    │   ├── kallisto_DESeq2.r
    │   ├── kallisto_edgeR.r
    │   ├── man
    │   │   ├── exportResults.DESeq2.Rd
    │   │   ├── exportResults.edgeR.Rd
    │   │   ├── loadTargetFile.Rd
    │   │   ├── nDiffTotal.Rd
    │   │   ├── run.DESeq2_corset.Rd
    │   │   ├── run.DESeq2_trans.Rd
    │   │   ├── run.edgeR.Rd
    │   │   ├── run.edgeR_trans.Rd
    │   │   ├── summarizeResults.DESeq2.Rd
    │   │   ├── summarizeResults.edgeR.Rd
    │   │   ├── tr2g_EnsDb.Rd
    │   │   ├── tr2g_ensembl.Rd
    │   │   ├── tr2g_fasta.Rd
    │   │   ├── tr2g_gff3.Rd
    │   │   ├── tr2g_GRanges.Rd
    │   │   ├── tr2g_gtf.Rd
    │   │   ├── tr2g_junction.Rd
    │   │   └── tr2g_TxDb.Rd
    │   ├── PlotDESEQ2.Rmd
    │   ├── PlotEDGER.Rmd
    │   ├── R
    │   │   ├── exportResults.DESeq2.R
    │   │   ├── exportResults.edgeR.R
    │   │   ├── load.TargetFile.R
    │   │   ├── nDiffTotal.r
    │   │   ├── run.DESeq2_corset.r
    │   │   ├── run.DESeq2_Salmon.r
    │   │   ├── run.DESeq2_trans.r
    │   │   ├── run.edgeR.r
    │   │   ├── run.edgeR_trans.r
    │   │   ├── summarizeResults.DESeq2.r
    │   │   ├── summarizeResults.edgeR.r
    │   │   ├── tr2g.R
    │   │   └── utils.R
    │   ├── README.md
    │   ├── rnaseqflow.Rproj
    │   ├── salmon_DESeq2.r
    │   ├── salmon_edgeR.r
    │   └── tx2gene.R
    ├── filter_rcorrector_pe.py
    ├── filter_rcorrector_se.py
    ├── FormatRockhopTrans.py
    ├── __init__.py
    ├── PlotDESEQ2.Rmd
    ├── PlotEDGER.Rmd
    └── utility.py

```
## Input files
###### bulkRNASeq pipeline script: `RNASeqPIPE.py`
###### --sampleMetadataFile `halomicronema_target.txt`

The user has to provide the `sampleMetadataFile` --a tab delimited file which describes the experiment. 
This file has one row per sample and is composed of at least three columns with headers. 
- **column 1 :** sample labels --must be unique. Header of column 1 is `lable`
- **column 2 :** name of the samples --must be unique. Header of column 2 **must be** `samples`
- **column 3 :** biological conditions associated with the samples. Header of column 3 is `conditions`
- **column 3 :** reference condition name is 'Control'
- **NOTE:** each sample must have at-least two replicates

**Example: sampleMetadataFile: halomicronema_target.txt**
```
lable	                  samples                 conditions
hexcentricum_s1_rep1	  hexcentricum_s1_rep1	  Control
hexcentricum_s1_rep2	  hexcentricum_s1_rep2	  Control
hexcentricum_s2_rep1	  hexcentricum_s2_rep1	  Treated
hexcentricum_s2_rep2	  hexcentricum_s2_rep2	  Treated
```

**Example: sampleMetadataFile: arthrinium_target.txt**
```
Lable	                  samples	          group
arthrinium_s1_rep1	  arthrinium_s1_rep1	  control
arthrinium_s1_rep2	  arthrinium_s1_rep2	  control
arthrinium_s2_rep1	  arthrinium_s2_rep1	  treated
arthrinium_s2_rep2	  arthrinium_s2_rep2	  treated
```

## Folder Structure
**raw_data:** The folder raw_data contains three sub-folders [1] RNASeq [2] genome and [3] transcriptome.

**raw_data/RNASeq** folder contains the raw-RNASeq reads. 
```
Reads must be gzipped. 
- For paired-end RNAseq reads, sample name must be suffixed with _R1.fastq.gz and _R2.fastq.gz. 
  Example: hexcentricum_s1_rep1_R1.fastq.gz 
           hexcentricum_s1_rep1_R2.fastq.gz
           
           where: hexcentricum_s1_rep1 is sample name
                  _R1.fastq.gz and _R2.fastq.gz are sufixes to sample name
       
- For single-end RNAseq reads sample name must be suffixed with .fastq.gz
  Example: yeast_sample_1_replicate_1.fastq.gz
  
           where: yeast_sample_1_replicate_1 is sample name
                  .fastq.gz is sufixes to sample name
```
**raw_data/genome** 
```
The sub-folder hexcentricum contain 
[1] The genome of the organism in fasta format with .fna extension (hexcentricum.fna) [MANDATORY for both prokaryotes or eukaryotes]
[2] The gene annotation of the organism in GTF format (hexcentricum.gtf) [MANDATORY for eukaryotes, OPTIONAL for prokaryotes]
NOTE: The sub-folder name, genome name and the GTF file name must be same

EXAMPLE: If the organism name is XYZ, then the sub-folder name must be XYZ, genome name must be XYZ.fna, annotation file name must be XYZ.gtf
```
**raw_data/transcriptome** 
```
The sub-folder hexcentricum contain 
[1] The transcriptome of the organism in fasta format with .ffn extension (hexcentricum.ffn) [MANDATORY for both prokaryotes or eukaryotes]
[2] The gene annotation of the organism in GTF format (hexcentricum.gtf) [MANDATORY for eukaryotes, OPTIONAL for prokaryotes]
NOTE: The sub-folder name, genome name and the GTF file name must be same

EXAMPLE: If the organism name is XYZ, then the sub-folder name must be XYZ, genome name must be XYZ.fna, annotation file name must be XYZ.g

## Test Installation

###### ARTHRINIUM

###### Prepare Project
python RNASeqPIPE.py prepareProject --projectName arthrinium --sampleMetadataFile arthrinium_target.txt --sampleHeader samples --groupHeader group --refCondName control --local-scheduler

###### Genome Guided Assembly and DE Analysis
python RNASeqPIPE.py genomeGuidedDEA --aligner hisat2 --projectName arthrinium --deaMethod deseq2 --resultTag treated_vs_control --readType paired --reportName GGTrinity-DESeq2 --factorInt group --refCond control --domain eukaryote --targetFile arthrinium_target.txt --genomeName arthrinium --local-scheduler 2>&1 | tee arthrinium_genome_guided_deseq2.log

###### Transcript based differential expression analysis
python RNASeqPIPE.py transcriptomeBasedDEA --projectName arthrinium --domain eukaryote --predTranscript no --readType paired --quantMethod salmon --genomeName arthrinium --transcriptName arthrinium --deaMethod deseq2 --factorInt group --refCond control --resultTag arthrinium_treated_vs_control --targetFile arthrinium_target.txt --local-scheduler 2>&1 | tee arthrinium_transcriptome_salmon_based_DEA_DESeq2.log

###### Genome based differential expression analysis
python RNASeqPIPE.py genomeBasedDEA --projectName arthrinium --domain eukaryote --aligner hisat2 --readType paired --genomeName arthrinium --deaMethod deseq2 --resultTag arthrinium_treated_vs_control --factorInt group --refCond control --targetFile arthrinium_target.txt --local-scheduler 2>&1 | tee arthrinium_genome_based_DEA_DESeq2.log


###### HALOMICRONEMA

###### Prepare Project
python RNASeqPIPE.py prepareProject --projectName `halomicronema` --sampleMetadataFile `halomicronema_target.txt` --sampleHeader `Samples` --groupHeader `Conditions` --refCondName `Control` --local-scheduler

###### Transcript based differential expression analysis
python RNASeqPIPE.py transcriptomeBasedDEA --projectName halomicronema/ --domain prokaryote --predTranscript yes --readType paired --quantMethod kallisto --genomeName hexcentricum --transcriptName hexcentricum --deaMethod deseq2 --factorInt Conditions --refCond Control --resultTag hexcentricum_treated_vs_control --targetFile halomicronema_target.txt --local-scheduler
