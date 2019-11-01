# bulkRNASeqPIPE
Differential Expression Analysis of bulk-RNASeq data


## Dependencies

Installation of bulkRNASeqPIPE, requires libboost-all-dev, libbz2-dev, and liblzma-dev installed in your system.


## Download and install

git clone https://github.com/sarangian/bulkRNASeqPIPE.git

cd bulkRNASeqPIPE

chmod 755 miniconda_x64_bulkRNASeq_install.sh

./miniconda_x64_bulkRNASeq_install.sh


## Post installation instruction

After successful installation. close the current terminal. In a new terminal please `source ~/.bashrc`
and activate bulkRNASeqPIPE environment using command `conda activate`

## Input files
###### --sampleMetadataFile `halomicronema_target.txt`

The user has to provide the `sampleMetadataFile` --a tab delimited file which describes the experiment. 
This file has one row per sample and is composed of at least three columns with headers. 
- column 1 : sample labels --must be unique. Header of column 1 is `Lable`
- column 2 : name of the samples --must be unique. Header of column 2 is `samples`
- column 3 : biological conditions associated with the samples. Header of column 3 is `group`
- column 3 : reference condition name is 'Control'
- NOTE: each sample must have two replicates
```
Lable	                  samples                 group
sample_1_replicate_1	  sample_1_replicate_1	  Control
sample_1_replicate_2	  sample_1_replicate_2	  Control
sample_2_replicate_1	  sample_2_replicate_1	  Treated
sample_2_replicate_2	  sample_2_replicate_2	  Treated
```
```
.
└── bulkRNASeqPIPE
    ├── bulkRNASeqWF.py
    ├── commands.txt
    ├── halomicronema_target.txt
    ├── LICENSE
    ├── miniconda_x64_bulkRNASeq_install.sh
    ├── raw_data
    │   ├── genome
    │   │   └── hexcentricum
    │   │       ├── hexcentricum.fna
    │   │       └── hexcentricum.gtf
    │   └── RNASeq
    │       ├── sample_1_replicate_1_R1.fastq.gz
    │       ├── sample_1_replicate_1_R2.fastq.gz
    │       ├── sample_1_replicate_2_R1.fastq.gz
    │       ├── sample_1_replicate_2_R2.fastq.gz
    │       ├── sample_2_replicate_1_R1.fastq.gz
    │       ├── sample_2_replicate_1_R2.fastq.gz
    │       ├── sample_2_replicate_2_R1.fastq.gz
    │       └── sample_2_replicate_2_R2.fastq.gz
    ├── README.md
    ├── setup.py
    ├── ubuntu_install.sh
    └── utility
        ├── adapters.fasta.gz
        ├── filter_rcorrector_pe.py
        ├── filter_rcorrector_se.py
        ├── FormatRockhopTrans.py
        ├── __init__.py
        ├── __init__.pyc
        ├── PlotDESEQ2.Rmd
        ├── PlotEDGER.Rmd
        ├── __pycache__
        │   ├── __init__.cpython-36.pyc
        │   └── utility.cpython-36.pyc
        ├── utility.py
        └── utility.pyc
```


## Test Installation

###### Create Project
python bulkRNASeqWF.py prepareProject --projectName `halomicronema` --sampleMetadataFile `halomicronema_target.txt` --sampleHeader `samples` --groupHeader `group` --refCondName `Control` --local-scheduler


###### Transcript based differential expression analysis
python bulkRNASeqWF.py transcriptomeBasedDEA --projectName halomicronema/ --domain prokaryote --predTranscript yes --readType paired --quantMethod kallisto --genomeName hexcentricum --transcriptName hexcentricum --deaMethod deseq2 --factorInt group --refCond Control --resultTag hexcentricum_treated_vs_control --targetFile halomicronema_target.txt --local-scheduler
