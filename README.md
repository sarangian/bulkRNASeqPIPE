# bulkRNASeqPIPE
Differential Expression Analysis of bulk-RNASeq data


# Dependencies

Installation of bulkRNASeqPIPE, requires libboost-all-dev, libbz2-dev, and liblzma-dev installed in your system.


# Download and install

git clone https://github.com/sarangian/bulkRNASeqPIPE.git

cd bulkRNASeqPIPE

chmod 755 miniconda_x64_bulkRNASeq_install.sh

./miniconda_x64_bulkRNASeq_install.sh


# Post installation instruction

After successful installation. close the current terminal. In a new terminal please source ~/.bashrc
and activate bulkRNASeqPIPE environment using command "conda activate"

# Test Installation

# Create Project
python3 bulkRNASeqWF.py prepareProject --projectName halomicronema --sampleMetadataFile halomicronema_target.txt --sampleHeader samples --groupHeader group --refCondName Control --local-scheduler


# Transcript based differential expression analysis
python3 bulkRNASeqWF.py transcriptomeBasedDEA --projectName halomicronema/ --domain prokaryote --predTranscript yes --readType paired --quantMethod kallisto --genomeName hexcentricum --transcriptName hexcentricum --deaMethod deseq2 --factorInt group --refCond Control --resultTag hexcentricum_treated_vs_control --targetFile halomicronema_target.txt --local-scheduler
