# bulkRNASeqPIPE
Differential Expression Analysis of bulk-RNASeq data

git clone https://github.com/sarangian/bulkRNASeqPIPE.git
cd bulkRNASeqPIPE
chmod 755 ubuntu_install.sh
./ubuntu_install.sh

After successful installation please source ~/.bashrc

Test Installation
-----------------

#Create Project

python3 bulkRNASeqWF.py prepareProject --projectName halomicronema --sampleMetadataFile halomicronema_target.txt --sampleHeader samples --groupHeader group --refCondName Control --local-scheduler


#Transcript based differential expression analysis

python3 bulkRNASeqWF.py transcriptomeBasedDEA --projectName halomicronema/ --domain prokaryote --predTranscript yes --readType paired --quantMethod kallisto --genomeName hexcentricum --transcriptName hexcentricum --deaMethod deseq2 --factorInt group --refCond Control --resultTag hexcentricum_treated_vs_control --targetFile halomicronema_target.txt --local-scheduler
