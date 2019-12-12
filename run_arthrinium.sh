##Prepare Project
#python bulkRNASeqWF.py prepareProject --projectName arthrinium --sampleMetadataFile arthrinium_target.txt --sampleHeader samples --groupHeader group --refCondName control --local-scheduler

#Genome Guided Assembly and DE Analysis
#python bulkRNASeqWF.py genomeGuidedDEA --aligner hisat2 --projectName arthrinium --deaMethod deseq2 --resultTag treated_vs_control --readType paired --reportName GGTrinity-DESeq2 --factorInt group --refCond control --domain eukaryote --targetFile arthrinium_target.txt --genomeName arthrinium --local-scheduler 2>&1 | tee arthrinium_genome_guided_deseq2.log

#Transcript based differential expression analysis
#python bulkRNASeqWF.py transcriptomeBasedDEA --projectName arthrinium --domain eukaryote --predTranscript no --readType paired --quantMethod salmon --genomeName arthrinium --transcriptName arthrinium --deaMethod deseq2 --factorInt group --refCond control --resultTag arthrinium_treated_vs_control --targetFile arthrinium_target.txt --local-scheduler 2>&1 | tee arthrinium_transcriptome_salmon_based_DEA_DESeq2.log

#Genome based differential expression analysis
python bulkRNASeqWF.py genomeBasedDEA --projectName arthrinium --domain eukaryote --aligner hisat2 --readType paired --genomeName arthrinium --deaMethod deseq2 --resultTag arthrinium_treated_vs_control --factorInt group --refCond control --targetFile arthrinium_target.txt --local-scheduler 2>&1 | tee arthrinium_genome_based_DEA_DESeq2.log
