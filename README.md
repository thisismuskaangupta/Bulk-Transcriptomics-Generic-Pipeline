# Bulk-Transcriptomics-Generic-Pipeline
this contains the code written for the MSc Bioinformatics course - RNA-seq and Next Generation Transcriptomics. grade: A5

in this assignment, generic fastq files were provided, and the following was done - 
1) read trimming using sickle and scythe
2) 'New Tuxedo' pipeline (hisat2 and stringtie) for read alignment and assembly
3) count matrices were generated from assembly using an external 'prepDE' python script
4) (based on information in instructions) an experiment design table was generated
5) filtering was done and DE analysis was run
6) dispersion plots were created
7) rlog-based PCA plots were created
8) “SD versus mean” plots were created
9) MA (MvA) plots were created based on null hypotheses provided
10) limma-based batch removal function was applied to correct rlog-transformed gene-counts
11) PCA plots were created with corrected rlog values.

Note - Unfortunately these fastq files are no longer available to test, but this can used as a skeleton to apply to a bulk transcriptomics experiment.
