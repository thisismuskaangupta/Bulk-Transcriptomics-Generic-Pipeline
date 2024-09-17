#! /bin/bash
#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=48:00:00
#PBS -N singleHPCscript
#PBS -d /export/home/biostuds/2873826g/rnaseq/part1assignment
#PBS -m abe
#PBS -M 2873826g@student.gla.ac.uk
#PBS -q bioinf-stud
#
# FILES AND DIRECTORY TO REFER TO WHILE RUNNING THE SCRIPT
illumina_adapter='/export/projects/polyomics/biostuds/data/illumina_adapter.fa'
hisat2index='/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2'
gtf_annotations='/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf'
data='.'
#
# MAKING SUBDIRECTORIES TO STORE THE RESULTS
hisat_results_dir='./hisat2'
stringtie_results_dir='./stringtie'
mkdir -p ${hisat_results_dir}
mkdir -p ${stringtie_results_dir}
#
gtflist='list.gtf.txt' # HARDCODING THE FILENAME OF THE FINAL GTF LIST rm -f 
${gtflist} # REMOVING IT IF IT EXISTS ALREADY
#
#EXECUTING THE PIPELINE
for sample in s1.c2 s2.c2 s3.c2 s4.c2 s5.c2 s6.c2 s7.c2 s8.c2 s9.c2 s10.c2 s11.c2 s12.c2
do
    raw_fastq="${data}/${sample}.fq"
    trim1_fastq="${data}/${sample}.t1.fq"
    trim2_fastq="${data}/${sample}.t2.fq"
    unsorted_bam="${hisat_results_dir}/${sample}.bam"
    sam_file="${hisat_results_dir}/${sample}.sam"
    sorted_bam="${hisat_results_dir}/${sample}.sort.bam"
    scythe -q sanger -a ${illumina_adapter} -o ${trim1_fastq} ${raw_fastq}
    sickle se -f ${trim1_fastq} -t sanger -o ${trim2_fastq} -q 10 -l 48
    hisat2 -p 4 --phred33 --rna-strandness R -x ${hisat2index} -U ${trim2_fastq} -S ${sam_file}
    samtools view -b -o ${unsorted_bam} ${sam_file}
    samtools sort -o ${sorted_bam} ${unsorted_bam}
    rm ${sam_file} ${unsorted_bam} #REMOVING ALL THE EXTRA FILES NOT NEEDED ANYMORE, EXCEPT RAW INPUT FILES AND SORTED BAM FILE.
    rm ${trim1_fastq} ${trim2_fastq}
    str_results_dir="${stringtie_results_dir}/${sample}"
    mkdir -p ${str_results_dir}
    sample_transcripts_gtf="${str_results_dir}/${sample}_transcripts.gtf"
    stringtie -p 4 -t -e -B --rf -G ${gtf_annotations} -o ${sample_transcripts_gtf} ${sorted_bam}
    gtfline="${sample} ${sample_transcripts_gtf}"
    echo ${gtfline} >> ${gtflist}
done
#
prepDE.py -i ${gtflist}
