#!/bin/bash

#--PBS JOB CONFIGURATION--
# Resource allocation for execution on a PBS-managed HPC cluster
# The original coursework analysis was executed using a university HPC environment
# Institution-specific paths and account details have been replaced with generic placeholders for this repository

#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=48:00:00
#PBS -N rna_seq_pipeline
#PBS -d /path/to/project
#PBS -m abe
# Email notifications were configured for the original HPC environment and have been omitted from this public version
#PBS -q bioinf-stud

# --RESOURCE FILES--
# Illumina adapter sequences used by Scythe for adapter trimming
adapter='/path/to/illumina_adapter.fa'

# HISAT2 index for chromosome 2 of the Mus musculus mm10 reference genome
hs2index='/path/to/mm10/Hisat2Index/chr2'

# GTF annotation file containing gene/transcript annotations for the reference genome
gtf='/path/to/mm10/annotations/chr2.gtf'

# Directory containing raw FASTQ files
fqpath='/path/to/assessment_data'


# --OUTPUT DIRECTORIES--
# Directory for input/read files and intermediate FASTQ files
data='../data'

# Directory for HISAT2 alignment and BAM files
hisat_dir='./hisat2'

# Directory for StringTie output
stringtie_dir='./stringtie'

mkdir -p ${data}
mkdir -p ${hisat_dir}
mkdir -p ${stringtie_dir}


# --PREPARE GTF FILE LIST--
# Filename for the list of sample-specific GTF files used by prepDE.py
gtflist='list.gtf.txt'

rm -f ${gtflist}


# --MAIN ANALYSIS LOOP--

# Process samples s1-s9.
for sample in s1 s2 s3 s4 s5 s6 s7 s8 s9

do
        # Define input FASTQ file and symbolic link
        # Path to the raw FASTQ file for the current sample
        raw_file="${fqpath}/${sample}.c2.fq"

        # Name of the symbolic link within the working directory
        link_name="${data}/${sample}.fq"

        # Remove existing symbolic link if present
        rm -f ${link_name}

        # Create symbolic link to the raw FASTQ file
        ln -s ${raw_file} ${link_name}

        # Define intermediate files
        # Input and trimmed FASTQ files
        fastq="${data}/${sample}.fq"
        trim1="${data}/${sample}.t1.fq"
        trim2="${data}/${sample}.t2.fq"

        # Alignment output files
        bam="${hisat_dir}/${sample}.bam"
        sam="${hisat_dir}/${sample}.sam"
        sorted_bam="${hisat_dir}/${sample}.sort.bam"

        # Read preprocessing
        # Remove Illumina adapter sequences using Scythe
        scythe -q sanger -a ${adapter} -o ${trim1} ${fastq}

        # Perform quality and length filtering using Sickle
        # Reads are retained with a minimum quality score of 10 and a minimum length of 50 bases
        sickle se -f ${trim1} -t sanger -o ${trim2} -q 10 -l 50


        # RNA-seq read alignment
        # Align processed RNA-seq reads to the Mus musculus mm10 chromosome 2 reference using HISAT2
        hisat2 -p 4 --phred33 --rna-strandness R --dta \
            -x ${hs2index} -U ${trim2} -S ${sam}

        # SAM/BAM processing
        # Convert SAM alignment file to BAM format
        samtools view -b -o ${bam} ${sam}

        # Sort BAM file by genomic coordinates
        samtools sort -o ${sorted_bam} ${bam}

        # Remove intermediate SAM and unsorted BAM files
        rm ${sam} ${bam}

        # Remove intermediate trimmed FASTQ files
        rm ${trim1} ${trim2}

        # Transcript assembly and quantification
        # Create a separate StringTie output directory for the current sample
        str_smp_dir="${stringtie_dir}/${sample}"
        mkdir -p ${str_smp_dir}

        # Define output transcript GTF file
        sample_tr_gtf="${str_smp_dir}/${sample}_transcripts.gtf"

        # Assemble and quantify transcripts using StringTie
        stringtie -p 4 --rf -e -B -G ${gtf} \
            -o ${sample_tr_gtf} ${sorted_bam}

        # Prepare output for downstream analysis
        # Create a line containing the sample name and corresponding transcript GTF file
        gtfline="${sample} ${sample_tr_gtf}"

        # Append information to the list used by prepDE.py.
        echo ${gtfline} >> ${gtflist}

done

# --GENERATE COUNT MATRICES--
# Generate gene- and transcript-level expression matrices from the StringTie GTF files using prepDE.py
prepDE.py -i ${gtflist}
