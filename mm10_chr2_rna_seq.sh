#!/bin/bash


#PBS JOB CONFIGURATION
#PBS -l nodes=1:ppn=4:centos7,cput=24:00:00,walltime=48:00:00
#PBS -N links
#PBS -d /export/biostuds/2658500b/rna_assignment/nondiscovery
#PBS -m abe
#PBS -M 2658500b@student.gla.ac.uk
#PBS -q bioinf-stud
#
#RESOURCE FILES
#illumina adapter sequences used by Scythe for adapter trimming
adapter='/export/projects/polyomics/biostuds/data/illumina_adapter.fa'

#HISAT2 index for chromosome 2 of Mus musculus mm10 reference genome
hs2index='/export/projects/polyomics/Genome/Mus_musculus/mm10/Hisat2Index/chr2'

#GTF annotation file containing gene/transcript annotations for reference genome
gtf='/export/projects/polyomics/Genome/Mus_musculus/mm10/annotations/chr2.gtf'

#directory containing raw FASTQ files
fqpath='/export/other/polyomics/miRNAResources/assessment_data'

#MAKE OUTPUT DIRECTORIES
#directory for input/read files and intermediate FASTQ files
data='../data'

#direcoty for HISAT2 alignment and BAM files
hisat_dir='./hisat2'

#directory for stringtie output
stringtie_dir='./stringtie'

mkdir -p ${data}
mkdir -p ${hisat_dir}
mkdir -p ${stringtie_dir}

#filename for final GTF list
gtflist='list.gtf.txt'
rm -f ${gtflist}

#MAIN ANALYSIS LOOP
#process samples s1-s9
for sample in s1 s2 s3 s4 s5 s6 s7 s8 s9

do
	raw_file="${fqpath}/${sample}.c2.fq"    #path to raw fastq file
        link_name="${data}/${sample}.fq"        #name it will appear as in directory
        #remove link if it exists already
        rm -f ${link_name}
        #create link
        ln -s ${raw_file} ${link_name}

        #create test file
        #test_file="${sample}.25K.fq"
        #head -n 100000 ${link_name} > ${test_file}

	#define intermediate files
	#input and trimmed FASTQ files
	fastq="${data}/${sample}.fq"
	trim1="${data}/${sample}.t1.fq"
	trim2="${data}/${sample}.t2.fq"

	#alignment output files
	bam="${hisat_dir}/${sample}.bam"
	sam="${hisat_dir}/${sample}.sam"
	sorted_bam="${hisat_dir}/${sample}.sort.bam"

	#trim reads
	#Scythe removes illumina adapter sequences
	scythe -q sanger -a ${adapter} -o ${trim1} ${fastq}

	#perform quality and length filtering using Sickle
	sickle se -f ${trim1} -t sanger -o ${trim2} -q 10 -l 50
	
	#align the processed RNA-seq reads to mm10 chromosome 2 reference with hisat2
	hisat2 -p 4 --phred33 --rna-strandness R --dta -x ${hs2index} -U ${trim2} -S ${sam}
	
	#convert SAM alignment file to BAM format
	samtools view -b -o ${bam} ${sam}
	#sort BAM file by genomic coordinates
	samtools sort -o ${sorted_bam} ${bam}
	#remove intermediate SAM and unsorted BAM
	rm ${sam} ${bam}
	#remove intermediate trimmed FASTQ files
	rm ${trim1} ${trim2}

	#transcript assembly
	#create seperate StringTie output directory for the sample
	str_smp_dir="${stringtie_dir}/${sample}"
	mkdir -p ${str_smp_dir}
	#define output transcript GTF file
	sample_tr_gtf="${str_smp_dir}/${sample}_transcripts.gtf"
	#run StringTie to assemble and quantify transcripts
	stringtie -p 4 --rf -e -B -G ${gtf} -o ${sample_tr_gtf} ${sorted_bam}
	
	#record output for downstream analysis
	#create a line containing sample name and corresponding transcript GTF file
	gtfline="${sample} ${sample_tr_gtf}"
	#append information to list used by prepDE.py
	echo ${gtfline} >> ${gtflist}
done

#generate expression matrices using StringTie GTF files
prepDE.py -i ${gtflist}