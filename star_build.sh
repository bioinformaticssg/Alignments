#!/bin/bash
#$ -N star_build
#$ -o /data/users/$USER/BioinformaticsSG/Getting-Ready-for-Analysis/MY_PROJECT/data/griffith_data/refs/star_build.out # the name of the output file. not to be confused with the results.
#$ -e /data/users/$USER/BioinformaticsSG/Getting-Ready-for-Analysis/MY_PROJECT/data/griffith_data/refs/star_build.err # name of the error file
#$ -q pub8i
#$ -pe openmp 8
#$ -m beas
#$ -ckpt blcr


# Make the script stop on any error.
set -uxeo pipefail
echo $HOSTNAME

#goal: build the STAR aligner specific indeces required prior to performing alignments


module load blcr
module load STAR/2.5.2a
module load enthought_python/7.3.2
module load samtools/1.3

DATA_DIR=/data/users/$USER/BioinformaticsSG/Getting-Ready-for-Analysis/MY_PROJECT/data/griffith_data

#reference genome file
REF_FASTA=${DATA_DIR}/refs/chr22.ERCC92.fa
#annotation file
REF_ANNOTATION=${DATA_DIR}/refs/chr22.ERCC92.gtf
#where to place the resulting indeces
OUTPUT_DIR=${DATA_DIR}/refs/
O=99 #this overhang ideally should be ReadLength-1.
P=8 #threads



STAR --runMode genomeGenerate --genomeDir $OUTPUT_DIR --genomeFastaFiles $REF_FASTA --runThreadN $P --sjdbGTFfile $REF_ANNOTATION --sjdbOverhang $O
