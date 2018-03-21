#!/bin/bash
#$ -N PE_staralign_local
#$ -o /data/users/$USER/BioinformaticsSG/Getting-Ready-for-Analysis/MY_PROJECT/alignments/staralign.out # the name of the output file. not to be confused with the results.
#$ -e /data/users/$USER/BioinformaticsSG/Getting-Ready-for-Analysis/MY_PROJECT/alignments/staralign.err # name of the error file
#$ -q pub8i
#$ -pe openmp 8
#$ -m beas
#$ -ckpt blcr

# set tells the program to stop if it runs into any issues, the flag descriptions are as follows:
# -e Exit immediately when a command fails.
# -u Treats unset/unbound variables as an error and exits immediately.
# -x Prints each command before executing it -- this is helpful for debugging.
# -o Sets the exit code (0 is successful) to the that of the rightmost command...continued on next line
# -o If the first command in a pipeline failed it will be carried through to the end of the pipeline and still exit instead of continuing through the script

set -euxo pipefail

module load blcr
module load STAR/2.5.2a
module load enthought_python/7.3.2
module load samtools/1.3

echo "You're aligning on $HOSTNAME"

#number of threads. MUST match the qsub header openmp value you chose
P=8
#where the fastq sample (reads) sequences are
#data is paired end
DATA_DIR=/data/users/$USER/BioinformaticsSG/Getting-Ready-for-Analysis/MY_PROJECT/data/griffith_data

#annotation file
REF_ANNOTATION=${DATA_DIR}/refs/chr22.ERCC92.gtf

#path to the index files for star. this replaces the reference genome because it is an indexed version of the fasta file
INDEX=${DATA_DIR}/refs/

#desired directory for output alignments
OUT_DIR=/data/users/$USER/BioinformaticsSG/Getting-Ready-for-Analysis/MY_PROJECT/alignments


#####################

# Iterate over each sample
for SAMPLE in HBR UHR;
do
    # Iterate over each of the replicates.
    for REPLICATE in 1 2 3;
    do
        # Build the name of the files.
        R1=${DATA_DIR}/reads/${SAMPLE}_${REPLICATE}_R1.fq.gz
        R2=${DATA_DIR}/reads/${SAMPLE}_${REPLICATE}_R2.fq.gz

        # Run the aligner.
		echo "Aligning ${R1} and ${R2}"


		STAR --chimSegmentMin 15 --chimJunctionOverhangMin 15 --outFilterMismatchNmax 3 --alignEndsType Local  --runThreadN $P \
		--outFilterMultimapNmax 1 --outBAMcompression 10 --outBAMsortingThreadN $P \
		--quantMode GeneCounts  TranscriptomeSAM  --quantTranscriptomeBAMcompression 10 \
		--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --alignSJDBoverhangMin 6 --alignIntronMax 300000 --genomeDir $INDEX \
		--sjdbGTFfile $REF_ANNOTATION --outFileNamePrefix ${OUT_DIR}/${SAMPLE}_${REPLICATE}.star \
		--readFilesCommand zcat --readFilesIn ${R1} ${R2}

		echo "Alignment of ${R1} and ${R2} complete"
    done
done

exit
########################################

#STAR --chimSegmentMin 15 --chimJunctionOverhangMin 15 \ #allows for circRNA detection
#--outFilterMismatchNmax 3 \ #maximum mismatches allowed
#--alignEndsType Local  --runThreadN $P \ #type of alignment and number of threads.Local is standard local alignment with soft-clipping allowed
#--quantMode GeneCounts  TranscriptomeSAM  --quantTranscriptomeBAMcompression 10 \ #performs gene counts like HTSEQ union mode. Just exons counts. max compression
#--outFilterMultimapNmax 1 \ # default is 10.setting this to 1 will limit the bam file to contain only TRULY uniquely mapped reads. Stats will still contain multimap info, but the definition of a multipmapper more correct by setting this to 1 instead of 10 $
#--outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate \
#--outBAMcompression 10 --outBAMsortingThreadN $P \ #these features increase compression and compression speed by multithreading more than the default 6 threads
#--alignSJDBoverhangMin 6 --alignIntronMax 300000 \ #default intron is 0. if this is 0 it turns off splicing detection.  default overhang is 5

#NOTES:
#--outSAMstrandField intronMotif needed if non stranded data
#.run this if error produced. input error amount --limitBAMsortRAM 30943606211
#note if you run out of ram, tactic is to take out BAM sorted option to yield unsorted then do it separately with samtools

