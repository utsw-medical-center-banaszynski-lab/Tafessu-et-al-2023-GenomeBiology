#!/bin/tcsh

set INDEX=genome.fa
set NUM_THREADS=30
set MISMATCH_PENALTY=8
set FILE=./fastq

module load BWA/0.7.5
module load samtools
module load picard/1.117
module load bedtools
module load UCSC_userApps
module load macs/1.4.2
module load igvtools/2.3.71

mkdir bam_files

foreach sample(\
	sample1\
	sample2\
)


bwa mem\
	 -t $NUM_THREADS\
	 -B $MISMATCH_PENALTY\
	 -M\
	 $INDEX\
	 $FILE/$sample.fastq.gz\
	 | samtools view -q10 -bS -o - -\
	 | samtools sort - -o ./bam_files/$sample.bam

    samtools sort -m 8000000000 ./bam_files/$sample.bam

    set MARK_DUP = /cm/shared/apps/picard/1.117/MarkDuplicates.jar

    mkdir nodup_files

    java -Xmx128g -jar $MARK_DUP\
    INPUT=./bam_files/$sample.bam\
    OUTPUT=./nodup_files/$sample.nodup.bam\
    METRICS_FILE=./nodup_files/metrics.$sample.txt\
    REMOVE_DUPLICATES=true\
    ASSUME_SORTED=true\
    TMP_DIR=temp_dir.$sample

    samtools index ./nodup_files/$sample.nodup.bam

    set CHROM_SIZE = genome.fa.fai

    mkdir bw_files

    bedtools genomecov -ibam ./nodup_files/$sample.nodup.bam\
			-bga\
			> ./bw_files/$sample.nodup.bedGraph

    sort -k1,1 -k2,2n ./bw_files/$sample.nodup.bedGraph > ./bw_files/$sample.nodup.sorted.bedGraph

    bedGraphToBigWig ./bw_files/$sample.nodup.sorted.bedGraph\
		     $CHROM_SIZE\
		     ./bw_files/$sample.nodup.bw
end
