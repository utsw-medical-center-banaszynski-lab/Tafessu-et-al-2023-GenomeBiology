#!/bin/bash

./feature-counts.sh -f WT_rep1.bam -s paired-end -a mm10.ncbiRefSeq.gtf -o ./WT_rep1
./feature-counts.sh -f WT_rep2.bam -s paired-end -a mm10.ncbiRefSeq.gtf -o ./WT_rep2

./feature-counts.sh -f KO_rep1.bam -s paired-end -a mm10.ncbiRefSeq.gtf -o ./KO_rep1
./feature-counts.sh -f KO_rep2.bam -s paired-end -a mm10.ncbiRefSeq.gtf -o ./KO_rep2
