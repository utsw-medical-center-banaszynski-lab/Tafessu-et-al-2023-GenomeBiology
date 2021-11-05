#!/bin/bash

module load macs

bams='bam1 bam2 bam3 bam4'

for bam in $bams

do

macs2 callpeak -t ${bam}.bam -g mm -f BAM -n ./${bam} -q 0.01

done
