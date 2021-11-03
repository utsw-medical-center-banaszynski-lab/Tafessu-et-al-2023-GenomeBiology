#!/bin/bash

module load deeptools



computeMatrix reference-point --referencePoint center -R bed_file_1.bed bed_file_2.bed -b 3000 -a 3000 --sortRegions no -S ATAC_H3.3WT.bw ATAC_H3.3KO.bw --skipZeros -o matrix_ATAC_1.gz --outFileSortedRegions matrix_ATAC_1.bed --numberOfProcessors max

plotProfile -m matrix_ATAC_1.gz -out ATAC_profile.pdf --perGroup --regionsLabel region_1 region_2 --samplesLabel WT H3.3KO --yMin 0 --color black red --startLabel center --plotTitle "ATAC"
