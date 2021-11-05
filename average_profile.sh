#!/bin/bash

module load deeptools

computeMatrix reference-point --referencePoint center -R bed_file_1.bed bed_file_2.bed -b 3000 -a 3000 --sortRegions no -S WT.bw KO.bw --skipZeros -o matrix.gz --outFileSortedRegions matrix.bed --numberOfProcessors max

plotProfile -m matrix.gz -out profile.pdf --perGroup --regionsLabel region_1 region_2 --samplesLabel WT KO --yMin 0 --color black red --startLabel center --plotTitle "TITLE"
