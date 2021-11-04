
module load deeptools


computeMatrix reference-point --referencePoint center -R bed_file.bed -b 3000 -a 3000 \
-S WT.bw \
KO.bw \
--missingDataAsZero -o matrix.gz \
--numberOfProcessors max
plotHeatmap -m matrix.gz -out heatmap.pdf --sortRegions no --whatToShow 'heatmap and colorbar' --colorMap OrRd OrRd --samplesLabel WT KO --plotTitle "Title"
