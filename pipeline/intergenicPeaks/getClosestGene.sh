#!/bin/bash

module load bedtools

##### sort the bed files 

## this is the bed file of all ensembl and refSeq transcripts. 

sort -k1,1 -k2,2n $BOUT/intergenicPeaks/allExons_refSeq_ensembl.bed > $BOUT/intergenicPeaks/allExons_refSeq_ensembl_sorted.bed

### this is the bed file of the most distal 3' position per gene. 

sort -k1,1 -k2,2n $BOUT/intergenicPeaks/toExtend_longestEnsembl_refSeq_n100.bed > $BOUT/intergenicPeaks/toExtend_longestEnsembl_refSeq_n100_sorted.bed


##### we want to calculate the distance between the most distal 3' end per gene and the next annotation (refSeq + ensembl), to prevent considering RNAseq singal coming from another annotation. 

bedtools closest -d -s -io -iu -D a -a $BOUT/intergenicPeaks/toExtend_longestEnsembl_refSeq_n100_sorted.bed -b $BOUT/intergenicPeaks/allExons_refSeq_ensembl_sorted.bed > $BOUT/intergenicPeaks/toExtend_longestEnsembl_refSeq_n100_sorted_distances.bed



