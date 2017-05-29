#!/bin/bash

### procedure to identify duplicates :

	# separate plus and minus refGenes.
	# sort the refGenes
	# completment the refGenes to get all non-gene regions
	# complement again to flatten the refGenes 
	# add identifier to each of the flattened annotations.
	# overlap multi transcript files identified in creatingWidths_1.R

### split refSeq genes based on strand

awk '/+/ { print }' /groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/annotation/mm10/refSeq_mm10_GRCm38_06-march-2017/refSeqGenes.bed | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/refSeqGenes_plus.bed
awk '/-/ { print }' /groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/annotation/mm10/refSeq_mm10_GRCm38_06-march-2017/refSeqGenes.bed  | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/refSeqGenes_minus.bed


cd /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/


cat  refSeqGenes_plus.bed | sort -k1,1 -k2,2n > refSeqGenes_plus1.bed
cat  refSeqGenes_minus.bed | sort -k1,1 -k2,2n > refSeqGenes_minus1.bed



### get the complement of the UTRs (all non-gene regions)

bedtools complement -i refSeqGenes_plus1.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeqGenes_plus_complement.bed
bedtools complement -i refSeqGenes_minus1.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeqGenes_minus_complement.bed


### complement of complement (all gene regions - Flattened)

bedtools complement -i refSeqGenes_plus_complement.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeqGenes_plus_complement_complement.bed
bedtools complement -i refSeqGenes_minus_complement.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeqGenes_minus_complement_complement.bed


### create an identifier to the flattened UTR annotation 

numeberLines_plus=`cat -n refSeqGenes_plus_complement_complement.bed | tail -n 1 | cut -f1`
numeberLines_minus=`cat -n refSeqGenes_minus_complement_complement.bed | tail -n 1 | cut -f1`

seq 1 $numeberLines_plus > refSeqGenes_plus_complement_complement_seq.bed
seq 1 $numeberLines_minus > refSeqGenes_minus_complement_complement_seq.bed


## paste the  created identifier to the flattened annotation

paste --delimiters='\t' refSeqGenes_plus_complement_complement.bed refSeqGenes_plus_complement_complement_seq.bed >refSeqGenes_plus_complement_complement_withNumber.bed
paste --delimiters='\t' refSeqGenes_minus_complement_complement.bed refSeqGenes_minus_complement_complement_seq.bed >refSeqGenes_minus_complement_complement_withNumber.bed

#### overlap with refSeq genes that have either duplicates or UTR introns


#### intersect with refSeq mRNAs to get the relationship between the flattened annotation and underlying transcripts.


bedtools intersect -b refSeqGenes_plus_complement_complement_withNumber.bed -a refSeq_multipleTranscripts_plus.bed -wo >intersectMultipleTranscripts_plus.bed

bedtools intersect -b refSeqGenes_minus_complement_complement_withNumber.bed -a refSeq_multipleTranscripts_minus.bed -wo >intersectMultipleTranscripts_minus.bed

