

#### steps :

# divide refSeq UTR annotations into plus and minus
# get complement -- everything other than annotated UTRs
# get another complement - this would be the flattened annotation -- all the UTRs get flattened.
# make an identifier for each flattened annotation
# overlap with flattened annotation and refSeq UTR annotations to get the relationship between the flattened annotations and the transcripts underlying them. 

### there are 2 classes of counting windows identied in the script createWidths_refSeq;
	# overlapping counting winows - <=5nts away from annotated refSeq annotations
	# non overlaping counting winows  >5nts away from annotated refSeq annotations. 



cd /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/
### divide refSeq annotation into plus and minus strand 

awk '/+/ { print }' /groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/annotation/mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed | sort -k1,1 -k2,2n >refSeq_mrna_utrsPresent_plus.bed

awk '/-/ { print }' /groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/annotation/mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed |  sort -k1,1 -k2,2n >refSeq_mrna_utrsPresent_minus.bed

### get the complement of the UTRs (all non-UTR regions)

bedtools complement -i refSeq_mrna_utrsPresent_minus.bed -g chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_minus_complement.bed
bedtools complement -i refSeq_mrna_utrsPresent_plus.bed -g chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_plus_complement.bed


### complement of complement (all UTR regions - Flattened)

bedtools complement -i refSeq_mrna_utrsPresent_minus_complement.bed -g chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_minus_complement_complement.bed
bedtools complement -i refSeq_mrna_utrsPresent_plus_complement.bed -g chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_plus_complement_complement.bed

### create an identifier to the flattened UTR annotation 

numeberLines_plus=`cat -n refSeq_mrna_utrsPresent_plus_complement_complement.bed | tail -n 1 | cut -f1`
numeberLines_minus=`cat -n refSeq_mrna_utrsPresent_minus_complement_complement.bed | tail -n 1 | cut -f1`

seq 1 $numeberLines_minus > refSeq_mrna_utrsPresent_minus_complement_complement_seq.bed
seq 1 $numeberLines_plus > refSeq_mrna_utrsPresent_plus_complement_complement_seq.bed


## paste the  created identifier to the flattened annotation

paste --delimiters='\t' refSeq_mrna_utrsPresent_plus_complement_complement.bed refSeq_mrna_utrsPresent_plus_complement_complement_seq.bed >refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed
paste --delimiters='\t' refSeq_mrna_utrsPresent_minus_complement_complement.bed refSeq_mrna_utrsPresent_minus_complement_complement_seq.bed >refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed


#### intersect with refSeq mRNAs to get the relationship between the flattened annotation and underlying transcripts.


bedtools intersect -a refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed -b refSeq_mrna_utrsPresent_minus.bed -wo >refSeq_mrna_utrsPresent_minus_intersect.bed

bedtools intersect  -a refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed -b refSeq_mrna_utrsPresent_plus.bed -wo >refSeq_mrna_utrsPresent_plus_intersect.bed


### (the nonoverlapping annotation is generated in creatingWidths_1)


### intersect nonOverlapping peaks and flattened annotations


awk '/+/ { print }' nonOverlappingCountingWindows.bed >nonOverlappingCountingWindows_plus.bed
awk '/-/ { print }' nonOverlappingCountingWindows.bed >nonOverlappingCountingWindows_minus.bed




bedtools intersect  -a nonOverlappingCountingWindows_plus.bed -b refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed -wo >refSeq_mrna_utrsPresent_plus_complement_complement_countingWindows.bed
bedtools intersect  -a nonOverlappingCountingWindows_minus.bed -b refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed -wo >refSeq_mrna_utrsPresent_minus_complement_complement_countingWindows.bed


bedtools intersect -v  -a nonOverlappingCountingWindows_plus.bed -b refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed -wo >nonOverlappingWithRefSeq_plus.bed
bedtools intersect -v  -a nonOverlappingCountingWindows_minus.bed -b refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed -wo >nonOverlappingWithRefSeq_minus.bed



#### overlap overlapping counting windows with the very end of refSeq annotations (1n) and +/- 5 nts of counting windows. This is done to get the exact annotation that the counting window 'belongs to'

bedtools intersect -s -a overlappingCountingWindows_plus.bed -b refSeqAnnotation_lastnt_plus.bed -wo > overlappingCountingWindows_refSeq_plus.bed

bedtools intersect -s -a overlappingCountingWindows_minus.bed -b refSeqAnnotation_lastnt_minus.bed -wo > overlappingCountingWindows_refSeq_minus.bed


