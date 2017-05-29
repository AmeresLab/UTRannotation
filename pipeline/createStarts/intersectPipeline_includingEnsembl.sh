cd /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/
mkdir /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/

### divide refSeq annotation into plus and minus strand 

###awk '/+/ { print }' /groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_plus.bed

####awk '/-/ { print }' /groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed |  sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_minus.bed

##### divide ensembl annotation into plus and minus strand

awk '/+/ { print }' /groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/proteinCoding_annotatedUTRs.bed | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/proteinCoding_annotatedUTRs_plus.bed

awk '/-/ { print }' /groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/proteinCoding_annotatedUTRs.bed | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/proteinCoding_annotatedUTRs_minus.bed

cd /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/

## cat the refSeq and encembl and refSeq annotations

cat  proteinCoding_annotatedUTRs_minus.bed | sort -k1,1 -k2,2n > refSeq_mrna_utrsPresent_minus1.bed
cat  proteinCoding_annotatedUTRs_plus.bed | sort -k1,1 -k2,2n > refSeq_mrna_utrsPresent_plus1.bed

ml bedtools

### get the complement of the UTRs (all non-UTR regions)

bedtools complement -i refSeq_mrna_utrsPresent_minus1.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_minus_complement.bed
bedtools complement -i refSeq_mrna_utrsPresent_plus1.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_plus_complement.bed


### complement of complement (all UTR regions - Flattened)

bedtools complement -i refSeq_mrna_utrsPresent_minus_complement.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_minus_complement_complement.bed
bedtools complement -i refSeq_mrna_utrsPresent_plus_complement.bed -g /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/chromosomeSizes_sorted.txt > refSeq_mrna_utrsPresent_plus_complement_complement.bed

### create an identifier to the flattened UTR annotation 

numeberLines_plus=`cat -n refSeq_mrna_utrsPresent_plus_complement_complement.bed | tail -n 1 | cut -f1`
numeberLines_minus=`cat -n refSeq_mrna_utrsPresent_minus_complement_complement.bed | tail -n 1 | cut -f1`

seq 1 $numeberLines_minus > refSeq_mrna_utrsPresent_minus_complement_complement_seq.bed
seq 1 $numeberLines_plus > refSeq_mrna_utrsPresent_plus_complement_complement_seq.bed


## paste the  created identifier to the flattened annotation

paste --delimiters='\t' refSeq_mrna_utrsPresent_plus_complement_complement.bed refSeq_mrna_utrsPresent_plus_complement_complement_seq.bed >refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed
paste --delimiters='\t' refSeq_mrna_utrsPresent_minus_complement_complement.bed refSeq_mrna_utrsPresent_minus_complement_complement_seq.bed >refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed


#### intersect with refSeq mRNAs to get the relationship between the flattened annotation and underlying transcripts.


bedtools intersect -a refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed -b refSeq_mrna_utrsPresent_minus1.bed -wo >refSeq_mrna_utrsPresent_minus_intersect.bed

bedtools intersect -a refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed -b refSeq_mrna_utrsPresent_plus1.bed -wo >refSeq_mrna_utrsPresent_plus_intersect.bed


### (the nonoverlapping annotation is generated in creatingWidths_1)


### intersect nonOverlapping peaks and flattened annotations


awk '/+/ { print }' /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/nonOverlappingCountingWindows_includingEnsembl.bed >nonOverlappingCountingWindows_plus_includingEnsembl.bed
awk '/-/ { print }' /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/nonOverlappingCountingWindows_includingEnsembl.bed >nonOverlappingCountingWindows_minus_includingEnsembl.bed




bedtools intersect -a nonOverlappingCountingWindows_plus_includingEnsembl.bed -b refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed -wo >refSeq_mrna_utrsPresent_plus_complement_complement_countingWindows.bed
bedtools intersect -a nonOverlappingCountingWindows_minus_includingEnsembl.bed -b refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed -wo >refSeq_mrna_utrsPresent_minus_complement_complement_countingWindows.bed

### the counting windows that do not intersect would be intergenic counting windows : amongst this theere can be intergenic counting windows + counting windows that only overlap with genes that do not have a UTR. 

bedtools intersect -a nonOverlappingCountingWindows_plus_includingEnsembl.bed -b refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed > duplicated_plus.bed
bedtools intersect -a nonOverlappingCountingWindows_minus_includingEnsembl.bed -b refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed > duplicated_minus.bed

bedtools intersect -v -a nonOverlappingCountingWindows_plus_includingEnsembl.bed -b refSeq_mrna_utrsPresent_plus_complement_complement_withNumber.bed -wo >refSeq_mrna_utrsPresent_plus_complement_complement_IntergeniccountingWindows.bed
bedtools intersect -v -a nonOverlappingCountingWindows_minus_includingEnsembl.bed -b refSeq_mrna_utrsPresent_minus_complement_complement_withNumber.bed -wo >refSeq_mrna_utrsPresent_minus_complement_complement_IntergeniccountingWindows.bed


### get those overlapping with refSeq genes 

bedtools intersect -s -a refSeq_mrna_utrsPresent_plus_complement_complement_IntergeniccountingWindows.bed -b //groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/annotation/mm10/refSeq_mm10_GRCm38_06-march-2017/refSeqGenes.bed > overlappingWithNonUTRContaining_plus.bed
bedtools intersect -s -a refSeq_mrna_utrsPresent_minus_complement_complement_IntergeniccountingWindows.bed -b //groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/annotation/mm10/refSeq_mm10_GRCm38_06-march-2017/refSeqGenes.bed > overlappingWithNonUTRContaining_minus.bed

#### overlap overlapping counting windows with the very end of refSeq annotations (1n) and +/- 5 nts of counting windows. This is done to get the exact annotation that the counting window 'belongs to'

bedtools intersect -a /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/overlappingCountingWindows_plus_includingEnsembl.bed -b /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/refSeqAnnotation_lastnt_plus_includingEnsembl.bed -wo > overlappingCountingWindows_refSeq_plus.bed

bedtools intersect -a /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/overlappingCountingWindows_minus_includingEnsembl.bed -b /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/refSeqAnnotation_lastnt_minus_includingEnsembl.bed -wo > overlappingCountingWindows_refSeq_minus.bed



