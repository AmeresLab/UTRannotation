### intersect nonOverlapping counting windows


############################
## nonOverlapping counting windows 
###########################

### reading in refSeq annotation

#refSeqAnnotation = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed",stringsAsFactors = F)
ensemblAnnotations = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/proteinCoding_annotatedUTRs.bed",stringsAsFactors = F)
### reading in the non-overlapping counting windows that have been intersected with flattened UTR annotations.

countingWindowOverlap_plus = read.table("/Volumes//groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_plus_complement_complement_countingWindows.bed")
countingWindowOverlap_minus = read.table("/Volumes//groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_minus_complement_complement_countingWindows.bed")


#### connecting trasncript to fragment - reading in the intersection between  flattened 'fragments'and UTR anntoation:



transcriptToFragment_minus= read.table("/Volumes//groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_minus_intersect.bed",sep="\t",stringsAsFactors = F)
transcriptToFragment_plus= read.table("/Volumes//groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_plus_intersect.bed",sep="\t",stringsAsFactors = F)
transcriptToFragment_plus$V4 = paste0("F",transcriptToFragment_plus$V4)
transcriptToFragment_minus$V4 = paste0("F",transcriptToFragment_minus$V4)

### just getting the fragment and the 'transcipt names'

transcriptToFragment_minus_connection = cbind.data.frame(transcriptToFragment_minus[,c("V4","V11")])
#transcriptToFragment_minus_connection$V4 = paste0("F",transcriptToFragment_minus_connection$V4)

### we also overlapped the counting windows with fragments. Adding 'F' to the fragment identifiers that overlap with counting windows

countingWindowOverlap_minus$V10 = paste0("F",countingWindowOverlap_minus$V10)

### just add these columns: 

countingWindowOverlap_minus$utrStr_chr = NA
countingWindowOverlap_minus$utrStr_start = NA
countingWindowOverlap_minus$utrStr_end = NA
countingWindowOverlap_minus$transcriptName = NA
### go through all the counting windows : for each counting window, check which fragment it overlaps with and which transcript the fragment is made of.



for(i in 1:nrow(countingWindowOverlap_minus)){
  ##check which fragment it overlaps with and which transcript the fragment is made of.
  transcriptThatOverlaps = transcriptToFragment_minus_connection[which(transcriptToFragment_minus_connection$V4 == countingWindowOverlap_minus$V10[i]),]$V11
  name_transcripts = paste(transcriptThatOverlaps,collapse = "|")
  ## check if this transcript is a part of any other fragment
  doesTranscriptExistInOtherFragment= transcriptToFragment_minus_connection[is.element(transcriptToFragment_minus_connection$V11,transcriptThatOverlaps),]
  
  #doesTranscriptExistInOtherFragment = transcriptToFragment_minus_connection[which(transcriptToFragment_minus_connection$V11 == transcriptThatOverlaps),]
  ## if the transctipt is only part of one fragment, just get the start and end of this fragment.
  ### if the transcrips is present in mpre than one fragment - 
  if(length(unique(doesTranscriptExistInOtherFragment$V4))==1){
    cat("transcript exists only in one fragment")
    countingWindowOverlap_minus[i,]$V3 = countingWindowOverlap_minus[i,]$V9
    countingWindowOverlap_minus[i,]$utrStr_chr = paste(countingWindowOverlap_minus[i,]$V7,collapse = "|")
    countingWindowOverlap_minus[i,]$utrStr_start = paste(countingWindowOverlap_minus[i,]$V8,collapse = "|")
    countingWindowOverlap_minus[i,]$utrStr_end = paste(countingWindowOverlap_minus[i,]$V9,collapse = "|")
    transcriptsInRefSeq = ensemblAnnotations[is.element(ensemblAnnotations$V7, doesTranscriptExistInOtherFragment$V11),]
    transcriptsInRefSeq = transcriptsInRefSeq[which.max(transcriptsInRefSeq$V3),]
    countingWindowOverlap_minus[i,]$transcriptName = transcriptsInRefSeq$V7
  
  }else{
    conuntingWindowsConsidered = transcriptToFragment_minus[is.element(transcriptToFragment_minus$V4,doesTranscriptExistInOtherFragment$V4),]
    #conuntingWindowsConsidered = countingWindowOverlap_plus[is.element(countingWindowOverlap_plus$V10,doesTranscriptExistInOtherFragment$V4),]
    conuntingWindowsConsidered =as.data.frame( conuntingWindowsConsidered %>% group_by(V4) %>% mutate(collapsed = paste(V11, collapse = "|")))
    a = conuntingWindowsConsidered[which.max(conuntingWindowsConsidered$V3),]
    conuntingWindowsConsidered = conuntingWindowsConsidered[which(conuntingWindowsConsidered$V11 == a$V11),]
    
    #conuntingWindowsConsidered$V11 = conuntingWindowsConsidered$collapsed
    conuntingWindowsConsidered$collapsed = NULL
    conuntingWindowsConsidered = conuntingWindowsConsidered[!duplicated(conuntingWindowsConsidered),]
    conuntingWindowsConsidered =  conuntingWindowsConsidered[which(conuntingWindowsConsidered$V8 == countingWindowOverlap_minus[i,]$V4),]
    countingWindowOverlap_minus[i,]$utrStr_chr = paste(conuntingWindowsConsidered$V5,collapse = "|")
    countingWindowOverlap_minus[i,]$utrStr_start = paste(conuntingWindowsConsidered$V6,collapse = "|")
    countingWindowOverlap_minus[i,]$utrStr_end = paste(conuntingWindowsConsidered$V7,collapse = "|")
    countingWindowOverlap_minus[i,]$V3 = max(conuntingWindowsConsidered$V3)
    countingWindowOverlap_minus[i,]$transcriptName = conuntingWindowsConsidered$V11[1]
    cat(i)
  }
  print(i)
  
}



#### plus strand 


transcriptToFragment_plus_connection = cbind.data.frame(transcriptToFragment_plus[,c("V4","V11")])
#transcriptToFragment_plus_connection$V4 = paste0("F",transcriptToFragment_plus_connection$V4)
countingWindowOverlap_plus$V10 = paste0("F",countingWindowOverlap_plus$V10)
countingWindowOverlap_plus$utrStr_chr = NA
countingWindowOverlap_plus$utrStr_start = NA
countingWindowOverlap_plus$utrStr_end = NA
countingWindowOverlap_plus$transcriptName = NA

for(i in 1:nrow(countingWindowOverlap_plus)){
  
  transcriptThatOverlaps = transcriptToFragment_plus_connection[which(transcriptToFragment_plus_connection$V4 == countingWindowOverlap_plus$V10[i]),]$V11
  name_transcripts = paste(transcriptThatOverlaps,collapse = "|")
  doesTranscriptExistInOtherFragment= transcriptToFragment_plus_connection[is.element(transcriptToFragment_plus_connection$V11,transcriptThatOverlaps),]
  
  #doesTranscriptExistInOtherFragment = transcriptToFragment_minus_connection[which(transcriptToFragment_minus_connection$V11 == transcriptThatOverlaps),]
  if(length(unique(doesTranscriptExistInOtherFragment$V4))==1){
    cat("transcript exists only in one fragment")
    countingWindowOverlap_plus[i,]$V2 = countingWindowOverlap_plus[i,]$V8
    countingWindowOverlap_plus[i,]$utrStr_chr = paste(countingWindowOverlap_plus[i,]$V7,collapse = "|")
    countingWindowOverlap_plus[i,]$utrStr_start = paste(countingWindowOverlap_plus[i,]$V8,collapse = "|")
    countingWindowOverlap_plus[i,]$utrStr_end = paste(countingWindowOverlap_plus[i,]$V9,collapse = "|")
    transcriptsInRefSeq = ensemblAnnotations[is.element(ensemblAnnotations$V7, doesTranscriptExistInOtherFragment$V11),]
    transcriptsInRefSeq = transcriptsInRefSeq[which.min(transcriptsInRefSeq$V2),]
    countingWindowOverlap_minus[i,]$transcriptName = transcriptsInRefSeq$V7

  }else{
    conuntingWindowsConsidered = transcriptToFragment_plus[is.element(transcriptToFragment_plus$V4,doesTranscriptExistInOtherFragment$V4),]
    #conuntingWindowsConsidered = countingWindowOverlap_plus[is.element(countingWindowOverlap_plus$V10,doesTranscriptExistInOtherFragment$V4),]
    conuntingWindowsConsidered = conuntingWindowsConsidered[!duplicated(conuntingWindowsConsidered),]
    conuntingWindowsConsidered =as.data.frame( conuntingWindowsConsidered %>% group_by(V4) %>% mutate(collapsed = paste(V11, collapse = "|")))
    a = conuntingWindowsConsidered[which.min(conuntingWindowsConsidered$V2),]
    conuntingWindowsConsidered = conuntingWindowsConsidered[which(conuntingWindowsConsidered$V11 == a$V11),]
    #conuntingWindowsConsidered$V11 = conuntingWindowsConsidered$collapsed
    conuntingWindowsConsidered$collapsed = NULL
    conuntingWindowsConsidered =  conuntingWindowsConsidered[which(conuntingWindowsConsidered$V8 == countingWindowOverlap_plus[i,]$V4),]
    countingWindowOverlap_plus[i,]$utrStr_chr = paste(conuntingWindowsConsidered$V5,collapse = "|")
    countingWindowOverlap_plus[i,]$utrStr_start = paste(conuntingWindowsConsidered$V6,collapse = "|")
    countingWindowOverlap_plus[i,]$utrStr_end = paste(conuntingWindowsConsidered$V7,collapse = "|")
    countingWindowOverlap_plus[i,]$V2 = min(conuntingWindowsConsidered$V2)
    countingWindowOverlap_plus[i,]$transcriptName = conuntingWindowsConsidered$V11[1]
    
    cat(i)
  }
  print(i)
  
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  


inaccurateAnnotationStarts = rbind(countingWindowOverlap_plus,countingWindowOverlap_minus)

write.table(inaccurateAnnotationStarts,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/inaccurateAnnotationStarts.bed",sep="\t",row.names = F,col.names = F,quote = F)

inaccurateAnnotationStarts_rearranged = cbind.data.frame(inaccurateAnnotationStarts[,c("V1","V2","V3","V4","V5","V6","utrStr_chr","utrStr_start","utrStr_end","transcriptName")])





write.table(inaccurateAnnotationStarts_rearranged,"/Users/pooja.bhat/Dropbox/UTRannotation/mESC/currentAnnotation_withThomas/annotationStarts//inaccurateAnnotationStarts_includingEnsembl.bed",sep="\t",row.names = F,col.names = F,quote = F)
inaccurateAnnotationStarts_rearranged = inaccurateAnnotationStarts_rearranged[is.finite(rowSums(inaccurateAnnotationStarts_rearranged[,c(2,3)])),]
write.table(inaccurateAnnotationStarts_rearranged,"/Users/pooja.bhat/Dropbox/UTRannotation/mESC/currentAnnotation_withThomas/annotationStarts//inaccurateAnnotationStarts_includingEnsembl_igv.bed",sep="\t",row.names = F,col.names = F,quote = F)



##### getting the starts now for overlapping ends

#refSeqAnnotation = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed",stringsAsFactors = F)
ensemblAnnotation = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/proteinCoding_annotatedUTRs.bed",stringsAsFactors = F)
refSeqAnnotation = rbind(ensemblAnnotation)


overlappingCountingWindows_plus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/overlappingCountingWindows_refSeq_plus.bed",stringsAsFactors = F)
overlappingCountingWindows_plus$V2 = overlappingCountingWindows_plus$V7
overlappingCountingWindows_plus$V3 = overlappingCountingWindows_plus$V8

overlappingCountingWindows_minus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/overlappingCountingWindows_refSeq_minus.bed",stringsAsFactors = F)
overlappingCountingWindows_minus$V2 = overlappingCountingWindows_minus$V7
overlappingCountingWindows_minus$V3 = overlappingCountingWindows_minus$V8

## checking if a counring window overlaps with more than one transcript

### splitting by the counting window

overlappingCountingWindows_plus_split = split(overlappingCountingWindows_plus,f = overlappingCountingWindows_plus$V9,drop = T)

for(i in 1:length(overlappingCountingWindows_plus_split)){
  a = which.min((overlappingCountingWindows_plus_split[[i]]$V17))
  ### going into the refSeq annotation and checking for the transcript name : if it is presrnt >1 times it is a duplication or UTR exon
  annotation_split = refSeqAnnotation[which(refSeqAnnotation$V7 == overlappingCountingWindows_plus_split[[i]]$V16[a]),]
  overlappingCountingWindows_plus_split[[i]]$V2 = min(annotation_split$V2)
  overlappingCountingWindows_plus_split[[i]]$utrStr_chr = paste(annotation_split$V1,collapse = "|")
  overlappingCountingWindows_plus_split[[i]]$utrStr_start = paste(annotation_split$V2,collapse = "|")
  overlappingCountingWindows_plus_split[[i]]$utrStr_end = paste(annotation_split$V3,collapse = "|")
  overlappingCountingWindows_plus_split[[i]]$V16 = paste(overlappingCountingWindows_plus_split[[i]]$V16,collapse = "|")
  overlappingCountingWindows_plus_split[[i]] = overlappingCountingWindows_plus_split[[i]][1,]
  overlappingCountingWindows_plus_split[[i]]$transcriptName = paste(annotation_split$V7,collapse = "|")
  
}

overlappingCountingWindows_plus = do.call(rbind, overlappingCountingWindows_plus_split)



overlappingCountingWindows_minus_split = split(overlappingCountingWindows_minus,f = overlappingCountingWindows_minus$V9,drop = T)

for(i in 1:length(overlappingCountingWindows_minus_split)){
  a = which.max((overlappingCountingWindows_minus_split[[i]]$V18))
  annotation_split = refSeqAnnotation[which(refSeqAnnotation$V7 == overlappingCountingWindows_minus_split[[i]]$V16[a]),]
  overlappingCountingWindows_minus_split[[i]]$V3 = max(annotation_split$V3)
  overlappingCountingWindows_minus_split[[i]]$utrStr_chr = paste(annotation_split$V1,collapse = "|")
  overlappingCountingWindows_minus_split[[i]]$utrStr_start = paste(annotation_split$V2,collapse = "|")
  overlappingCountingWindows_minus_split[[i]]$utrStr_end = paste(annotation_split$V3,collapse = "|")
  overlappingCountingWindows_minus_split[[i]]$V16 = paste(overlappingCountingWindows_minus_split[[i]]$V16,collapse = "|")
  overlappingCountingWindows_minus_split[[i]] = overlappingCountingWindows_minus_split[[i]][1,]
  overlappingCountingWindows_minus_split[[i]]$transcriptName = paste(annotation_split$V7,collapse = "|")
}

overlappingCountingWindows_minus = do.call(rbind, overlappingCountingWindows_minus_split)

overlappingCountingWindows = rbind(overlappingCountingWindows_plus,overlappingCountingWindows_minus)
overlappingCountingWindows_rearranged = cbind.data.frame(overlappingCountingWindows[,c("V1","V2","V3","V4","V5","V6","utrStr_chr","utrStr_start","utrStr_end","transcriptName")])
write.table(overlappingCountingWindows_rearranged,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/overlappingCountinngWindows_final.bed",sep="\t",quote = F,row.names = F)

write.table(overlappingCountingWindows_rearranged,"/Users/pooja.bhat/Dropbox/UTRannotation/mESC/currentAnnotation_withThomas/annotationStarts/overlappingCountinngWindows_final_includingEnsembl.bed",sep="\t",quote = F,row.names = F)

totalEnsembl = rbind(overlappingCountingWindows_rearranged,inaccurateAnnotationStarts_rearranged)
write.table(totalEnsembl,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/starts_ensemblOverlappers.bed",quote = F,row.names = F,col.names = F,sep="\t")
##### getting in the intergenic counting windows. This file also consists of some counting windows that are added some transcripts that do not contain UTR annotations
### this will be removed now based on distance

refSeqAnnotation = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/refSeq_mm10_GRCm38_06-march-2017/processed/refSeq_mrna_utrsPresent.bed",stringsAsFactors = F)
ensemblAnnotation = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/incoming/3PrimeEndAnnotation_OLD/annotations_mm10/ensembl_mm10_Ensembl_Genes_87_06-march-2017/proteinCoding_annotatedUTRs.bed",stringsAsFactors = F)
allAnnotation = rbind(refSeqAnnotation,ensemblAnnotation)
allAnnotation_plus = allAnnotation[which(allAnnotation$V6 == "+"),]
allAnnotation_minus = allAnnotation[which(allAnnotation$V6 == "-"),]



intergenic_plus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_plus_complement_complement_IntergeniccountingWindows.bed",stringsAsFactors = F)
intergenic_minus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/refSeq_mrna_utrsPresent_minus_complement_complement_IntergeniccountingWindows.bed",stringsAsFactors = F)
overlappingNonUTRcontaining_minus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/overlappingWithNonUTRContaining_minus.bed")
overlappingNonUTRcontaining_plus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/includingEnsembl/overlappingWithNonUTRContaining_plus.bed")

overlappingNonUTRcontaining_minus = overlappingNonUTRcontaining_minus[!duplicated(overlappingNonUTRcontaining_minus),]
overlappingNonUTRcontaining_plus = overlappingNonUTRcontaining_plus[!duplicated(overlappingNonUTRcontaining_plus),]


overlappingNonUTRcontaining_plus$id = paste0(overlappingNonUTRcontaining_plus$V1,overlappingNonUTRcontaining_plus$V2,overlappingNonUTRcontaining_plus$V3,overlappingNonUTRcontaining_plus$V4,overlappingNonUTRcontaining_plus$V6)
overlappingNonUTRcontaining_minus$id = paste0(overlappingNonUTRcontaining_minus$V1,overlappingNonUTRcontaining_minus$V2,overlappingNonUTRcontaining_minus$V3,overlappingNonUTRcontaining_minus$V4,overlappingNonUTRcontaining_minus$V6)


intergenic_plus$id = paste0(intergenic_plus$V1,intergenic_plus$V2,intergenic_plus$V3,intergenic_plus$V4,intergenic_plus$V6)
intergenic_minus$id = paste0(intergenic_minus$V1,intergenic_minus$V2,intergenic_minus$V3,intergenic_minus$V4,intergenic_minus$V6)

intergenic_plus = intergenic_plus[!is.element(intergenic_plus$id,overlappingNonUTRcontaining_plus$id),]
intergenic_minus = intergenic_minus[!is.element(intergenic_minus$id,overlappingNonUTRcontaining_minus$id),]
intergenic_plus$id = NULL
intergenic_minus$id = NULL

intergenic_plus$utrStr_chr = NA
intergenic_plus$utrStr_start = NA
intergenic_plus$utrStr_end = NA
intergenic_plus$transcriptName=NA
for(i in 1:nrow(intergenic_plus)){
  a = allAnnotation_plus[which(allAnnotation_plus$V4 == intergenic_plus[i,]$V4),]
  chosenTranscripts = a[which(a$V3 == max(a$V3)),]
  chosenTranscripts = chosenTranscripts[which.min(chosenTranscripts$V2),]
  intergenic_plus[i,]$V2 = chosenTranscripts$V2
    # chosenTranscripts = which.max(a$V3)
  # intergenic_plus[i,]$V2 = a[which.max(a$V3),]$V2
  chosenId = chosenTranscripts$V7
  fromAnnotaiton = allAnnotation_plus[which(allAnnotation_plus$V7 == chosenId),]
  intergenic_plus$utrStr_chr[i] = paste0(fromAnnotaiton$V1,collapse = "|")
  intergenic_plus$utrStr_start[i] = paste0(fromAnnotaiton$V2,collapse="|")
  intergenic_plus$utrStr_end[i] = paste0(fromAnnotaiton$V3,collapse="|")
  intergenic_plus$transcriptName[i]=  paste0(chosenTranscripts$V7,collapse = "|")
  print(i)
}


### minus strand 


intergenic_minus$utrStr_chr = NA
intergenic_minus$utrStr_start = NA
intergenic_minus$utrStr_end = NA

for(i in 1:nrow(intergenic_minus)){
  a = allAnnotation_minus[which(allAnnotation_minus$V4 == intergenic_minus[i,]$V4),]
  intergenic_minus[i,]$V2 -  a$V2
  chosenTranscripts = a[which(a$V2 == min(a$V2)),]
  chosenTranscripts = chosenTranscripts[which.max(chosenTranscripts$V3),]
  intergenic_minus[i,]$V3 = chosenTranscripts$V3
  # chosenTranscripts = which.min(a$V2)
  # intergenic_minus[i,]$V3 = a[which.min(a$V2),]$V3
  chosenId = chosenTranscripts$V7
  fromAnnotaiton = allAnnotation_minus[which(allAnnotation_minus$V7 == chosenId),]
  intergenic_minus$utrStr_chr[i] = paste0(fromAnnotaiton$V1,collapse = "|")
  intergenic_minus$utrStr_start[i] = paste0(fromAnnotaiton$V2,collapse="|")
  intergenic_minus$utrStr_end[i] = paste0(fromAnnotaiton$V3,collapse="|")
  intergenic_minus$transcriptName[i]=  paste0(chosenTranscripts$V7,collapse = "|")
  print(i)
}
totalIntergenic = rbind(intergenic_plus,intergenic_minus)
#totalIntergenic  = totalIntergenic[-which((totalIntergenic$V3 - totalIntergenic$V2 )<=0),]
write.table(totalIntergenic,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/intergenicStarts.bed",sep="\t",quote = F,row.names = F,col.names = F)
