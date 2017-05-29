### reading in refSeq overlappingCounting windows 


refSeqOverlappers = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/starts_refSeqOverlappers.bed",stringsAsFactors = F,sep="\t")
ensemblOverlappers = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/starts_ensemblOverlappers.bed",stringsAsFactors = F,sep="\t")
intergenicOverlappers = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/intergenicStarts.bed",stringsAsFactors = F,sep="\t")


allStarts = rbind(refSeqOverlappers,ensemblOverlappers,intergenicOverlappers)
allStarts$id = paste0("I",c(1:nrow(allStarts)))
allStarts_plus = allStarts[which(allStarts$V6 == "+"),]
allStarts_minus = allStarts[which(allStarts$V6 == "-"),]

#### processing the plusStrand

grep("|",allStarts_plus$V7)
splitCoords = strsplit(x = allStarts_plus$V8,split = "|",fixed = T)
isSplit = unlist(lapply(splitCoords,function(x) length(x)))
allStarts_plus_multiple = allStarts_plus[which(isSplit>1),]
allStarts_plus_single = allStarts_plus[which(isSplit == 1),]

split_matrix = as.data.frame(matrix(data = NA,ncol = 8,nrow=1))
### I want to get the correct end - now the ends are based on refSeq annotations, i want to make it absed on identified ends

for(i in 1:nrow(allStarts_plus_multiple)){
  
  starts_tmp = strsplit(allStarts_plus_multiple[i,8],split = "|",fixed = T)
  ends_tmp = strsplit(allStarts_plus_multiple[i,9],split = "|",fixed = T)
  starts_ends = cbind.data.frame(as.numeric(starts_tmp[[1]]),as.numeric(ends_tmp[[1]]))
  colnames(starts_ends) = c("V1","V2")
  starts_ends = starts_ends[order(starts_ends$V2,decreasing = F),]
  starts_ends$V2[nrow(starts_ends)] = allStarts_plus_multiple[i,3]
  
  ### for the plus strand the end has to be the end of the counting window. 
  #names = paste(allStarts_plus_multiple[i,]$V4,c(1:length(ends_tmp[[1]])),sep=".")
  names = paste(allStarts_plus_multiple[i,]$V4)
  id = allStarts_plus_multiple[i,]$id
  chr = allStarts_plus_multiple[i,]$V1
  strand = allStarts_plus_multiple[i,]$V6
  transcriptNames = allStarts_plus_multiple[i,]$V10
  #id = allStarts_plus_multiple[i,]$id
  rearranged = cbind.data.frame(chr,starts_ends,names,0,strand,id,transcriptNames)
  colnames(rearranged) = paste0("V",c(1:8))
  split_matrix = rbind.data.frame(split_matrix,rearranged)
  
}

### similar for minus strand

splitCoords = strsplit(x = allStarts_minus$V8,split = "|",fixed = T)
isSplit = unlist(lapply(splitCoords,function(x) length(x)))
allStarts_minus_multiple = allStarts_minus[which(isSplit>1),]
allStarts_minus_single = allStarts_minus[which(isSplit == 1),]

### I want to get the correct start - now the start are based on refSeq annotations, i want to make it absed on identified starts

split_matrix_minus = as.data.frame(matrix(data = NA,ncol = 8,nrow=1))


for(i in 1:nrow(allStarts_minus_multiple)){
  
  starts_tmp = strsplit(allStarts_minus_multiple[i,8],split = "|",fixed = T)
  ends_tmp = strsplit(allStarts_minus_multiple[i,9],split = "|",fixed = T)
  starts_ends = cbind.data.frame(as.numeric(starts_tmp[[1]]),as.numeric(ends_tmp[[1]]))
  colnames(starts_ends) = c("V1","V2")
  starts_ends = starts_ends[order(starts_ends$V2,decreasing = F),]
  starts_ends$V1[1] = allStarts_minus_multiple[i,2]
  ### for the plus strand the end has to be the end of the counting window. 
  #names = paste(allStarts_minus_multiple[i,]$V4,c(1:length(ends_tmp[[1]])),sep=".")
  names = paste(allStarts_minus_multiple[i,]$V4)
  id = allStarts_minus_multiple[i,]$id
  chr = allStarts_minus_multiple[i,]$V1
  strand = allStarts_minus_multiple[i,]$V6
  transcriptNames = allStarts_minus_multiple[i,]$V10
 # id = allStarts_plus_multiple[i,]$id
  # for(i in 1:length(starts_tmp[[1]])){
  #   starts_numeric = as.numeric(starts_tmp[[1]])
  #   ends_numeric = as.numeric(ends_tmp[[1]])
  #   starts_ends = cbind.data.frame(starts_numeric,ends_numeric)
  # }
  
  rearranged = cbind.data.frame(chr,starts_ends,names,0,strand,id,transcriptNames)
  colnames(rearranged) = paste0("V",c(1:8))
  split_matrix_minus = rbind.data.frame(split_matrix_minus,rearranged)
  
}
split_matrix_minus = split_matrix_minus[-1,]
split_matrix = split_matrix[-1,]

singleBlocks = rbind(allStarts_minus_single,allStarts_plus_single)
singleBlocks = singleBlocks[complete.cases(singleBlocks),]
singleBlocks = singleBlocks[,c("V1","V2","V3","V4","V5","V6","id","V10")]
colnames(singleBlocks) = paste0("V",c(1:ncol(singleBlocks)))
multipleBlocks = rbind(split_matrix_minus,split_matrix)
total = rbind(singleBlocks,multipleBlocks)
difference_total = total$V3 - total$V2
total = total[which(difference_total>0),]
### while running getFasta i get this error :
##Feature (chr18:150248555-150248873) beyond the length of chr18 size (90702639 bp).  Skipping.
##Feature (chr10:134537255-134538105) beyond the length of chr10 size (130694993 bp).  Skipping.
## manually removing these two

#total = total[-which(total$V1 =="chr18" & total$V2 == "150248555"),]
#total = total[-which(total$V1 =="chr10" & total$V2 == "134537255"),]

write.table(total,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/getSeqsFromIntervals.bed",sep="\t",row.names = F,col.names = F,quote = F)

### run bedtools getfasta -s -fi /groups/ameres/bioinformatics/references/mmu/mm10/mmu_mm10_whole_genome.fa -bed getSeqsFromIntervals.bed >sequencesFrom_geteqsFromIntervals.bed


##### reading in the fasta sequences 

fasta_total = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/sequencesFrom_geteqsFromIntervals.bed",stringsAsFactors = F)
fasta_total = fasta_total[seq(2,nrow(fasta_total),by = 2),]
fasta_total = as.character(fasta_total)
total_includingFasta = cbind.data.frame(total,fasta_total,stringsAsFactors=F)
library(dplyr)

a = as.data.frame(total_includingFasta %>% group_by(V7) %>% mutate(collapsedFasta = paste0(fasta_total,collapse="")),stringsAsFactors=F)
total_includingFasta_short = a[,c("V7","collapsedFasta")]
total_includingFasta_short = total_includingFasta_short[!duplicated(total_includingFasta_short),]
mergedStarts_withFasta = merge(allStarts,total_includingFasta_short,by.x="id",by.y="V7")
id= mergedStarts_withFasta$id
mergedStarts_withFasta$id = NULL
mergedStarts_withFasta$id = id
write.table(mergedStarts_withFasta,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/startsWithFastaSequence.bed",sep="\t",quote = F,row.names = F,col.names = F)

mergedStarts_withFasta$collapsedFasta = NULL
write.table(mergedStarts_withFasta,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/startsWithoutSeq_igv.bed",sep="\t",quote = F,row.names = F,col.names = F)

sessionInfo()


###


mESCcountingWindows = read.table("/Volumes/groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/mESC/output/final90percent/allAnnotations.bed",sep="\t")
mESCcountingWindows$number = paste0("CW",c(1:nrow(mESCcountingWindows)))
mESCcountingWindows_plus = mESCcountingWindows[which(mESCcountingWindows$V6 == "+"),]
mESCcountingWindows_minus = mESCcountingWindows[which(mESCcountingWindows$V6 == "-"),]
mESCcountingWindows_plus$Uid = paste0(mESCcountingWindows_plus$V3,mESCcountingWindows_plus$V4,mESCcountingWindows_plus$V6)
mESCcountingWindows_minus$Uid = paste0(mESCcountingWindows_minus$V2,mESCcountingWindows_minus$V4,mESCcountingWindows_minus$V6)

mESCcountingWindows = rbind.data.frame(mESCcountingWindows_plus, mESCcountingWindows_minus)
write.table(mESCcountingWindows,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/countingWindowsWithId.bed",sep="\t",quote = F,col.names = F,row.names = F)
mESC_starts = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/startsWithFastaSequence.bed",stringsAsFactors = F)

mESC_starts_plus = mESC_starts[which(mESC_starts$V6 == "+"),]
mESC_starts_minus = mESC_starts[which(mESC_starts$V6 == "-"),]

mESC_starts_plus$Uid = paste0(mESC_starts_plus$V3,mESC_starts_plus$V4,mESC_starts_plus$V6)
mESC_starts_minus$Uid = paste0(mESC_starts_minus$V2,mESC_starts_minus$V4,mESC_starts_minus$V6)

mESC_starts = rbind.data.frame(mESC_starts_plus, mESC_starts_minus)

mESC_starts  = merge(mESC_starts,mESCcountingWindows,by="Uid",all.x=T)
mESC_starts = mESC_starts[,c(2:13,20)]
colnames(mESC_starts) = paste0("V",c(1:ncol(mESC_starts)))
write.table(mESC_starts,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/startsWithFastaSequence_idIncluded.bed",quote = F,row.names = F,col.names = F,sep="\t")
##### 

mESCbed = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/getSeqsFromIntervals.bed")
mESCbed  =merge(mESCbed,mESC_starts,by.x="V7",by.y="V12")
mESCbed = mESCbed[,c(1:8,20)]
mESCbed_rearranged = mESCbed[,c("V1.x","V2.x","V3.x","V4.x","V5.x","V6.x","V8.x","V7","V13")]
colnames(mESCbed_rearranged)= paste0("V",c(1:ncol(mESCbed_rearranged)))
write.table(mESCbed_rearranged,"/Volumes/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/finalStarts/getSeqsFromIntervals_idIncluded.bed",quote = F,row.names = F,col.names = F)

