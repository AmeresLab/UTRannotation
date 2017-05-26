#! /usr/bin/ Rscript
args = commandArgs(trailingOnly=T)

library(ggplot2)
library(reshape)
library(GenomicFeatures)
options(scipen=999)
#peaks_plus = read.table("//groups/ameres/Pooja/Projects/wholePipeline/singleread100processing/42075_GGTATA_C9NM7ANXX_2_20160919B_20160919_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped_filtered_plusStrand_filtered.bed")
# peaks_minus = read.table("/groups/ameres/Pooja/Projects/wholePipeline/singleread100processing/42075_GGTATA_C9NM7ANXX_2_20160919B_20160919_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped_filtered_minusStrand_filtered.bed")
# peaks_plus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/getMaxOfpeak/plus_changedCoords.bed")
# peaks_minus = read.table("/Volumes/groups/ameres/Pooja/Projects/wholePipeline/getMaxOfpeak/minus_changedCoords.bed")


peaks_plus = read.table(args[1])
peaks_minus = read.table(args[2])
peaks_minus$V4=NULL
peaks_plus$V4=NULL
colnames(peaks_plus) = c("chromosome_name","X3_utr_start","X3_utr_end","count")
colnames(peaks_minus) = c("chromosome_name","X3_utr_start","X3_utr_end","count")

peaks_plus$strand = "1"
peaks_minus$strand = "-1"

peaks_total = rbind(peaks_minus,peaks_plus)

get120bps = function(UTRannotation){
  UTRannotation_positive = UTRannotation[which(UTRannotation[,5]=="1"| UTRannotation[,5]=="+"),]
  UTRannotation_negative = UTRannotation[which(UTRannotation[,5]=="-1"| UTRannotation[,5]=="-"),]
  
  ###### processing the UTRs on the positive strand first :
  
  UTRannotation_positive$start_original = UTRannotation_positive[,2]
  UTRannotation_positive$end_original = UTRannotation_positive[,3]
  
  UTRannotation_positive[,2] = UTRannotation_positive$end_original - 60
  UTRannotation_positive[,3] = UTRannotation_positive$end_original  + 60
  
  ###### processing the UTRs on the negative strand
  
  UTRannotation_negative$start_original = UTRannotation_negative[,2]
  UTRannotation_negative$end_original = UTRannotation_negative[,3]
  
  UTRannotation_negative[,2] = UTRannotation_negative$start_original - 60
  UTRannotation_negative[,3] = UTRannotation_negative$start_original + 60
  
  UTRannotation_total = rbind(UTRannotation_positive,UTRannotation_negative)
  
  
  return(UTRannotation_total)
}

# get120bps(UTRannotation = peaks_total)
# get120bps = function(UTRannotation){
#   UTRannotation_positive = UTRannotation[which(UTRannotation$strand=="1"),]
#   UTRannotation_negative = UTRannotation[which(UTRannotation$strand=="-1"),]
#   
#   ###### processing the UTRs on the positive strand first :
#   
#   UTRannotation_positive$start_original = UTRannotation_positive$X3_utr_start
#   UTRannotation_positive$end_original = UTRannotation_positive$X3_utr_end
#   
#   UTRannotation_positive$X3_utr_start = UTRannotation_positive$end_original - 60
#   UTRannotation_positive$X3_utr_end = UTRannotation_positive$end_original + 60
#   
#   ###### processing the UTRs on the negative strand
#   
#   UTRannotation_negative$start_original = UTRannotation_negative$X3_utr_start
#   UTRannotation_negative$end_original = UTRannotation_negative$X3_utr_end
#   
#   UTRannotation_negative$X3_utr_start = UTRannotation_negative$start_original - 60
#   UTRannotation_negative$X3_utr_end = UTRannotation_negative$start_original + 60
#   
#   
#   UTRannotation_positive_bed = cbind.data.frame(UTRannotation_positive$chromosome_name,UTRannotation_positive$X3_utr_start,UTRannotation_positive$X3_utr_end)
#   UTRannotation_positive_bed$orient = "forward"
#   UTRannotation_positive_bed$num = 1
#   UTRannotation_positive_bed$strand = "+"
#   UTRannotation_positive_bed = cbind.data.frame(UTRannotation_positive_bed,UTRannotation_positive$start_original,UTRannotation_positive$end_original,UTRannotation_positive$count)
#   colnames(UTRannotation_positive_bed) = c("chr","start","end","orient","num","strand","start_original","end_original","count")
#   
#   UTRannotation_negative_bed = cbind.data.frame(UTRannotation_negative$chromosome_name,UTRannotation_negative$X3_utr_start,UTRannotation_negative$X3_utr_end)
#   UTRannotation_negative_bed$orient = "reverse"
#   UTRannotation_negative_bed$num = -1
#   UTRannotation_negative_bed$strand = "-"
#   
#   UTRannotation_negative_bed = cbind.data.frame(UTRannotation_negative_bed, UTRannotation_negative$start_original,UTRannotation_negative$end_original, UTRannotation_negative$count)
#   
#   colnames(UTRannotation_negative_bed) = c("chr","start","end","orient","num","strand","start_original","end_original","count")
#   
#   UTRannotation_modified = rbind(UTRannotation_positive_bed,UTRannotation_negative_bed)
#   return(UTRannotation_modified)  
# }
# 

peaks_total_modified = get120bps(UTRannotation = peaks_total)

peaks_total_modified$strand[which(peaks_total_modified$strand == "1")]<-"+"
peaks_total_modified$strand[which(peaks_total_modified$strand == "-1")]<-"-"

peak_name = paste("peak",c(1:nrow(peaks_total_modified)),sep="_")
peaks_total_modified_rearranged = cbind.data.frame(peaks_total_modified$chromosome_name,peaks_total_modified$X3_utr_start,peaks_total_modified$X3_utr_end,peak_name,peaks_total_modified$count,peaks_total_modified$strand,peaks_total_modified$start_original,peaks_total_modified$end_original)
colnames(peaks_total_modified_rearranged) = c("chr","start","end","peak_name","count","strand","start_original","end_original")
write.table(peaks_total_modified_rearranged,args[3],sep="\t",quote = F,row.names = F,col.names = F)

