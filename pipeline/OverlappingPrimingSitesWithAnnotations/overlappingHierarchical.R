#### removing duplicated in annotation

library(dplyr)
removeDuplicates_annotationGene = function(annotation){
  
  ################## collapsing anotations that are exactly overlapping : 
  
  annotation$id = paste(annotation$V1,annotation$V2,annotation$V3,annotation$V4,annotation$V6,sep="_")
  
  annotation = annotation %>% group_by(id) %>% mutate(collapsed = paste(V7, collapse = ","))
  
  ### removing the duplicated ids
  
  annotation =  annotation[!duplicated(annotation$id),]
  annotation$id = NULL
  annotation$V7 = annotation$collapsed
  annotation$collapsed = NULL
  return(as.data.frame(annotation, stringsAsFactors = F))
}

removeDuplicates_annotationAll = function(annotation){
  
  ################## collapsing anotations that are exactly overlapping : 
  
  annotation$id = paste(annotation$V1,annotation$V2,annotation$V3,annotation$V6,sep="_")
  
  annotation = annotation %>% group_by(id) %>% mutate(collapsed = paste(V7, collapse = ","))
  
  ### removing the duplicated ids
  
  annotation =  annotation[!duplicated(annotation$id),]
  annotation$id = NULL
  annotation$V7 = annotation$collapsed
  annotation$collapsed = NULL
  return(as.data.frame(annotation, stringsAsFactors = F))
}









### function for overlapping peaks of with annotations for the nuclotide profiles






#### overlapAnnotations - overlap between peak file and annotation. 




overlapAnnotations = function(peakFile,annotation){

#    peakFile$start1 <- apply(peakFile, 1, function(x) { if ((as.numeric(x["V3"]) - as.numeric(x["V2"])) > 0) {return(as.numeric(x["V2"]) + 1)} else {return(as.numeric(x["V2"]))}})
#    annotation$start1 <- apply(annotation, 1, function(x) { if ((as.numeric(x["V3"]) - as.numeric(x["V2"])) > 0) {return(as.numeric(x["V2"]) + 1)} else {return(as.numeric(x["V2"]))}})
    
  overlappingNonOverlapping = vector("list",2)
  names(overlappingNonOverlapping) = c("overlappers","nonOverlappers")
  peakFile_positive = peakFile[which(peakFile$V6 == "+"),]
  peakFile_negative= peakFile[which(peakFile$V6 == "-"),]
  peakFile_positive_grange = with(peakFile_positive,GRanges(V1, IRanges(start = V2,end = V3),strand = V6))
  peakFile_negative_grange = with(peakFile_negative,GRanges(V1,IRanges(start = V2,end = V3),strand = V6))

    #0->1based
    start(peakFile_positive_grange) <- start(peakFile_positive_grange) + 1
    start(peakFile_negative_grange) <- start(peakFile_negative_grange) + 1
    
  #### now creating the ranges for the annotation file 
  #annotation  = annotation[!duplicated(annotation[,c(1:3,6)]),]
  annotation_positive = annotation[which(annotation$V6 == "+"),]
  annotation_negative = annotation[which(annotation$V6 == "-"),]
  
  annotation_positive_granges =  with(annotation_positive,GRanges(V1, IRanges(start = V2,end = V3),strand = V6))
  annotation_negative_granges =  with(annotation_negative,GRanges(V1, IRanges(start = V2,end = V3),strand = V6))

    #0->1based
    start(annotation_positive_granges) <- start(annotation_positive_granges) + 1
    start(annotation_negative_granges) <- start(annotation_negative_granges) + 1

    
  overlappingPeak_positive = findOverlaps(peakFile_positive_grange,annotation_positive_granges,select="all")
  overlappingPeak_negative = findOverlaps(peakFile_negative_grange, annotation_negative_granges,select = "all")
  
  subjectQueryOverlap_positive = cbind(peakFile_positive[queryHits(overlappingPeak_positive),], annotation_positive[subjectHits(overlappingPeak_positive),])
#  subjectQueryOverlap_positive = subjectQueryOverlap_positive[with(subjectQueryOverlap_positive, order(subjectQueryOverlap_positive[,10], subjectQueryOverlap_positive[,12],decreasing = F)), ]
  
  subjectQueryOverlap_negative = cbind(peakFile_negative[queryHits(overlappingPeak_negative),], annotation_negative[subjectHits(overlappingPeak_negative),])
#  subjectQueryOverlap_negative = subjectQueryOverlap_negative[with(subjectQueryOverlap_negative, order(subjectQueryOverlap_negative[,10], subjectQueryOverlap_negative[,11],decreasing = T)), ]

    #remove start1 cols
    #subjectQueryOverlap_positive <- subjectQueryOverlap_positive[,colnames(subjectQueryOverlap_positive)!="start1"]
    #subjectQueryOverlap_negative <- subjectQueryOverlap_negative[,colnames(subjectQueryOverlap_negative)!="start1"]
    
  colnames(subjectQueryOverlap_positive) = c(paste0("V",c(1:6)),"sequences_polyApeaks_120bps","peakName","downstreamSeq","totalAs",paste0("V",c(7:13)))
  colnames(subjectQueryOverlap_negative) = c(paste0("V",c(1:6)),"sequences_polyApeaks_120bps","peakName","downstreamSeq","totalAs",paste0("V",c(7:13)))

  
    subjectQueryOverlap_positive  = subjectQueryOverlap_positive %>% group_by(V4) %>% mutate(collapsed = paste(V13,collapse="_"))
    subjectQueryOverlap_positive  = subjectQueryOverlap_positive %>% group_by(V4) %>% mutate(collapsedGene = paste(V10,collapse="_"))
  subjectQueryOverlap_positive = subjectQueryOverlap_positive[!duplicated(subjectQueryOverlap_positive$V4),]
  
  
    subjectQueryOverlap_negative  = subjectQueryOverlap_negative %>% group_by(V4) %>% mutate(collapsed = paste(V13,collapse="_"))
    subjectQueryOverlap_negative  = subjectQueryOverlap_negative %>% group_by(V4) %>% mutate(collapsedGene = paste(V10,collapse="_"))
  subjectQueryOverlap_negative = subjectQueryOverlap_negative[!duplicated(subjectQueryOverlap_negative$V4),]


    subjectQueryOverlap_positive <- as.data.frame(subjectQueryOverlap_positive, stringsAsFactors = F)
    subjectQueryOverlap_negative <- as.data.frame(subjectQueryOverlap_negative, stringsAsFactors = F)
    
 # subjectQueryOverlap_positive = subjectQueryOverlap_positive[!duplicated(subjectQueryOverlap_positive[,c(1:3)]),]
  # subjectQueryOverlap_negative = subjectQueryOverlap_negative[!duplicated(subjectQueryOverlap_negative[,c(1:3)]),]
  
  total_subjectQuery = rbind(subjectQueryOverlap_positive, subjectQueryOverlap_negative)
  
  ### non overlapping peaks 
  
  
  NonoverlappingPeak_positive_df = peakFile_positive[-queryHits(overlappingPeak_positive),]
  NonoverlappingPeak_negative_df = peakFile_negative[-queryHits(overlappingPeak_negative),]
  
  NonoverlappingPeaks = rbind(NonoverlappingPeak_positive_df, NonoverlappingPeak_negative_df)
  
  overlappingNonOverlapping[[1]] = total_subjectQuery
  overlappingNonOverlapping[[2]] = NonoverlappingPeaks 
  return(overlappingNonOverlapping)
}


### dividing based on the fraction of As. this function takes as input data frames produced in the previpus step 



getPASfractions_thresholded_binned = function(QueryData){
  nonOverlappinguery = QueryData
  outputList = vector("list",length = length(seq(0.12,0.96,0.12)))
  names(outputList) = c("0-0.12","0.12-0.24","0.24-0.36","0.36-0.48","0.48-0.60","0.60-0.72","0.72-0.84","0.84-1.00")
  threshold = c(seq(0.12,0.84,0.12), 1.01)
  for(i in 1: length(threshold)){
    query_threshold = QueryData[which(QueryData$totalAs< threshold[i]),]
    outputList[[i]] = query_threshold
    QueryData = QueryData[-which(QueryData$totalAs< threshold[i]),]
    
    cat(i)
  }
  return(outputList)
}


### checking for the presence of PAS signals

# 
# checkPAS = function(query_threshold)
# {
#   
#   query_threshold$upstreamSequences = substr(query_threshold$sequences_polyApeaks_120bps,start = 20,stop = 55)
#   query_threshold$upstreamSequences = toupper(query_threshold$upstreamSequences)
#   toMatch <- c("TATAAA", "AGTAAA", "AATACA","CATAAA","AATATA","GATAAA","AATGAA","AAGAAA","ACTAAA","AATAGA","AATAAT","AACAAA","ATTACA","ATTATA","AACAAG","AATAAG")
#   
#   ###### there ae if - else loops included here as, R throws an error, if the pattern is not found. 
#   
#   query_AATAAA = query_threshold[grep("AATAAA",query_threshold$upstreamSequences),]
#   if(nrow(query_AATAAA)>0){
#     query_minusAATAAA = query_threshold[-grep("AATAAA",query_threshold$upstreamSequences),]
#   }else{
#     #query_minusAATAAA = query_threshold
#     query_minusAATAAA = query_AATAAA
#   }
#   
#   query_ATTAAA=  query_minusAATAAA[grep("ATTAAA",  query_minusAATAAA$upstreamSequences),]
#   if(nrow(query_ATTAAA)>0){
#     query_minusATTAAA = query_minusAATAAA[-grep("ATTAAA",query_minusAATAAA$upstreamSequences),]
#   }else{
#     query_minusATTAAA = query_minusAATAAA
#   }
#   
#   
#   
#   query_APA = query_minusATTAAA[grep(paste(toMatch,collapse="|"),query_minusATTAAA$upstreamSequences),]
#   
#   
#   if(nrow(query_APA)>0){
#     query_minusAPA = query_minusATTAAA[-grep(paste(toMatch,collapse="|"),query_minusATTAAA$upstreamSequences),]
#   }else{
#     query_minusAPA = query_minusATTAAA
#   }
#   
#   
#   
#   queryStats = list(query_AATAAA,query_ATTAAA, query_APA, query_minusAPA)
#   names(queryStats) = c("AATAAA","ATTAAA","APA","noPAS")
#   return(queryStats)
#   
# }




checkPAS = function(query_threshold)
{
  
  query_threshold$upstreamSequences = substr(query_threshold$sequences_polyApeaks_120bps,start = 20,stop = 55)
  query_threshold$upstreamSequences = toupper(query_threshold$upstreamSequences)
  toMatch <- c("TATAAA", "AGTAAA", "AATACA","CATAAA","AATATA","GATAAA","AATGAA","AAGAAA","ACTAAA","AATAGA","AATAAT","AACAAA","ATTACA","ATTATA","AACAAG","AATAAG")
  
  query_AATAAA = filter(query_threshold,grepl("AATAAA",upstreamSequences))
  query_threshold = filter(query_threshold,!grepl("AATAAA",upstreamSequences))
  query_ATTAAA = filter(query_threshold,grepl("ATTAAA",upstreamSequences))
  query_threshold = filter(query_threshold,!grepl("ATTAAA",upstreamSequences))
  query_APA = filter( query_threshold,grepl(paste(toMatch,collapse="|"),upstreamSequences))
  query_minusAPA= filter( query_threshold,!grepl(paste(toMatch,collapse="|"),upstreamSequences))
  queryStats = list(query_AATAAA,query_ATTAAA, query_APA, query_minusAPA)
  names(queryStats) = c("AATAAA","ATTAAA","APA","noPAS")
  return(queryStats)
  
}


##### plotting the nucleotide profiles



gettingNucleotideComposition_binned = function(queryData){
  
  
  nucleotidePerseq = mapply(function(x) unlist(strsplit(toupper(x), "",fixed = T)),queryData) #convert all to upper case
  if(class(nucleotidePerseq)=="list"){
    names(nucleotidePerseq)=c(1:length(nucleotidePerseq))
    nucleotidePerseq=do.call(cbind,nucleotidePerseq)
  }
  nucleotideTablePerSeq = apply(nucleotidePerseq,1,table)
  
  if(class(nucleotideTablePerSeq)!="list"){ range_samples = c(-60:60)
  range_samples = range_samples[-which(range_samples == 0)]
  
  colnames(nucleotideTablePerSeq) = range_samples
  
  #no small acgt needed since toupper
  nucleotidePerseq_A = nucleotideTablePerSeq["A",]
  nucleotidePerseq_C = nucleotideTablePerSeq["C",]
  nucleotidePerseq_G = nucleotideTablePerSeq["G",]
  nucleotidePerseq_T = nucleotideTablePerSeq["T",]
  
  nucleotideDistribution = rbind(nucleotidePerseq_A , nucleotidePerseq_C , nucleotidePerseq_G , nucleotidePerseq_T )
  rownames(nucleotideDistribution) = c("A","C","G","T")
  nucleotideDistribution_sum = colSums(nucleotideDistribution)
  nucleotideDistribution = (nucleotideDistribution/nucleotideDistribution_sum)*100
  
  #nucleotidePerseq_A_plus_nonoverlapping = nucleotideDistribution[1,]
  #nucleotidePerseq_T_plus_nonoverlapping = nucleotideDistribution[4,]
  
  nucleotideDistribution.melt = melt(t(nucleotideDistribution))
  p = ggplot(nucleotideDistribution.melt,aes(x=X1,y=value,group=X2)) + geom_line(aes(col=X2),size=1.5) + theme_bw() + xlab("Nucleotide position (bps)") + ylab("Nucleotide composition (%)") + ggtitle("Nucleotide composition") + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm"))
  p = p + scale_x_continuous(breaks=seq(-60,60,5))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p = p+ scale_colour_manual(values = c("A" = "green4","C" = "blue","G" = "orange","T" = "red"))
  
  }else{a = t(nucleotidePerseq)
  range_samples = c(-60:60)
  range_samples = range_samples[-which(range_samples == 0)]
  Acontent = apply(a,MARGIN = 2,function(x) length(which(x=="A")))
  acontent = apply(a,MARGIN = 2,function(x) length(which(x=="a")))
  Tcontent = apply(a,MARGIN = 2,function(x) length(which(x=="T")))
  tcontent = apply(a,MARGIN = 2,function(x) length(which(x=="t")))
  Gcontent = apply(a,MARGIN = 2,function(x) length(which(x=="G")))
  gcontent = apply(a,MARGIN = 2,function(x) length(which(x=="g")))
  Ccontent = apply(a,MARGIN = 2,function(x) length(which(x=="C")))
  ccontent = apply(a,MARGIN = 2,function(x) length(which(x=="c")))
  
  Acontent_total = Acontent +acontent
  Tcontent_total = Tcontent + tcontent
  Ccontent_total = Ccontent + ccontent
  Gcontent_total = Gcontent + gcontent
  total_nucleotides = Acontent_total + Tcontent_total + Ccontent_total + Gcontent_total
  
  total_nucleotides = total_nucleotides[1]
  
  allNucleotides = list(Acontent_total,Tcontent_total,Ccontent_total,Gcontent_total)
  allNucleotides = lapply(allNucleotides, function(x) ((x*100)/total_nucleotides))
  names(allNucleotides) = c("A","T","C","G")
  allNucleotides.melt = melt(allNucleotides)
  allNucleotides.melt$position = range_samples
  
  p = ggplot(allNucleotides.melt,aes(x=position,y=value,group=L1))+geom_line(aes(col=L1),size=1.5)+ theme_bw() + xlab("Nucleotide position (bps)") + ylab("Nucleotide composition (%)") + ggtitle("Nucleotide composition") + theme(legend.title=element_blank(),axis.text.x = element_text(margin=margin(10,15,10,15,"pt")),axis.text.y = element_text(margin=margin(5,15,10,5,"pt")), axis.ticks.length = unit(-0.25 , "cm"))
  p = p + scale_x_continuous(breaks=seq(-60,60,5))+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
  p = p+ scale_colour_manual(values = c("A" = "green4","C" = "blue","G" = "orange","T" = "red"))
  
  }
  return(p)
  
}


### plotting




                     
