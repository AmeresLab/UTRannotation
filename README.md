# UTRannotation

Contains scripts for UTR annotation pipeline version 0.0.1
This pipeline can be used to refine the 3' end annotations for mRNA transcripts using QuantSeq 3' end sequencing datasets. 
The different steps for this include: 

  1.	Preprocessing - trimming and quality control of reads.
  2.	Identification of priming sites - finding positions in the genome where priming occurs
  3.	Differentiating internal priming evnets - finding a trhreshold based on UTR nucleotide signatures and genomic A content to distinguish internal priming events.
  4.	Identification of intergenic ends - using RNAseq to identify potentially unannotated 3' ends
  5.	Identification of high confidence 3’ ends - removing background and retining sites where most priming occurs
  6.	Creation of counting windows - creating counting blocks to count QuantSeq signal.
