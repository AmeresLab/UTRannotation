#$ -S /bin/bash
#$ -q public.q
#$ -cwd
# -pe smp 6-24
#$ -l hostname=compute-*
#$ -v QUANT
#$ -v QUANT_ALIGN

#### this is a script that can be used to obtain polyA reads from multiple fastQ files. 
### a parameter file is required that has information about the file and one column with the index with the file 


# 8. Parse parameter file to get variables.
#module load bedtools
#module load  bam2fastq

#set variables
ADAPTER="AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG" #defined by libary prep
INPUT=$QUANT
PARAMETER="$QUANT_ALIGN/sampleInfo.txt"

#set dirs
OUTDIR=$QUANT_ALIGN
LOG=$QUANT_ALIGN/logs/
mkdir -p $LOG

#get input file (index)
COUNTER=0
for index in `cat $PARAMETER`; do #parameter file = list of files
    COUNTER=$((COUNTER+1))
    if [ $COUNTER = $SGE_TASK_ID ]; then
	break
    fi
done

#init files
touch "$OUTDIR"/logFile.txt ## initialise a log file
touch "$OUTDIR"/readStats.txt ## initialise a file to store your read stats at different processing steps
touch "$OUTDIR"/processingSteps.txt ### just another file to list the processing steps 

touch "$LOG"/stderr_"$index".txt
touch "$LOG"/stdo_"$index".txt


date >"$OUTDIR"/logFile.txt
 date >>"$LOG"/stderr_"$index".txt
 date >>"$LOG"/stdo_"$index".txt

#################
#ADAPTER trimming
#################
module load cutadapt >> "$LOG"/stdo_"$index".txt

cutadapt_version=$(cutadapt --version)
echo running cutadapt version "$cutadapt_version" >> "$LOG"/stdo_"$index".txt

## the step below trims the adapter : you should replace this based on your library prep method

initialrReads=$(zcat $INPUT/"$index" | echo $((`wc -l`/4)))
echo the number of reads in initial fastq is "$initialrReads" >>"$LOG"/stdo_"$index".txt

cutadapt -a $ADAPTER -o "$OUTDIR"/"$index"_trimmed.fastq  --trim-n $INPUT/"$index" 2>"$LOG"/stderr_"$index".txt 1>>"$LOG"/stdo_"$index".txt
adapterTrimming=$(cat  "$OUTDIR"/"$index"_trimmed.fastq | echo $((`wc -l`/4)))

echo the number of reads after adapter trimming is "$adapterTrimming" >>"$LOG"/stdo_"$index".txt


#################
#Quality control
#################
module load fastqc 

fastqc "$OUTDIR"/"$index"_trimmed.fastq 

#################
# remove 12 nt from the 5' end
#################
module load fastx-toolkit

echo fastx_trimmer to remove 12 nts from the 5 end of the reads >>"$LOG"/stdo_"$index".txt
echo version of fastx_trimmer used is Version: `which fastx_trimmer` >>"$LOG"/stdo_"$index".txt

fastx_trimmer -Q33 -f 13 -m 1 -i "$OUTDIR"/"$index"_trimmed.fastq  > "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq 2>"$LOG"/stderr_"$index".txt 

fivePrimeTrimmed=$(cat  "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq | echo $((`wc -l`/4)))

echo number of reads after 5 trimmig "$fivePrimeTrimmed" >>"$LOG"/stdo_"$index".txt

#################
# getting reads that have A's at the 3' end :
#################
echo "retaining reads that have >=5 As the the 3 end" >>"$LOG"/stdo_"$index".txt

#AAAAA$ could theoretically also happen in PHRED score
egrep -A2 -B1 'AAAAA$' "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq  | sed '/^--$/d' > "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq  2>"$LOG"/stderr_"$index".txt

polyAreads=$(cat "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq | echo $((`wc -l`/4)))
echo readd withs polyA "$polyAreads"  >>"$LOG"/stdo_"$index".txt

#################
# remove polyA
#################

echo remove the polyAs at the end of polyA reads >>"$LOG"/stdo_"$index".txt
echo this is done using cutadapat. At this step we also perform a size filteirng of minimum 18 nucleotides >>"$LOG"/stdo_"$index".txt

#cut super long polyA to avoid internal polyA cutting of 5 or more As
cutadapt --no-indels -m 18 -e 0 -a "A{1000}"  -o "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq  "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq 1>>"$LOG"/stdo_"$index".txt 2>"$LOG"/stderr_"$index".txt

finalReads=$(cat "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq | echo $((`wc -l`/4)))
echo final reads after size filtering "$finalReads"  >>"$LOG"/stdo_"$index".txt

#################
# write stats
#################

touch "$LOG"/preProcessingNumbers_"$index".txt

echo initialFile:"$initialrReads" >>"$LOG"/preProcessingNumbers_"$index".txt
echo adapterTrimmed:"$adapterTrimming" >>"$LOG"/preProcessingNumbers_"$index".txt
echo fivePrimeTrimming:"$fivePrimeTrimmed" >>"$LOG"/preProcessingNumbers_"$index".txt
echo polyAcontaining:"$polyAreads" >>"$LOG"/preProcessingNumbers_"$index".txt
echo finalFile:"$finalReads" >>"$LOG"/preProcessingNumbers_"$index".txt


echo the pre processing before mapping has been completed.  >>"$LOG"/stdo_"$index".txt


echo initialFile:"$initialrReads"  >>"$LOG"/stdo_"$index".txt
echo adapterTrimmed:"$adapterTrimming"  >>"$LOG"/stdo_"$index".txt
echo fivePrimeTrimming:"$fivePrimeTrimmed"  >>"$LOG"/stdo_"$index".txt
echo polyAcontaining:"$polyAreads"  >>"$LOG"/stdo_"$index".txt
echo finalFile:"$finalReads"  >>"$LOG"/stdo_"$index".txt


#################
# getting the read length distribution 
#################

awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq >"$OUTDIR"/"$index"_lengthDistribution.txt

