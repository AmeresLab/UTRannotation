cd /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/

## split into + and - strand 

awk '/+/ { print }' /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/ncRNAsTotal.bed | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/ncRNAsTotal_plus.bed
awk '/-/ { print }' /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/ncRNAsTotal.bed  | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts/ncRNAsTotal_minus.bed


##### separate the counting windows into different strands :


awk '/+/ { print }' /groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/mESC/output/final90percent/allAnnotations.bed | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_plus.bed
awk '/-/ { print }' /groups/bioinfo/shared/Ameres.GRP/Veronika.Herzog/SLAMannotation/mESC/output/final90percent/allAnnotations.bed | sort -k1,1 -k2,2n >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_minus.bed


ml bedtools

bedtools intersect -v -a /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_plus.bed -b /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//ncRNAsTotal_plus.bed -wo >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_plus_ncRNARemoved.bed
bedtools intersect -v -a /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_minus.bed -b /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//ncRNAsTotal_minus.bed -wo >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_minus_ncRNARemoved.bed

cat /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_plus_ncRNARemoved.bed /groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//allAnnotations_minus_ncRNARemoved.bed >/groups/ameres/Pooja/Projects/wholePipeline/pipeline_1/createStarts//countingWindows_ncRNAsRemoved.bed 
