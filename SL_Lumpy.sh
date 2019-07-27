#!/bin/bash -l

#SBATCH --job-name=lumpy-sv # Job name

#sbatch -N 1 -n 24 -p med -t 40:00:00 SL_Lumpy.sh  sample sample.R1.fastq sample.R2.fastq control control.R1.fastq control.R2.fastq

SAMPLE=$1
FILE=$2
FILE2=$3
CONTROL=$4

# modules
module load bio samtools bedtools speedseq perl
# KNOWN_SNPS must be in reference order
# REF must end in .fa or .fasta
REF="/home/jdyuzon/Pr_genome2017_phased/ND886_haps_alt.fa"

echo "speedseq"
# Align the data
RG="@RG\tID:id\tSM:"$SAMPLE"\tLB:lib"
echo $RG
/share/apps/speedseq-2017-10-23/bin/speedseq align -R "@RG\tID:id\tSM:"$SAMPLE"\tLB:lib" -o "$SAMPLE" \
$REF $FILE $FILE2

samtools view -b -o t.$SAMPLE.bam $SAMPLE.bam
samtools index t.$SAMPLE.bam

echo "samtools view"

samtools view -r id $SAMPLE.bam \
    | /share/apps/speedseq-2017-10-23/bin/pairend_distro.py \
    -r 101 \
    -X 4 \
    -N 10000 \
    -o $SAMPLE.lib1.histo


echo  "lumpy express"
##Run LUMPY Express on a tumor-normal pair
/share/apps/bio/bio/bin/lumpyexpress \
    -B t.$SAMPLE.bam,t.$CONTROL.bam \
    -S $SAMPLE.splitters.bam,$CONTROL.splitters.bam \
    -D $SAMPLE.discordants.bam,$CONTROL.discordants.bam \
    -o ${SAMPLE}_${CONTROL}.vcf

echo "lumpy"
lumpy -e -P -b -mw 4 -tt 0 \
-pe id:$SAMPLE,bam_file:$SAMPLE.discordants.bam,histo_file:$SAMPLE.lib1.histo,mean:500,stdev:50,read_length:101,min_non_overlap:101,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
-sr id:$SAMPLE,bam_file:$SAMPLE.splitters.bam,back_distance:10,weight:1,min_mapping_threshold:20 \
    > Lumpy.${SAMPLE}_${CONTROL}.bedpe

echo "svtyper"
/share/apps/speedseq-2017-10-23/bin/svtyper \
    -B $SAMPLE.bam \
    -S $SAMPLE.splitters.bam \
    -i ${SAMPLE}_${CONTROL}.vcf \
    > ${SAMPLE}_${CONTROL}.gt.vcf

echo "svtyper if truncated error"
/share/apps/speedseq-2017-10-23/bin/svtyper \
    -B t.$SAMPLE.bam \
    -S $SAMPLE.splitters.bam \
    -i ${SAMPLE}_${CONTROL}.vcf \
    > ${SAMPLE}_${CONTROL}.gt.vcf

echo "format bedpe to bed"
awk '{gsub(/^[ \t]+|[ \t]+$/,"");print}' Lumpy.${SAMPLE}_${CONTROL}.bedpe|grep -e '^contig_'|grep --invert-match -e 'INTERCHROM'|grep -e 'TYPE'|cut -f 1,2,6,11|awk '{print $1 "\t" $2 "\t" $3 "\t" sample "_" control "\t" $4}' sample=$SAMPLE control=$CONTROL >Lumpy.${SAMPLE}_${CONTROL}.tmp.bed #no interchrom
awk '{gsub(/^[ \t]+|[ \t]+$/,"");print}' Lumpy.${SAMPLE}_${CONTROL}.bedpe|grep -e '^contig_'|grep -e 'INTERCHROM'|cut -f 1,2,3,11|awk '{print $1 "\t" $2 "\t" $3 "\t" sample "_" control "\t" $4}' sample=$SAMPLE control=$CONTROL> Lumpy.${SAMPLE}_${CONTROL}.tmp2.bed #get interchrom, first chromposition
awk '{gsub(/^[ \t]+|[ \t]+$/,"");print}' Lumpy.${SAMPLE}_${CONTROL}.bedpe|grep -e '^contig_'|grep -e 'INTERCHROM'|cut -f 4,5,6,11|awk '{print $1 "\t" $2 "\t" $3 "\t" sample "_" control "\t" $4}' sample=$SAMPLE control=$CONTROL> Lumpy.${SAMPLE}_${CONTROL}.tmp3.bed #get interchrom,  second chromposition
cat Lumpy.${SAMPLE}_${CONTROL}.tmp3.bed Lumpy.${SAMPLE}_${CONTROL}.tmp2.bed Lumpy.${SAMPLE}_${CONTROL}.tmp.bed \
|sed 's/TYPE://g'|sed 's/INTERCHROM/BND/g'|sed 's/DELETION/DEL/g' |sed 's/DUPLICATION/DUP/g'|sed 's/INSERTION/INS/g'|sed 's/INVERSION/INV/g'|bedtools sort -i stdin -faidx $REF.sizes>Lumpy.${SAMPLE}_${CONTROL}.fix.bed

rm Lumpy.${SAMPLE}_${CONTROL}.tmp3.bed Lumpy.${SAMPLE}_${CONTROL}.tmp2.bed Lumpy.${SAMPLE}_${CONTROL}.tmp.bed


bedtools intersect -header -v -a ${SAMPLE}_${CONTROL}.gt.vcf -b ${CONTROL}_${CONTROL}.gt.vcf > ${SAMPLE}_${CONTROL}.tmp.gt.vcf
bedtools intersect -header -v -a ${SAMPLE}_${CONTROL}.tmp.gt.vcf -b Lumpy.${CONTROL}_${CONTROL}.fix.bed> ${SAMPLE}_${CONTROL}.tmp1.gt.vcf

tumor has to have 0.7 or more SupportingEvidence SU than control FORMAT field tumorSU >=1.5* controlSU
cat ${SAMPLE}_${CONTROL}.tmp1.gt.vcf| java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "((GEN[0].SU >= 1.5*GEN[1].SU)&(GEN[0].DP >= 10)&(GEN[0].SU >= 0.2*GEN[0].DP))">${SAMPLE}_${CONTROL}.flt.gt.vcf
rm ${SAMPLE}_${CONTROL}.tmp1.gt.vcf ${SAMPLE}_${CONTROL}.tmp.gt.vcf


bedtools sort -faidx $REF.fai -i Lumpy.${SAMPLE}_${CONTROL}.fix.bed |bedtools intersect -wa -wb -a stdin -b ${SAMPLE}_${CONTROL}.flt.gt.vcf|awk -F\; '{print $1}'|awk -F '\t' '{print $1"\t"$2"\t"$3"\t"$11$2$3}'|awk '{gsub(/\SVTYPE=/,"");}1'|bedtools sort -faidx $REF.fai -i stdin >${SAMPLE}_${CONTROL}.flt.gt.bed

