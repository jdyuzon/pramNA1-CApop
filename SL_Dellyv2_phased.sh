#!/bin/bash

#SBATCH --job-name=delly # Job name

#sbatch -N 1 -n 24 -p med -t 40:00:00 SL_Delly.sh  sample (isolate name) PREFIX SUFFIX control1 control2

begin=`date +%s`

SAMPLE=$1
PREFIX=$2
SUFFIX=$3

#"ND886Aug2017"
#"Pr1556Aug2017"
#"ND886Jul2017"
#"Pr1556Jul2017"
#"ND886Jun2016"
#"Pr1556Jun2015"
#"Pr1556unusual"


# modules
module load bio

# KNOWN_SNPS must be in reference order
# REF must end in .fa or .fasta
REF="/home/jdyuzon/Pr_genome2017_phased/ND886_haps_alt.fa"

for entry in ${SAMPLE}_*.tsv
do
Stem=`echo $entry|awk -F. '{print $1}'`  #bwaMK649a
echo $Stem
CONTROL=`echo $Stem|awk -F_ '{print $2}'`
echo $CONTROL

echo $PREFIX$SAMPLE$SUFFIX 
echo $PREFIX$CONTROL$SUFFIX
echo $PREFIX$CONTROL2$SUFFIX
delly call -t DEL -o ${SAMPLE}_${CONTROL}.DEL.bcf -g $REF $PREFIX$SAMPLE$SUFFIX \
$PREFIX$CONTROL$SUFFIX 

delly call -t DUP -o ${SAMPLE}_${CONTROL}.DUP.bcf -g $REF $PREFIX$SAMPLE$SUFFIX  \
$PREFIX$CONTROL$SUFFIX 

delly call -t INV -o ${SAMPLE}_${CONTROL}.INV.bcf -g $REF $PREFIX$SAMPLE$SUFFIX  \
$PREFIX$CONTROL$SUFFIX 

delly call -t BND -o ${SAMPLE}_${CONTROL}.BND.bcf -g $REF $PREFIX$SAMPLE$SUFFIX  \
$PREFIX$CONTROL$SUFFIX 

delly call -t INS -o ${SAMPLE}_${CONTROL}.INS.bcf -g $REF $PREFIX$SAMPLE$SUFFIX  \
$PREFIX$CONTROL$SUFFIX 

echo "Filter to get only PASS; Then produce two files with Precise only and Precise/Imprecise"


bcftools view ${SAMPLE}_${CONTROL}.DEL.bcf > ${SAMPLE}_${CONTROL}.DEL.vcf
bcftools view ${SAMPLE}_${CONTROL}.DUP.bcf > ${SAMPLE}_${CONTROL}.DUP.vcf
bcftools view ${SAMPLE}_${CONTROL}.INV.bcf > ${SAMPLE}_${CONTROL}.INV.vcf
bcftools view ${SAMPLE}_${CONTROL}.BND.bcf > ${SAMPLE}_${CONTROL}.BND.vcf
bcftools view ${SAMPLE}_${CONTROL}.INS.bcf > ${SAMPLE}_${CONTROL}.INS.vcf

bgzip -f ${SAMPLE}_${CONTROL}.DEL.vcf
bgzip -f ${SAMPLE}_${CONTROL}.DUP.vcf
bgzip -f ${SAMPLE}_${CONTROL}.INV.vcf
bgzip -f ${SAMPLE}_${CONTROL}.BND.vcf
bgzip -f ${SAMPLE}_${CONTROL}.INS.vcf

tabix -f ${SAMPLE}_${CONTROL}.DEL.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.DUP.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.INV.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.BND.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.INS.vcf.gz

echo "bcftools concat vcf"
bcftools concat -a -O b -o ${SAMPLE}_${CONTROL}.mergedSV.bcf ${SAMPLE}_${CONTROL}.DEL.vcf.gz ${SAMPLE}_${CONTROL}.DUP.vcf.gz ${SAMPLE}_${CONTROL}.INV.vcf.gz ${SAMPLE}_${CONTROL}.BND.vcf.gz ${SAMPLE}_${CONTROL}.INS.vcf.gz
echo "bcftools concat bcf"
bcftools concat -a -O b -o ${SAMPLE}_${CONTROL}.mergedSV.bcf ${SAMPLE}_${CONTROL}.DEL.bcf ${SAMPLE}_${CONTROL}.DUP.bcf ${SAMPLE}_${CONTROL}.INV.bcf ${SAMPLE}_${CONTROL}.BND.bcf ${SAMPLE}_${CONTROL}.INS.bcf
echo "bcftools view"
bcftools view -o ${SAMPLE}_${CONTROL}.mergedSV.vcf  ${SAMPLE}_${CONTROL}.mergedSV.bcf
echo "svprops make bed file"
~/BioSoft/svprops/src/svprops ${SAMPLE}_${CONTROL}.mergedSV.bcf | tail -n +2 | cut -f 1,2,4,5 > ${SAMPLE}_${CONTROL}.mergedSV.bed



delly filter -m 15 -a 0.1 -t DEL -f somatic -o ${SAMPLE}_${CONTROL}.DEL.somatic.bcf -s ${SAMPLE}_${CONTROL}.tsv ${SAMPLE}_${CONTROL}.DEL.bcf
delly filter -m 15 -a 0.1 -t DUP -f somatic -o ${SAMPLE}_${CONTROL}.DUP.somatic.bcf -s ${SAMPLE}_${CONTROL}.tsv ${SAMPLE}_${CONTROL}.DUP.bcf
delly filter -m 15 -a 0.1 -t INV -f somatic -o ${SAMPLE}_${CONTROL}.INV.somatic.bcf -s ${SAMPLE}_${CONTROL}.tsv ${SAMPLE}_${CONTROL}.INV.bcf
delly filter -m 15 -a 0.1 -t BND -f somatic -o ${SAMPLE}_${CONTROL}.BND.somatic.bcf -s ${SAMPLE}_${CONTROL}.tsv ${SAMPLE}_${CONTROL}.BND.bcf
delly filter -m 15 -a 0.1 -t INS -f somatic -o ${SAMPLE}_${CONTROL}.INS.somatic.bcf -s ${SAMPLE}_${CONTROL}.tsv ${SAMPLE}_${CONTROL}.INS.bcf

bcftools view ${SAMPLE}_${CONTROL}.DEL.somatic.bcf > ${SAMPLE}_${CONTROL}.DEL.somatic.vcf
bcftools view ${SAMPLE}_${CONTROL}.DUP.somatic.bcf > ${SAMPLE}_${CONTROL}.DUP.somatic.vcf
bcftools view ${SAMPLE}_${CONTROL}.INV.somatic.bcf > ${SAMPLE}_${CONTROL}.INV.somatic.vcf
bcftools view ${SAMPLE}_${CONTROL}.BND.somatic.bcf > ${SAMPLE}_${CONTROL}.BND.somatic.vcf
bcftools view ${SAMPLE}_${CONTROL}.INS.somatic.bcf > ${SAMPLE}_${CONTROL}.INS.somatic.vcf

bgzip -f ${SAMPLE}_${CONTROL}.DEL.somatic.vcf
bgzip -f ${SAMPLE}_${CONTROL}.DUP.somatic.vcf
bgzip -f ${SAMPLE}_${CONTROL}.INV.somatic.vcf
bgzip -f ${SAMPLE}_${CONTROL}.BND.somatic.vcf
bgzip -f ${SAMPLE}_${CONTROL}.INS.somatic.vcf

tabix -f ${SAMPLE}_${CONTROL}.DEL.somatic.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.DUP.somatic.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.INV.somatic.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.BND.somatic.vcf.gz
tabix -f ${SAMPLE}_${CONTROL}.INS.somatic.vcf.gz

echo "bcftools concat"
bcftools concat -a -O b -o ${SAMPLE}_${CONTROL}.dellySV.bcf ${SAMPLE}_${CONTROL}.DEL.somatic.vcf.gz ${SAMPLE}_${CONTROL}.DUP.somatic.vcf.gz ${SAMPLE}_${CONTROL}.INV.somatic.vcf.gz ${SAMPLE}_${CONTROL}.BND.somatic.vcf.gz ${SAMPLE}_${CONTROL}.INS.somatic.vcf.gz
echo "bcftools concat"
bcftools concat -a -O b -o ${SAMPLE}_${CONTROL}.dellySV.bcf ${SAMPLE}_${CONTROL}.DEL.somatic.bcf ${SAMPLE}_${CONTROL}.DUP.somatic.bcf ${SAMPLE}_${CONTROL}.INV.somatic.bcf ${SAMPLE}_${CONTROL}.BND.somatic.bcf ${SAMPLE}_${CONTROL}.INS.somatic.bcf
echo "bcftools view"
bcftools view -o ${SAMPLE}_${CONTROL}.dellySV.vcf  ${SAMPLE}_${CONTROL}.dellySV.bcf

echo "svprops make bed file"
~/BioSoft/svprops/src/svprops ${SAMPLE}_${CONTROL}.dellySV.bcf | tail -n +2 | cut -f 1,2,4,5 > ${SAMPLE}_${CONTROL}.dellySV.bed

echo "make bed file"
~/BioSoft/svprops/src/svprops ${SAMPLE}_${CONTROL}.BND.somatic.bcf| tail -n +2 | awk '{print $1"\t"($2)"\t"($2+1)"\t"$5"L";}'|perl -pi -e 's/-//g' > ${SAMPLE}_${CONTROL}.BND.somatic.bed
~/BioSoft/svprops/src/svprops ${SAMPLE}_${CONTROL}.BND.somatic.bcf| tail -n +2 |awk '{print $3"\t"($4)"\t"($4+1)"\t"$5"R";}'|perl -pi -e 's/-//g' >> ${SAMPLE}_${CONTROL}.BND.somatic.bed
awk -F '\t' '$1~/^#/ || !(int($2)==0)' ${SAMPLE}_${CONTROL}.dellySV.vcf>${SAMPLE}_${CONTROL}.dellySV.1.vcf
cat ${SAMPLE}_${CONTROL}.dellySV.1.vcf|java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "(SVTYPE[0]!='BND')">${SAMPLE}_${CONTROL}.dellytmp.bcf
~/BioSoft/svprops/src/svprops ${SAMPLE}_${CONTROL}.dellytmp.bcf| tail -n +2 | cut -f 1,2,4,5>${SAMPLE}_${CONTROL}.dellytmp.bed
echo ${SAMPLE}_${CONTROL}.dellytmp.bed
cat ${SAMPLE}_${CONTROL}.BND.somatic.bed ${SAMPLE}_${CONTROL}.dellytmp.bed|awk '{print $0, '\t' $3-$2}'|bedtools sort -faidx $REF.fai -i stdin > ${SAMPLE}_${CONTROL}.dellySV.bed
echo ${SAMPLE}_${CONTROL}.mergedSV.bed

echo "filter dellySV.vcf;none for controls"
cat ${SAMPLE}_${CONTROL}.dellySV.1.vcf| java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter \
"(MAPQ>=20)" >${SAMPLE}_${CONTROL}.flt0.delly.vcf
cat ${SAMPLE}_${CONTROL}.flt0.delly.vcf| java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter \
"((GEN[0].DV+GEN[0].DR+GEN[0].RR+GEN[0].RV) >= 10)">${SAMPLE}_${CONTROL}.flt1.delly.vcf   #read depth is at least 10 reads
cat ${SAMPLE}_${CONTROL}.flt1.delly.vcf|java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter \
"((GEN[0].DV+GEN[0].RV) >= 4)">${SAMPLE}_${CONTROL}.flt.delly.vcf
rm ${SAMPLE}_${CONTROL}.BND.somatic.bed ${SAMPLE}_${CONTROL}.dellytmp.bcf ${SAMPLE}_${CONTROL}.dellytmp.bed ${SAMPLE}_${CONTROL}.flt0.delly.vcf ${SAMPLE}_${CONTROL}.flt1.delly.vcf ${SAMPLE}_${CONTROL}.dellySV.1.vcf

done
