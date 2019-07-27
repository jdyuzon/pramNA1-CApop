#!/bin/bash -l
# NOTE the -l flag!
# If you need any help, please email help@cse.ucdavis.edu
# Name of the job - You'll probably want to customize this.
#SBATCH -J bwa
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o bwa-%j.output
#SBATCH -e bwa-%j.output
#SBATCH --partition=serial
#SBATCH --job-name=chimeragenome
# no -n here, the user is expected to provide that on the command line.
# The useful part of your job goes below
#chmod u+rwx SL_BWA_batch.sh

#6/4/2013
#sckPr102_L.fastq
module load samtools picardtools/1.107 R maven java GATK bamtools hmmcopy bcftools vcftools plink bwa bedtools

n="45"


echo "BWA"

entry1=$1 #~/Obj1Fastq/*fastq
entry2=$2
echo $entry1 $entry2
#contig="/home/jdyuzon/Pr_genome2017_phased/ND886_haps_alt.fa"
contig="/home/jdyuzon/Pr_genome2017_phased/falcon_consensus.fasta"
#for entry in *.fastq
	isolate=`echo $entry1|awk -F/ '{print $5}'|awk -F_ '{print $1}'` 
	echo $isolate
outputName1=bbduk${isolate}_L.fastq
outputName2=bbduk${isolate}_R.fastq
       echo $outputName1
       echo $outputName2

echo "remember to unzip and convert FASTQ quality scores K to J"


echo "BBDUK"
/share/apps/bbmap-37-50/bbduk.sh -Xmx1g in1=$entry1 in2=$entry2 out1=$outputName1 out2=$outputName2 ref=~/BioSoft/adaptersv2.fa ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=rl trimq=10


echo "BWA ALN"
bwa aln -n 0.04 -o 1 -e -1 -i 5 -d 10 -l 32 -k 2 -M 3 -O 11 -E 4 -R 30 -q 0 -B 0 \
        $contig $outputName1 >bwa${isolate}_L.sai
bwa aln -n 0.04 -o 1 -e -1 -i 5 -d 10 -l 32 -k 2 -M 3 -O 11 -E 4 -R 30 -q 0 -B 0 \
        $contig $outputName2 >bwa${isolate}_R.sai
bwa sampe $contig bwa${isolate}_L.sai bwa${isolate}_R.sai $outputName1 $outputName2 > bwa${isolate}.PE.sam

echo "Sort"
#for entry in bwaTKTK002A_L.sai  #*L.sai  #test
        bwaStem=`echo bwa${isolate}.PE.sam|awk -F. '{print $1}'` #bwaTKTK002A
	echo "bwaStem" $bwaStem
#samtools:SAM to BAM and Sort
samtools view -Sb -o bwa${isolate}.PE.bam bwa${isolate}.PE.sam
echo "samtools sort"  bwa${isolate}.PE.bam bwa${isolate}.merged.sorted.bam
samtools sort  -o bwa${isolate}.merged.sorted.bam -O bam -T bwa${isolate}.PE.bam.sorted bwa${isolate}.PE.bam

#HMMCopy and generating wig files #also generating normal wig files
bamtools index -in  bwa${isolate}.merged.sorted.bam 
/share/apps/hmmcopy-0.1.1/bin/readCounter bwa${isolate}.merged.sorted.bam > ${isolate}.wig


echo "Bsc Bcf"
#Bsc
        #first, bam to sam
#        samtools view -h -o ${bwaStem}.merged.sorted.sam ${bwaStem}.merged.sorted.bam
        outputName=${bwaStem}.merged.sorted.sam
        echo $outputName
        outputName1=${bwaStem}.merged.sorted_scALL
 #       awk -F "\t" 'substr($3,0,8)=="scaffold" {print $3, "\t",$4}' $outputName > $outputName1

  #      perl /home/jdyuzon/BioSoft/SuperScaff.pl $outputName1 |sort -nk 1,1> bsc${isolate}.seq

#Bcf
   #     samtools mpileup -uf /home/jdyuzon/Pr_genome2017_phased/falcon_consensus.fasta ${bwaStem}.merged.sorted.bam|bcftools view -bvcg - > $isolate.raw.bcf
    #    bcftools view $isolate.raw.bcf | /home/jdyuzon/BioSoft/vcfutils.pl varFilter -D100 > $isolate.flt.vcf



echo "MarkDuplicates"
#Mark Duplicates
java -jar /share/apps/picardtools-1.107/MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT \
    AS=true REMOVE_DUPLICATES=true I= ${bwaStem}.merged.sorted.bam \
    O=${bwaStem}_sorted.markdup.bam M=${bwaStem}_sorted.markdup.bam.metrics

#AddOrReplaceReadGroups
java -jar /share/apps/picardtools-1.107/AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT \
    I=${bwaStem}_sorted.markdup.bam O=${bwaStem}_sorted.rg.bam \
    RGID=$isolate RGLB=$isolate RGPL=illumina RGPU=run RGSM=$isolate

#Index BAM file for easy navigation
samtools index ${bwaStem}_sorted.rg.bam 


echo "samtools mpileup"
isolate_stem=${bwaStem#bwa}
#RealignerTargetCreator
java -jar /share/apps/GATK-3.3/target/GenomeAnalysisTK.jar -T RealignerTargetCreator \
     -nt 4 -R $contig -I ${bwaStem}_sorted.rg.bam \
     -o ${isolate_stem}.intervals

#IndelRealigner
java -jar /share/apps/GATK-3.3/target/GenomeAnalysisTK.jar -T IndelRealigner \
    -R $contig -I ${bwaStem}_sorted.rg.bam  \
    -targetIntervals ${isolate_stem}.intervals -o ${isolate_stem}.realigned.bam

#samtools index bwa${isolate}.merged.sorted.bam
samtools index bwa${isolate}_sorted.markdup.bam


