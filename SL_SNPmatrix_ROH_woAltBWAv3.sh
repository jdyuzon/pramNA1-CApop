#!/bin/bash -l
# NOTE the -l flag!
# If you need any help, please email help@cse.ucdavis.edu
# Name of the job - You'll probably want to customize this.
#SBATCH -J bwa
# Standard out and Standard Error output files with the job number in the name.
#SBATCH -o bwa-%j.output
#SBATCH -e bwa-%j.output
#SBATCH --partition=serial
#SBATCH --job-name=SNPwoAlt
# no -n here, the user is expected to provide that on the command line.
# The useful part of your job goes below
#chmod u+rwx SL_BWA_batch.sh

#6/4/2013
#sckPr102_L.fastq
module load samtools picardtools/1.107 R maven java GATK bamtools hmmcopy bcftools vcftools plink bwa/0.7.17.r1188 bedtools bcftools
#Genome="/home/jdyuzon/Pr_genome2017_phased/ND886_haps_alt.fa"
Genome="/home/jdyuzon/Pr_genome2017_phased/falcon_consensus.fasta"
n="45"


#For normHetSNPs.vcf to run TitanCNV

/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP  bwaPr1556Aug2017.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >Pr1556Aug2017_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP  bwaPr1556Jul2017.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >Pr1556Jul2017_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP  bwaPr1556Jun2015.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >Pr1556Jun2015_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP  bwaPr1556unusual.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >Pr1556unusual_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP  bwaPr1556Jun2015.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >Pr1556Jun2015_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP  bwaPr1556unusual.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >Pr1556unusual_het.vcf


/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP bwaTKTK001C.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >TKTK001C_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP bwaPr745.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >Pr745_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP bwaPR710.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >PR710_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP bwaND886Jul2017.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >ND886Jul2017_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP bwaND886Aug2017.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >ND886Aug2017_het.vcf
/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome  -I -D -S -g -t DP,AD,ADF,ADR,SP bwaND886Jun2016.merged.sorted.bam | bcftools call -f GQ -m | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "isHet( GEN[0] )" >ND886Jun2016_het.vcf

#Obtain GC content for bins (precomputed once for the same reference and bin size)
/share/apps/hmmcopy-0.1.1/bin/gcCounter $Genome > gc.wig
#Getting a Mappability File
/share/apps/hmmcopy-0.1.1/util/mappability/generateMap.pl -b $Genome
/share/apps/hmmcopy-0.1.1/util/mappability/generateMap.pl $Genome
#Obtain average mappability for bins (precomputed once for the same reference and bin size)
/share/apps/hmmcopy-0.1.1/bin/mapCounter $Genome.map.bw > map.wig


/share/apps/samtools-1.3.1/bin/samtools mpileup -Q15 -uf $Genome -g -t DP,AD,ADF,ADR,SP \
ND886Aug2017.realigned.bam \
Pr1556Aug2017.realigned.bam \
ND886Jul2017.realigned.bam \
Pr1556Jul2017.realigned.bam \
ND886Jun2016.realigned.bam \
Pr1556Jun2015.realigned.bam \
Pr1556unusual.realigned.bam \
PR710.realigned.bam \
TKTK001C.realigned.bam \
Pr745.realigned.bam \
BS2016-32.realigned.bam \
CS1.realigned.bam \
CS2.realigned.bam \
FOP340-5299.realigned.bam \
JLSP04-4325.realigned.bam \
MR176.realigned.bam \
Pr1612.realigned.bam \
Pr223.realigned.bam \
Pr237.realigned.bam \
Pr52.realigned.bam \
ANN14-997-L2.realigned.bam \
BS2016-10.realigned.bam \
CS26-T2-1.realigned.bam \
HMG2017-2.realigned.bam \
LTS-UMCA1.realigned.bam \
MMWD2015-74.realigned.bam \
MMWD2015-91.realigned.bam \
Pr120.realigned.bam \
Pr1620.realigned.bam \
Pr93.realigned.bam \
Pr106.realigned.bam \
BS2014-584.realigned.bam \
Pr218.realigned.bam \
Pr438.realigned.bam \
Pr451.realigned.bam \
Pr455.realigned.bam \
Pr458.realigned.bam \
Pr467.realigned.bam \
Pr472.realigned.bam \
Pr486.realigned.bam \
Pr1537.realigned.bam \
Pr1652.realigned.bam \
BS96.realigned.bam \
MK649a.realigned.bam \
MK79j.realigned.bam \
TKTK001A.realigned.bam \
Pr102.realigned.bam \
|bcftools call -f GQ -m -A > m.PR_initialSNPset.0.bcf
/share/apps/bcftools-1.2/bin/vcfutils.pl varFilter -D8000 -a 0  m.PR_initialSNPset.0.bcf > m.sm.PR_initialSNPset.vcf

cat *tracks*bed |bedtools sort -faidx /home/jdyuzon/Pr_genome2017_phased/falcon_consensus.fasta.sizes  -i stdin | mergeBed -i stdin>CNV_Meta.mask.bed
/share/apps/bedtools-2.25.0/bin/subtractBed -a m.sm.PR_initialSNPset.vcf -b CNV_Meta.mask.bed -header >PR_initialSNPset.vcf

bgzip -c PR_initialSNPset.vcf >PR_initialSNPset.vcf.gz
tabix PR_initialSNPset.vcf.gz
vcf-contrast -n +ND886Aug2017 -ND886Jul2017,ND886Jun2016 PR_initialSNPset.vcf.gz > ND886filter.vcf
vcf-contrast -n +Pr1556Aug2017 -Pr1556Jul2017,Pr1556Jun2015,Pr1556unusual PR_initialSNPset.vcf.gz > Pr1556filter.vcf
bgzip Pr1556filter.vcf
bgzip ND886filter.vcf
tabix  Pr1556filter.vcf.gz
tabix ND886filter.vcf.gz
vcf-merge Pr1556filter.vcf.gz   ND886filter.vcf.gz  >Pr1556_ND886_filter.vcf
bedtools intersect -v -header -a PR_initialSNPset.vcf -b Pr1556_ND886_filter.vcf > PR_initialSNPset.min2.flt.vcf  
cat PR_initialSNPset.min2.flt.vcf | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter --inverse "( (GEN[ND886Jul2017].GT='./.')&(GEN[ND886Jun2016].GT='./.')&(GEN[Pr1556Jul2017].GT='./.' ) & (GEN[Pr1556Jun2015].GT='./.')& (GEN[Pr1556unusual].GT='./.'))" >  PR_initialSNPset.min2.flt2.vcf

### filters
cat PR_initialSNPset.min2.flt2.vcf | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "( QUAL >= 22 )" > m.PR_initialSNPset.QUAL.vcf

cat m.PR_initialSNPset.QUAL.vcf | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "( MQ >= 22 )" > m.PR_initialSNPset.MQ.0.vcf
java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar rmInfo m.PR_initialSNPset.MQ.0.vcf "DP, AD, ADF ADR" > m.PR_initialSNPset.MQ.vcf
/share/apps/vcftools-0.1.13/bin/vcftools --vcf m.PR_initialSNPset.MQ.vcf --recode --out m.DP --minDP 15 --maxDP 700 --remove-indels 


cat m.DP.recode.vcf | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "((GEN[*].GT='0/0')  & (GEN[*].AD[0]>=6))| ((GEN[*].GT='1/1')  & (GEN[*].AD[1]>=6))| ((GEN[*].GT='2/2')  & (GEN[*].AD[2]>=6))| ((GEN[*].GT='0/1')& (GEN[*].AD[0]>=6) &(GEN[*].AD[1]>=6)) |((GEN[*].GT='0/2')  & (GEN[*].AD[0]>=6)&(GEN[*].AD[2]>=6)) |((GEN[*].GT='1/2') & (GEN[*].AD[1]>=6) & (GEN[*].AD[2]>=6)) " > m.PR_initialSNPset.AD.0.vcf
#######genotype can't have more than 10 conflicting reads
bcftools filter -e '(GT="0/0" & (FMT/AD[1] + FMT/AD[2])>10) ' m.PR_initialSNPset.AD.0.vcf  | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.AD.00.vcf
bcftools filter -e '(GT="1/1" & (FMT/AD[0] + FMT/AD[2])>10) ' m.PR_initialSNPset.AD.00.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.AD.11.vcf
bcftools filter -e '(GT="2/2" & (FMT/AD[0] + FMT/AD[1])>10) ' m.PR_initialSNPset.AD.11.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.AD.22.vcf
bcftools filter -e '(GT="1/2" & (FMT/AD[0]>10)) ' m.PR_initialSNPset.AD.22.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.AD.12.vcf
bcftools filter -e '(GT="0/2" & (FMT/AD[1]>10)) ' m.PR_initialSNPset.AD.12.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.AD.02.vcf
#####genotype can't have 10% conflicting reads
bcftools filter -e ' (GT="0/0"  & ((FMT/AD[2]+FMT/AD[1])/(FMT/AD[0]+FMT/AD[1]+FMT/AD[2])>0.1)) ' m.PR_initialSNPset.AD.02.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.ADF.00.vcf
bcftools filter -e '  (GT="1/1"  & ((FMT/AD[2]+FMT/AD[0])/(FMT/AD[0]+FMT/AD[1]+FMT/AD[2])>0.1))' m.PR_initialSNPset.ADF.00.vcf| /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.ADF.11.vcf
bcftools filter -e '  (GT="2/2"  & ((FMT/AD[0]+FMT/AD[1])/(FMT/AD[0]+FMT/AD[1]+FMT/AD[2])>0.1)) ' m.PR_initialSNPset.ADF.11.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.ADF.22.vcf
bcftools filter -e '  (GT="1/2"  & ((FMT/AD[0])/(FMT/AD[0]+FMT/AD[1]+FMT/AD[2])>0.1)) ' m.PR_initialSNPset.ADF.22.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.ADF.12.vcf
bcftools filter -e '  (GT="0/2"  & ((FMT/AD[1])/(FMT/AD[0]+FMT/AD[1]+FMT/AD[2])>0.1)) ' m.PR_initialSNPset.ADF.12.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.ADF.02.vcf
bcftools filter -e '  (GT="0/1"  & ((FMT/AD[2])/(FMT/AD[0]+FMT/AD[1]+FMT/AD[2])>0.1)) ' m.PR_initialSNPset.ADF.02.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.PR_initialSNPset.ADF.01.vcf

cat m.PR_initialSNPset.ADF.01.vcf | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "( GEN[*].PL[*] = 0 ) " > m.PR_initialSNPset.PL.vcf
/share/apps/vcftools-0.1.13/bin/vcftools --vcf m.PR_initialSNPset.PL.vcf --recode --out m.hiconf.0 --minGQ 50
bcftools filter -e ' ((GT="0/0" | GT="1/1" | GT="2/2") & FMT/GQ<80 )'  m.hiconf.0.recode.vcf  |/share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > m.hiconf.recode.0.vcf ###????
/share/apps/vcftools-0.1.13/bin/vcftools --vcf m.hiconf.recode.0.vcf --recode --out m.hiconf.mis --max-missing-count 42
bcftools filter -i '(GT="1/1"| GT="0/0"| GT="2/2")' m.hiconf.mis.recode.vcf | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > SVfreeSNPs.poly.0.vcf
bcftools view SVfreeSNPs.poly.0.vcf --max-af 0.99 --exclude-types indels | /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > SVfreeSNPs.poly.1.vcf
bcftools filter -i '(GT="0/0"|GT="1/1"|GT="2/2")' SVfreeSNPs.poly.1.vcf| /share/apps/bcftools-1.2/bin/vcfutils.pl varFilter > SVfreeSNPs.poly.2.vcf 
echo $n " genotypes that occur 47 or 48 times/fixed will be removedy"
####this produces a splitstree with even number of branches
cat SVfreeSNPs.poly.2.vcf |java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "(countHom()>=countHet())&& (countHet()>=1)">SVfreeSNPs.min00.vcf 
######cat SVfreeSNPs.poly.2.vcf | java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter " ( ((( countHet(GEN[*].GT = '0/0') >=2 ) & ( countHet(GEN[*].GT = '0/0') <=$n ))  & ((((countHet(GEN[*].GT = '0/1') >=2) & (countHet(GEN[*].GT = '0/1') <=$n)))  | (((countHet(GEN[*].GT = '0/2') >=2) & (countHet(GEN[*].GT = '0/2') <=$n))) |  (((countHet(GEN[*].GT = '1/1') >=2) & (countHet(GEN[*].GT = '1/1') <=$n))) |  (((countHet(GEN[*].GT = '1/2') >=2) & (countHet(GEN[*].GT = '1/2')<=$n)))|(((countHet(GEN[*].GT = '2/2') >=2) & (countHet(GEN[*].GT = '2/2') <=$n))) ) )| (((countHet(GEN[*].GT = '0/1') >=2 ) & ( countHet(GEN[*].GT = '0/1') <=$n )) & (   (((countVariant(GEN[*].GT = '0/0') >=2) & (countHet(GEN[*].GT = '0/0') <=$n))) | (((countHet(GEN[*].GT = '0/2') >=2)& (countHet(GEN[*].GT = '0/2') <=$n))) |  (((countHet(GEN[*].GT = '1/1') >=2)& (countHet(GEN[*].GT = '1/1') <=$n))) |  (((countHet(GEN[*].GT = '1/2') >=2)& (countHet(GEN[*].GT = '1/2') <=$n))) |  (((countHet(GEN[*].GT = '2/2') >=2) & (countHet(GEN[*].GT = '2/2') <=$n))) ) ) | (((countHet(GEN[*].GT = '0/2') >=2 ) & ( countHet(GEN[*].GT = '0/2') <=$n )) & (   (((countVariant(GEN[*].GT = '0/0') >=2) & (countHet(GEN[*].GT = '0/0') <=$n)))| (((countHet(GEN[*].GT = '0/1') >=2) & (countHet(GEN[*].GT = '0/1') <=$n)))   |  (((countHet(GEN[*].GT = '1/1') >=2) & (countHet(GEN[*].GT = '1/1') <=$n))) |  (((countHet(GEN[*].GT = '1/2') >=2) & (countHet(GEN[*].GT = '1/2') <=$n))) |  (((countHet(GEN[*].GT = '2/2') >=2) & (countHet(GEN[*].GT = '2/2') <=$n))) ) ) | (((countHet(GEN[*].GT = '1/1') >=2 ) & (countHet(GEN[*].GT = '1/1') <=$n )) & (   (((countVariant(GEN[*].GT = '0/0') >=2) & (countHet(GEN[*].GT = '0/0') <=$n))) | (((countHet(GEN[*].GT = '0/1') >=2) & (countHet(GEN[*].GT = '0/1') <=$n)))  | (((countHet(GEN[*].GT = '0/2') >=2) & (countHet(GEN[*].GT = '0/2') <=$n))) | (((countHet(GEN[*].GT = '1/2') >=2) & (countHet(GEN[*].GT = '1/2') <=$n))) |  (((countHet(GEN[*].GT = '2/2') >=2) & (countHet(GEN[*].GT = '2/2') <=$n))) ) )  | (((countHet(GEN[*].GT = '1/2') >=2 ) & ( countHet(GEN[*].GT = '1/2') <=$n )) & (   (((countVariant(GEN[*].GT = '0/0') >=2) & (countHet(GEN[*].GT = '0/0') <=$n))) | (((countHet(GEN[*].GT = '0/1') >=2) & (countHet(GEN[*].GT = '0/1') <=$n)))  | (((countHet(GEN[*].GT = '0/2') >=2) & (countHet(GEN[*].GT = '0/2') <=$n))) |  (((countHet(GEN[*].GT = '1/1') >=2) & (countHet(GEN[*].GT = '1/1') <=$n)))  |   (((countHet(GEN[*].GT = '2/2') >=2) & (countHet(GEN[*].GT = '2/2') <=$n))) ) ) | (((countHet(GEN[*].GT = '2/2') >=2 ) & ( countHet(GEN[*].GT = '2/2') <=$n )) & (   (((countVariant(GEN[*].GT = '0/0') >=2) & (countHet(GEN[*].GT = '0/0') <=$n))) | (((countHet(GEN[*].GT = '0/1') >=2) & (countHet(GEN[*].GT = '0/1') <=$n)))  | (((countHet(GEN[*].GT = '0/2') >=2) & (countHet(GEN[*].GT = '0/2') <=$n)))  |  (((countHet(GEN[*].GT = '1/1') >=2) & (countHet(GEN[*].GT = '1/1') <=$n)))  |  (((countHet(GEN[*].GT = '1/2') >=2) & (countHet(GEN[*].GT = '1/2') <=$n)))  ))  )" > SVfreeSNPs.min0.vcf
##/share/apps/vcftools-0.1.13/bin/vcftools --vcf SVfreeSNPs.min0.vcf --recode --out SVfreeSNPs.min2 --max-maf 0.48
######/share/apps/vcftools-0.1.13/bin/vcftools --vcf SVfreeSNPs.min2.recode.vcf --012 --out SVfreeSNPs.matrix
##only ~80SNPs
###cat m.hiconf.mis.recode.vcf |java -jar /share/apps/bio/bio/share/snpsift/SnpSift.jar filter "(countHom()>=countHet())&& (countHet()>=2) ">test2.vcf

module load bedtools bcftools

bedtools subtract  -a SVfreeSNPs.min00.vcf -b SNPstooclose.bed -header >SVfreeSNPs.min0.vcf 

plink --vcf SVfreeSNPs.min0.vcf  --recode --out SVfreeSNPs.min2.flt --allow-extra-chr
tr -d '[:blank:]' < SVfreeSNPs.min2.flt.ped >test.ped
sed 's/^/>/' test.ped|sed -e $'s/000-9/\\\n/g' >SVfreeSNPs.min2.flt.matrix.ped
sed -i '/>/! s/0/-/g' SVfreeSNPs.min2.flt.matrix.ped
sed -e '/>/ s/-//g' SVfreeSNPs.min2.flt.matrix.ped > tmp.ped
##perl convert.sh tmp.ped > SVfreeSNPs.min2.raxml.fasta
##raxmlHPC -x 12345 -p 12345 -# autoMRE -m GTRGAMMA -n TEST -s SVfreeSNPs.min2.raxml.fasta
##echo "For Splitstree SVfreeSNPs.min2.flt.matrix.ped"
##echo "For RAxML SVfreeSNPs.min2.raxml.fasta"

sed -i -e 's/>ANN14-997-L2ANN14-997-L2/>ANN14-997-L2_Sonoma_2007_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>BS2014-584BS2014-584/>BS2014-584_Monterey_2014_Notholithocarpus-densiflorus_bark/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>BS2016-10BS2016-10/>BS2016-10_Monterey_2016_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>BS2016-32BS2016-32/>BS2016-32_Monterey_2016_Notholithocarpus-densiflorus_twig/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>BS96BS96/>BS96_Monterey_2003_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>CS1CS1/>CS1_San-Mateo_2016_Notholithocarpus-densiflorus_leaves/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>CS2CS2/>CS2_San-Mateo_2016_Umbellularia-californica_leaves/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>CS26-T2-1CS26-T2-1/>CS26-T2-1_Sonoma_2008_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>FOP340-5299FOP340-5299/>FOP340-5299_Sonoma_2010_Umbellularia-californica_leaves/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>HMG2017-2HMG2017-2/>HMG2017-2_Santa-Cruz_2017_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>JLSP04-4325JLSP04-4325/>JLSP04-4325_Sonoma_2010_Umbellularia-californica_leaves/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>LTS-UMCA1LTS-UMCA1/>LTS-UMCA1_Santa-Cruz_2017_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>TKTK001CTKTK001C/>MK548_San-Mateo_2008_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>MK649aMK649a/>MK649a_San-Mateo_2008_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>TKTK001ATKTK001A/>MK649b_San-Mateo_2008_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>MK79jMK79j/>MK79j_San-Mateo_2008_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>MMWD2015-74MMWD2015-74/>MMWD2015-74_Marin_2015_Notholithocarpus-densiflorus_twig/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>MMWD2015-91MMWD2015-91/>MMWD2015-91_Marin_2015_Notholithocarpus-densiflorus_twig/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>MR176MR176/>MR176_Marin_2005_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>ND886Aug2017ND886Aug2017/>ND886Aug2017_Marin_2004_Camellia-sp_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>ND886Jul2017ND886Jul2017/>ND886Jul2017_Marin_2004_Camellia-sp_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>ND886Jun2016ND886Jun2016/>ND886Jun2016_Marin_2004_Camellia-sp_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr52Pr52/>Pr52_Santa-Cruz_2001_Rhododendron_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr93Pr93/>Pr93_Santa-Cruz_2001_Rhododendron-sp_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr102Pr102/>Pr102_Marin_2001_Quercus_bark/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr106Pr106/>Pr106_Sonoma_2001_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr120Pr120/>Pr120_Mendocino_2001_Notholithocarpus-densiflorus_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr218Pr218/>Pr218_Sonoma_2002_Frangula-californica_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr223Pr223/>Pr223_Humboldt_2002_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr237Pr237/>Pr237_Monterey_2002_Lysimachia-latifolia_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr438Pr438/>Pr438_Mendocino_2006_Arbutus-menziesii_stem/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr451Pr451/>Pr451_Alameda_2004_Sequoia-sempervirens_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr455Pr455/>Pr455_Sonoma_2005_Osmorhiza-berteroi_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr458Pr458/>Pr458_Sonoma_2005_Adiantum-jordanii_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr467Pr467/>Pr467_Mendocino_2006_Corylus-cornuta_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr472Pr472/>Pr472_Marin_2006_Choisya-ternata_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr486Pr486/>Pr486_Yolo_2006_Camellia-sp_unknown/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>PR710PR710/>PR710_Santa-Clara_2009_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr745Pr745/>Pr745_Santa-Clara_2010_rainwater_water/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1537Pr1537/>Pr1537_Mendocino_2012_Abies-grandis_twig/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1556Aug2017Pr1556Aug2017/>Pr1556Aug2017_Santa-Clara_2012_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1556Jul2017Pr1556Jul2017/>Pr1556Jul2017_Santa-Clara_2012_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1556Jun2015Pr1556Jun2015/>Pr1556Jun2015_Santa-Clara_2012_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1556unusualPr1556unusual/>Pr1556unusual_Santa-Clara_2012_Umbellularia-californica_leaf/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1612Pr1612/>Pr1612_Mendocino_2013_Notholithocarpus-densiflorus_twig/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1620Pr1620/>Pr1620_Sonoma_2013_Notholithocarpus-densiflorus_twig/g' SVfreeSNPs.min2.flt.matrix.ped
sed -i -e 's/>Pr1652Pr1652/>Pr1652_Humboldt_2014_stream_water/g' SVfreeSNPs.min2.flt.matrix.ped


