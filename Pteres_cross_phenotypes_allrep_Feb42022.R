library(ggpubr)
library(tidyverse)
library(grid)
library(car)

####read in phenotype data
phenotypes<-read.csv("PterespopsonKombar.csv", header = TRUE)
####label column showing mean percent leaf coverage for all three reps
colnames(phenotypes)[8]<-"mean.percent.leaf.coverage"
####label the Net and Spot lesions
phenotypes$Net.blotch.form[phenotypes$Net.blotch.form=='N']<-'Net Lesion'
phenotypes$Net.blotch.form[phenotypes$Net.blotch.form=='S']<-'Spot Lesion'
####reorder factors of lesion type and cross type
phenotypes$Net.blotch.form <- factor(phenotypes$Net.blotch.form,
                                levels = c('Net Lesion','Spot Lesion'),ordered = TRUE)
phenotypes$Cross.Type <- factor(phenotypes$Cross.Type,
                                levels = c('Parents','Intra','Inter'),ordered = TRUE)
phenotypes$CrossSubspecies <- factor(phenotypes$CrossSubspecies,
                                levels = c('Parents','PTTxPTT','PTMxPTM','PTTxPTM'),ordered = TRUE)
####set highest percent leaf coverage as fitness=1 for all reps
phenotypes$percent.leaf.coverage<-
  phenotypes$percent.leaf.coverage/max(na.omit(phenotypes$percent.leaf.coverage))
####set highest percent leaf coverage as fitness=1 for average of all reps
phenotypes$mean.percent.leaf.coverage<-phenotypes$mean.percent.leaf.coverage/max(na.omit(phenotypes$mean.percent.leaf.coverage))


my_comparisons <- list(c("PTMxPTM","PTTxPTM"), 
                       c("PTTxPTT","PTTxPTM"), 
                       c("Parents","PTTxPTM"), 
                       c("PTTxPTT","PTMxPTM"), 
                       c("Parents","PTTxPTT"))

a<-ggboxplot(phenotypes, x = "CrossSubspecies", y = "percent.leaf.coverage", 
             palette =c("red","blue"),
             add = "jitter",add.params = list(color = "Net.blotch.form"))+ 
  facet_wrap(~Rep)+
  stat_compare_means(comparisons = my_comparisons, size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  xlab("")+ylab("Fitness (% Leaf Coverage)")+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) 




pdf("~/Dropbox/PteresGeneticIncompatibilitiesManuscript/SharedData/Phenotype_percentleaf_allrep.pdf")
plot(a)
graphics.off()

##### Test of homogeneity of variances #######
leveneTest(percent.leaf.coverage ~ as.factor(Rep), data = phenotypes)
#Df F value Pr(>F)
#group   3  2.0016 0.1127
#541   
fligner.test(percent.leaf.coverage ~ as.factor(Rep), data = phenotypes)
#data:  percent.leaf.coverage by as.factor(Rep)
#Fligner-Killeen:med chi-squared = 6.6358, df = 3, p-value = 0.08446


####Medians-percent.leaf.coverage
###Parents-Net
median(phenotypes$percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
             phenotypes$Cross.Type=='Parents']) #0.9

###Intraspecies-Net
median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                                       phenotypes$Cross.Type=='Intra'])) #0.6315789

###Interspecies-Net
median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Net Lesion' &
                                                       phenotypes$Cross.Type=='Inter'])) #0.3473684

###Parents-Spot
median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                                       phenotypes$Cross.Type=='Parents'])) #0.5263158

###Intraspecies-Spot
median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                                       phenotypes$Cross.Type=='Intra'])) #0.4421053

###Interspecies-Spot
median(na.omit(phenotypes$mean.percent.leaf.coverage[phenotypes$Net.blotch.form=='Spot Lesion' &
                                                       phenotypes$Cross.Type=='Inter'])) #0.2526316


###################################################################
##################### Control for Lesion Type #####################
###################################################################

####Pooled data (all reps)

a<-ggboxplot(phenotypes, x = "CrossSubspecies", y = "mean.percent.leaf.coverage", add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons, size=2, 
                     method.args = list(alternative = "greater"))+ # Add pairwise comparisons p-value
  xlab("")+ylab("Fitness (% Leaf Coverage)")+
  theme(legend.position="None", strip.text.x = element_text(color = "white"))+
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) 


pdf("~/Dropbox/PteresGeneticIncompatibilitiesManuscript/Figure2_Phenotype_percentleaf_pooled")
plot(a)
graphics.off()


