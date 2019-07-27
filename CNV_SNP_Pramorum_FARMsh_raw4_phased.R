#!/usr/bin/Rscript

library(TitanCNA)
library(HMMcopy)

args <- commandArgs(trailingOnly = TRUE) 
Isolate <-as.character(args[1])
Normal <- as.character(args[2])
norm_reads.wig<-paste0(Normal, ".wig")
gc.wig<-paste0("gc.wig")
map.wig<-paste0("map.wig")
positions <-paste0(Normal,"_het.vcf")

bamFile <-paste0("bwa",Isolate, ".merged.sorted.bam")
bamIndex <-paste0("bwa",Isolate, ".merged.sorted.bam.bai")
tum_reads.wig<-paste0( Isolate, ".wig")

print(Isolate)

tum_reads.wig<-paste0(Isolate, ".wig")
print(tum_reads.wig)

#Correct readcount and estimate copy number
#tum_uncorrected_reads <- wigsToRangedData(tum_reads.wig,"gc.wig", "map.wig")
#norm_uncorrected_reads <- wigsToRangedData(norm_reads.wig,"gc.wig", "map.wig")
#tum_corrected_copy <- correctReadcount(tum_uncorrected_reads)
#norm_corrected_copy <- correctReadcount(norm_uncorrected_reads)
#Segmentation
#tum_corrected_copy$copy <- tum_corrected_copy$copy - norm_corrected_copy$copy
# Segmenting
#param <- HMMsegment(tum_corrected_copy, getparam = TRUE) # retrieve converged parameters via EM
#param$mu <- log(c(1, 1.4, 2, 2.7, 3, 4.5) / 2, 2)
#param$m <- param$mu
#segmented_copy <- HMMsegment(tum_corrected_copy, param) # perform segmentation via Viterbi
# Export to SEG format for CNAseq segmentation
#rangedDataToSeg(tum_corrected_copy, file = paste0(Isolate,"_tum_corrected_copy.seg"))########change so isolate# shows up

#filename<-paste0(Isolate, "plots.pdf")
#pdf(file=filename, onefile=TRUE) ####change so isolate # shows up
#plotBias(tum_corrected_copy)
#plotCorrection(tum_corrected_copy)
#plotSegments(tum_corrected_copy, segmented_copy)
#dev.off()

numClusters <-1 #2’ or higher treats the tumour data as being subclonal.

##CHANGE ONCE CALCULATE PERCENT OF CNV/LOH
params <- loadDefaultParameters(copyNumber=8,numberClonalClusters=numClusters) #maximum number of copies TitanCNV can consider
params$normalParams$n_0 <- 0   #If the sample has absolutely no normal contamination, then assign nParams$n_0 <- 0 and use argument normalEstimateMethod="fixed". #set initial normal proportion to 0.5

params$ploidyParams$phi_0 <- 2   #Though mapped to both haplotypes, best to set ploidy to 2 so LogRatio baseline is 2
				 #tried ploidy to 1 but plotClonalFrequency_JYedit sets LogRatio baseline to -1


#####################################################################
#####################################################################

## THE FOLLOWING SCRIPT ONLY WORKS WITH Rsamtools (>= 1.17.11)
library(Rsamtools)
library(VariantAnnotation)
#####################################################################
################ USERS - PLEASE MODIFY THESE PATHS #################
#####################################################################
tumbamFile <- paste0("bwa",Isolate, ".merged.sorted",".bam")
tumIndexFile <- paste(tumbamFile,".bai",sep = "")
vcfFile <- paste0(Normal,"_het.vcf")
outFile <- paste0(Isolate,"_",Normal,"_tumAlleleCounts.tsv")

#####################################################################
#################### LOAD HET POSITIONS VCF FILE ####################
#####################################################################
## read in vcf file of het positions
fl <- positions
vcf <- readVcf(fl, genome = "/home/jdyuzon/Pr_genome2017_phased/ND886_haps_alt.fa")

#####################################################################
####################### SETUP PILEUP PARAMETERS #####################
#####################################################################
## setup PileupParam using sequence read filters
pp <- PileupParam(min_base_quality = 10, min_mapq = 20, 
                  min_nucleotide_depth = 10, max_depth = 200, 
                  distinguish_strands = FALSE, 
                  distinguish_nucleotides = TRUE)
## setup the positions of interest to generate the pileup for
which <- GRanges(as.character(seqnames(vcf)), 
                 IRanges(start(ranges(vcf)), width = 1))
## setup addition BAM filters, such as excluding duplicate reads
sbp <- ScanBamParam(flag = scanBamFlag(isDuplicate = FALSE), which = which)

#####################################################################
########################## GENERATE PILEUP ##########################
#####################################################################
## generate pileup using function (Rsamtools >= 1.17.11)
## this step can take a while, # to get counts of tumour bam
tumbamObj <- BamFile(tumbamFile, index = tumIndexFile)
counts <- pileup(tumbamObj, scanBamParam = sbp,  pileupParam = pp)

## set of command to manipulate the "counts" data.frame of output
##     by pileup() such that multiple nucleotides are in a single
##     row rather than in multiple rows.
countsMerge <- xtabs(count ~ which_label + nucleotide, counts)
label <- do.call(rbind, strsplit(rownames(countsMerge), ":"))
posn <- do.call(rbind, strsplit(label[, 2],"-"))
#CHROM<-as.list(label[, 1])
countsMerge <- cbind(chr = label[, 1], position = posn[, 1], countsMerge)
countsMg<-as.data.frame(countsMerge)
#mode(countsMerge) <- "numeric"

#####################################################################
############### GET REFERENCE AND NON-REF READ COUNTS ###############
#####################################################################
## this block of code is used to match up the reference and 
##   non-reference nucleotide when assigning read counts #of tumour bam
##   final output data.frame is "countMat"
## setup output data.frame
v<-read.table(positions)
alt<-v$V5
countMat <- data.frame(chr = as.data.frame(seqnames(vcf))$value, #67499
                       position = posn[,1], #67499
                       ref = as.data.frame(fixed(vcf)["REF"])$REF, refCount = 0,  #81997
                       Nref = alt, NrefCount = 0, #67714
                       stringsAsFactors = FALSE)

## match rows with vcf positions of interest
countMat <- merge(countMat, countsMerge, by = c("chr","position"), 
                  sort = FALSE)



## assign the flattened table of nucleotide counts to ref, Nref
## note that non-reference (Nref) allele is sum of other bases
##    that is not matching the ref.
NT <- c("A", "T", "C", "G")
for (n in 1:length(NT)){	
  indRef <- countMat$ref == NT[n]
  countMat[indRef, "refCount"] <- as.vector(countMat[indRef, NT[n]])#####!!!!!!!:)
  foo<-countMat[indRef, NT[-n]]#####!!!!!!!:)
  bar <- data.frame(x=as.numeric(as.matrix(as.vector(foo[1]))), y=as.numeric(as.matrix(as.vector(foo[2]))), z=as.numeric(as.matrix(as.vector(foo[3]))), stringsAsFactors = FALSE) #####!!!!!!!:)
  rowSums(bar)#####!!!!!!!:)
  countMat[indRef, "NrefCount"] <- rowSums(bar)#####!!!!!!!:)
}
## remove "chr" string from chromosome
countMat$chr <- gsub("chr","",countMat$chr)	
## only use autosomes and sex chrs
#countMat <- countMat[countMat$chr %in% c(as.character(1:22),"X","Y"),]
## only use first 6 columns for TitanCNA
countMat <- countMat[,1:6]

#####################################################################
####################### OUTPUT TO TEXT FILE #########################
#####################################################################
## output text file will have the same format required by TitanCNA
write.table(countMat, file = outFile, row.names = FALSE, 
            col.names = TRUE, quote = FALSE, sep = "\t")

#####################################################################
#####################################################################
#####################################################################
loadAlleleCounts_JYedit <- function(inCounts, symmetric = TRUE, 
                                    genomeStyle = "NCBI", sep = "\t", header = TRUE) {
  if (is.character(inCounts)){
    #### LOAD INPUT READ COUNT DATA ####
    message("titan: Loading data ", inCounts)
    #data <- read.delim(inCounts, header = header, stringsAsFactors = FALSE, 
    #        sep = sep)
    intable <- read.delim(inCounts, header = header, stringsAsFactors = FALSE, 
                          sep = sep)
    #if (typeof(data[,2])!="integer" || typeof(data[,4])!="integer" || 
    #  typeof(data[,6])!="integer"){
    # stop("loadAlleleCounts: Input counts file format does not 
    #    match required specifications.")		
    # }
    if (typeof(intable[,2])!="integer" || typeof(intable[,4])!="integer" || 
        typeof(intable[,6])!="integer"){
      stop("loadAlleleCounts: Input counts file format does not 
           match required specifications.")		
    }
    #}else if (is.data.frame(inCounts)){  #inCounts is a data.frame
    #data <- inCounts
    }else if (is.data.frame(inCounts)){  #inCounts is a data.frame
      intable <- inCounts  
      #}else{
      # stop("loadAlleleCounts: Must provide a filename or data.frame 
      #     to inCounts")
      #}
    }else{
      stop("loadAlleleCounts: Must provide a filename or data.frame 
           to inCounts")
    }
  ## use GenomeInfoDb
  #require(GenomeInfoDb)
  # convert to desired genomeStyle and only include autosomes, sex chromosomes
  #data[, 1] <- setGenomeStyle(data[, 1], genomeStyle)
  
  ## sort chromosomes
  #indChr <- orderSeqlevels(as.character(data[, 1]), X.is.sexchrom = TRUE)
  indSc<-orderSeqlevels(as.character(intable[, 1]), X.is.sexchrom = TRUE)
  #data <- data[indChr, ]
  intable<-intable[indSc,]
  ## sort positions within each chr
  #for (x in unique(data[, 1])){
  # ind <- which(data[, 1] == x)
  #data[ind, ] <- data[ind[sort(data[ind, 2], index.return = TRUE)$ix], ]
  #}
  for (x in unique(intable[, 1])){
    ind <- which(intable[, 1] == x)
    intable[ind, ] <- intable[ind[sort(intable[ind, 2], index.return = TRUE)$ix], ]
  }
  
  #refOriginal <- as.numeric(data[, 4])
  refOriginal <- as.numeric(intable[, 4])
  #nonRef <- as.numeric(data[, 6])
  nonRef <- as.numeric(intable[, 6])
  tumDepth <- refOriginal + nonRef
  if (symmetric) {
    ref <- apply(cbind(refOriginal, nonRef), 1, max, na.rm = TRUE)
  } else {
    ref <- refOriginal
  }
  
  # return(list(chr = data[, 1], posn = data[, 2], ref = ref, 
  #            refOriginal = refOriginal, nonRef = nonRef, 
  #           tumDepth = tumDepth))
  return(list(chr = intable[, 1], posn = intable[, 2], ref = ref, 
              refOriginal = refOriginal, nonRef = nonRef, 
              tumDepth = tumDepth))
    }

#####################################################################
#####################################################################
#####################################################################
correctReadDepth_JYedit <- function(tumWig, normWig, gcWig, mapWig, 
                                    genomeStyle = "NCBI", targetedSequence = NULL) {
  #require(HMMcopy)
  
  message("Reading GC and mappability files")
  gc <- wigToRangedData(gcWig)
  map <- wigToRangedData(mapWig)
  
  ### LOAD TUMOUR AND NORMAL FILES ###
  message("Loading tumour file:", tumWig)
  tumour_reads <- wigToRangedData(tumWig)
  message("Loading normal file:", normWig)
  normal_reads <- wigToRangedData(normWig)
  
  ### set the genomeStyle: NCBI or UCSC
  #require(GenomeInfoDb)
  #names(gc) <- setGenomeStyle(names(gc), genomeStyle)
  #names(map) <- setGenomeStyle(names(map), genomeStyle)
  #names(tumour_reads) <- setGenomeStyle(names(tumour_reads), genomeStyle)
  #names(normal_reads) <- setGenomeStyle(names(normal_reads), genomeStyle)
  
  ### make sure tumour wig and gc/map wigs have same
  ### chromosomes
  gc <- gc[gc$space %in% tumour_reads$space, ]
  map <- map[map$space %in% tumour_reads$space, ]
  samplesize <- 50000
  
  ### for targeted sequencing (e.g.  exome capture),
  ### ignore bins with 0 for both tumour and normal
  ### targetedSequence = RangedData (IRanges object)
  ### containing list of targeted regions to consider;
  ### 3 columns: chr, start, end
  if (!is.null(targetedSequence)) {
    message("Analyzing targeted regions...")
    targetIR <- RangedData(ranges = IRanges(start = targetedSequence[, 2], 
                                            end = targetedSequence[, 3]), space = targetedSequence[, 1])
    
    #keepInd <- unlist(as.list(findOverlaps(tumour_reads, targetIR, select = "first")))
    #keepInd <- !is.na(keepInd)
    hits <- findOverlaps(query = tumour_reads, subject = targetIR)
    keepInd <- unique(queryHits(hits))    
    
    # ind <- tumour_reads$value>10 &
    # normal_reads$value>10 tumThres <-
    # quantile(tumour_reads$value[ind],1/4) normThres
    # <- quantile(normal_reads$value[ind],1/4) keepInd
    # <- which(ind & !is.na(tumour_reads$value) &
    # !is.na(normal_reads$value) &
    # tumour_reads$value>=tumThres &
    # normal_reads$value>=normThres)
    tumour_reads <- tumour_reads[keepInd, ]
    normal_reads <- normal_reads[keepInd, ]
    gc <- gc[keepInd, ]
    map <- map[keepInd, ]
    samplesize <- min(ceiling(nrow(tumour_reads) * 
                                0.1), samplesize)
  }
  
  ### add GC and Map data to IRanges objects ###
  tumour_reads$gc <- gc$value
  tumour_reads$map <- map$value
  colnames(tumour_reads) <- c("reads", "gc", "map")
  normal_reads$gc <- gc$value
  normal_reads$map <- map$value
  colnames(normal_reads) <- c("reads", "gc", "map")
  
  ### CORRECT TUMOUR DATA FOR GC CONTENT AND
  ### MAPPABILITY BIASES ###
  message("Correcting Tumour")
  tumour_copy <- correctReadcount(tumour_reads,mappability = 0.5, samplesize = samplesize)
  
  ### CORRECT NORMAL DATA FOR GC CONTENT AND
  ### MAPPABILITY BIASES ###
  message("Correcting Normal")
  normal_copy <- correctReadcount(normal_reads, mappability = 0.5, samplesize = samplesize)
  
  ### COMPUTE LOG RATIO ###
  message("Normalizing Tumour by Normal")
  tumour_copy$copy <- tumour_copy$copy - normal_copy$copy
  rm(normal_copy)
  
  ### PUTTING TOGETHER THE COLUMNS IN THE OUTPUT ###
  temp <- cbind(chr = as.character(space(tumour_copy)), 
                start = start(tumour_copy), end = end(tumour_copy), 
                logR = tumour_copy$copy)
  temp <- as.data.frame(temp, stringsAsFactors = FALSE)
  mode(temp$start) <- "numeric"
  mode(temp$end) <- "numeric"
  mode(temp$logR) <- "numeric"
  return(temp)
}
#####################################################################
#####################################################################
#####################################################################

#####################################################################
############### Change Plotting Parameters###########################
#####################################################################

getGenomeWidePositions <- function(chrs,posns){  
  #create genome coordinate scaffold
  positions <- as.numeric(posns)
  chrsNum <- unique(chrs)
  chrBkpt <- rep(0,length(chrsNum)+1)
  for (i in 2:length(chrsNum)){
    chrInd <- which(chrs==chrsNum[i])
    prevChrPos <- positions[chrInd[1]-1]      
    chrBkpt[i] = prevChrPos
    positions[chrInd] = positions[chrInd] + prevChrPos
  }
  chrBkpt[i+1] <- positions[length(positions)]
  return(list(posns=positions,chrBkpt=chrBkpt))
}
plotChrLines_JYedit2 <- function(chrs, chrBkpt, yrange) {
  # plot vertical chromosome lines
  for (j in 1:length(chrBkpt)) {
    lines(rep(chrBkpt[j], 2), yrange, type = "l", 
          lty = 2, col = "black", lwd = 0.75)
  }
  numLines <- length(chrBkpt)
  mid <- (chrBkpt[1:(numLines - 1)] + chrBkpt[2:numLines])/2
  chrs[chrs == "X"] <- 23
  chrs[chrs == "Y"] <- 24
  chrsToShow <- sort(unique(as.numeric(chrs)))
  chrsToShow[chrsToShow == 23] <- "X"
  chrsToShow[chrsToShow == 24] <- "Y"
  axis(side = 1, at = mid, labels = c(gsub("contig_|_alt", "", unique(chrs))), 
       cex.axis = 0.7, tick = FALSE)
}
plotChrLines_JYedit <- function(chrs, chrBkpt, yrange) {
  # plot vertical chromosome lines
  for (j in 1:length(chrBkpt)) {
    lines(rep(chrBkpt[j], 2), yrange, type = "l", 
          lty = 2, col = "black", lwd = 0.75)
  }
  numLines <- length(chrBkpt)
  mid <- (chrBkpt[1:(numLines - 1)] + chrBkpt[2:numLines])/2
  chrs[chrs == "X"] <- 23
  chrs[chrs == "Y"] <- 24
  chrsToShow <- sort(unique(as.numeric(chrs)))
  chrsToShow[chrsToShow == 23] <- "X"
  chrsToShow[chrsToShow == 24] <- "Y"
  axis(side = 1, at = mid, labels = c(gsub("contig_|_alt", "", unique(chrs))), 
       cex.axis = 0.7, tick = FALSE)
}

plotCNlogRByChr_JYedit <- function(dataIn, chr = NULL, segs = NULL, geneAnnot = NULL, 
                                   ploidy = NULL, normal = NULL, spacing = 4, alphaVal = 1, xlim = NULL, ...) {
  # color coding
  alphaVal <- ceiling(alphaVal * 255)
  class(alphaVal) = "hexmode"
  cnCol <- c("#00FF00", "#006400", "#0000FF", "#880000", 
             "#BB0000", "#CC0000", "#DD0000", "#EE0000", "#FF0000")
  cnCol <- paste(cnCol, alphaVal, sep = "")
  # cnCol <-
  # col2rgb(c('green','darkgreen','blue','darkred','red','brightred'))
  names(cnCol) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8")
  
  ## adjust logR values for ploidy ##
  if (!is.null(ploidy)) {
    if (is.null(normal)){
      stop("plotCNlogRByChr: Please provide \"normal\" contamination estimate.")
    }
    dataIn[, "LogRatio"] <- as.numeric(unlist(dataIn[, "LogRatio"])) + log2(((1-normal)*ploidy+normal*2)/2)
    
    if (!is.null(segs)){
      segs[, "Median_logR"] <- segs[, "Median_logR"] + log2(((1-normal)*ploidy+normal*2) / 2)
    }
  }
  
  if (!is.null(chr)) {
    for (i in chr) {
      dataByChr <- dataIn[dataIn[, "Chr"] == 
                            i, ]
      dataByChr <- dataByChr[dataByChr[, "TITANcall"] != "OUT", ]
      # plot the data if (outfile!=''){
      # pdf(outfile,width=10,height=6) }
      par(mar = c(spacing, 8, 2, 2))
      # par(xpd=NA)
      if (missing(xlim)) {
        xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), 
                                          "Position"]))
      }
      coord <- as.numeric(dataByChr[, "Position"])
      plot(coord, as.numeric(dataByChr[, "LogRatio"]), 
           col = cnCol[as.character(dataByChr[, 
                                              "CopyNumber"])], pch = 16, xaxt = "n", 
           las = 1, ylab = "Copy Number (log ratio)", 
           xlim = xlim, ...)
      lines(xlim, rep(0, 2), type = "l", col = "grey", lwd = 0.75)
      if (!is.null(segs)){
        segsByChr <- segs[segs[,"Chromosome"]==as.character(i), , drop=FALSE]
        tmp <- apply(segsByChr, 1, function(x){
          lines(x[c("Start_Position.bp.","End_Position.bp.")], 
                rep(x["Median_logR"], 2), col = "green", lwd = 3)
        })
      }
      
      
      if (!is.null(geneAnnot)) {
        plotGeneAnnotation(geneAnnot, i)
      }
    }
  } else {
    # plot for all chromosomes
    coord <- getGenomeWidePositions(as.matrix(dataIn[, "Chr"]), 
                                    as.matrix(dataIn[, "Position"]))
    plot(coord$posns, as.numeric(unlist(dataIn[, "LogRatio"])), 
         col = cnCol[as.character(unlist(dataIn[, "CopyNumber"]))], 
         pch = 16, xaxt = "n", las = 1, bty = "n", 
         ylab = "Copy Number (log ratio)", ...)
    lines(as.numeric(c(1, coord$posns[length(coord$posns)])), 
          rep(0, 2), type = "l", col = "grey", lwd = 2)
    plotChrLines_JYedit2(as.matrix(dataIn[, "Chr"]), coord$chrBkpt, par("yaxp")[1:2])
    #plot segments
    if (!is.null(segs)){
      coordEnd <- getGenomeWidePositions(segs[, "Chromosome"], segs[, "End_Position.bp."])
      coordStart <- coordEnd$posns - (segs[, "End_Position.bp."] - segs[, "Start_Position.bp."] + 1)
      xlim <- as.numeric(c(1, coordEnd$posns[length(coordEnd$posns)]))
      col <- cnCol[as.character(segs[, "Copy_Number"])]
      value <- as.numeric(segs[, "Median_logR"])
      mat <- as.data.frame(cbind(coordStart, coordEnd$posns, value, col))
      rownames(mat) <- 1:nrow(mat)
      tmp <- apply(mat, 1, function(x){
        lines(x[1:2], rep(x[3], 2), col = x[4], lwd = 3)
      })
    }
    
  }
  
}

plotAllelicRatio_JYedit <- function(dataIn, chr = NULL, geneAnnot = NULL, 
                                    spacing = 4,  xlim = NULL, ...) {
  # color coding alphaVal <- ceiling(alphaVal * 255);
  # class(alphaVal) = 'hexmode'
  lohCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000", 
              "#006400", "#BEBEBE", "#FF0000", "#BEBEBE", 
              "#FF0000")
  # lohCol <- paste(lohCol,alphaVal,sep='') lohCol <-
  # col2rgb(c('green','darkgreen','blue','darkgreen','grey','red'))
  names(lohCol) <- c("HOMD", "DLOH", "NLOH", "GAIN", 
                     "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")
  
  
  if (!is.null(chr)) {
    for (i in chr) {
      dataByChr <- dataIn[dataIn[, "Chr"] == 
                            i, ]
      dataByChr <- dataByChr[dataByChr[, "TITANcall"] != 
                               "OUT", ]
      # plot the data if (outfile!=''){
      # pdf(outfile,width=10,height=6) }
      par(mar = c(spacing, 8, 2, 2))
      # par(xpd=NA)
      if (missing(xlim)) {
        xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), 
                                          "Position"]))
      }
      plot(dataByChr[, "Position"], dataByChr[, 
                                              "AllelicRatio"], col = lohCol[dataByChr[, 
                                                                                      "TITANcall"]], pch = 16, xaxt = "n", 
           las = 1, ylab = "Allelic Ratio", xlim = xlim, 
           ...)
      lines(as.numeric(c(1, dataByChr[nrow(dataByChr), 
                                      "Position"])), rep(0.5, 2), type = "l", 
            col = "grey", lwd = 3)
      
      if (!is.null(geneAnnot)) {
        plotGeneAnnotation(geneAnnot, i)
      }
    }
  } else {
    # plot for all chromosomes
    coord <- getGenomeWidePositions(as.matrix(dataIn[, "Chr"]), 
                                    as.matrix(dataIn[, "Position"]))
    plot(coord$posns, as.numeric(unlist(dataIn[, "AllelicRatio"])), 
         col = lohCol[unlist(dataIn[, "TITANcall"])], pch = 16, 
         xaxt = "n", bty = "n", las = 1, ylab = "Allelic Ratio", xlab="Scaffold",  
         ...)
    lines(as.numeric(c(1, coord$posns[length(coord$posns)])), 
          rep(0.5, 2), type = "l", col = "grey", 
          lwd = 3)
    plotChrLines_JYedit(unique(as.matrix(dataIn[, "Chr"])), coord$chrBkpt, 
                        c(-0.1, 1.1))
    
  }
}

plotClonalFrequency_JYedit <- function(dataIn, chr = NULL, 
                                       normal = NULL, geneAnnot = NULL, spacing = 4, xlim = NULL, ...) {
  # color coding
  lohCol <- c("#00FF00", "#006400", "#0000FF", "#8B0000", 
              "#006400", "#BEBEBE", "#FF0000", "#FF0000", 
              "#FF0000")
  names(lohCol) <- c("HOMD", "DLOH", "NLOH", "GAIN", 
                     "ALOH", "HET", "ASCNA", "BCNA", "UBCNA")
  
  # get unique set of cluster and estimates table:
  # 1st column is cluster number, 2nd column is
  # clonal freq
  clusters <- unique(dataIn[, c("ClonalCluster", 
                                "CellularPrevalence")])
  clusters <- clusters[!is.na(unlist(clusters[, 1])), , drop = FALSE]  #exclude NA
  if (!is.null(normal)) {
    clusters[, 2] <- (as.numeric(clusters[, 2])) * 
      (1 - as.numeric(normal))
  }
  
  dataToUse <- dataIn[unlist(dataIn[, "TITANcall"]) != "OUT", ]
  dataToUse[unlist(dataToUse[, "CellularPrevalence"]) == 
              "NA" | is.na(unlist(dataToUse[, "CellularPrevalence"])), 
            c("ClonalCluster", "CellularPrevalence")] <- c(NA, NA)
  # extract clonal info
  clonalFreq <- cbind(as.matrix(dataToUse[, "ClonalCluster"]), 
                      as.matrix(dataToUse[, "CellularPrevalence"]))
  # mode(clonalFreq) <- 'numeric' clonalFreq[,2] <- 1
  # - clonalFreq[,2]
  if (!is.null(normal)) {
    clonalFreq[, 2] <- clonalFreq[, 2] * (1 - normal)
  }
  clonalFreq[is.na(clonalFreq[, 2]) | clonalFreq[, 
                                                 2] == "0" | clonalFreq[, 2] == "NA", 2] <- 0
  
  # plot per chromosome
  if (!is.null(chr)) {
    for (i in chr) {
      ind <- dataToUse[, "Chr"] == as.character(i)
      dataByChr <- dataToUse[ind, ]
      clonalFreq <- clonalFreq[ind, ]
      # plot the data
      par(mar = c(spacing, 8, 2, 2), xpd = NA)
      # par(xpd=NA)
      
      # PLOT CLONAL FREQUENCIES
      if (missing(xlim)) {
        xlim <- as.numeric(c(1, dataByChr[nrow(dataByChr), 
                                          "Position"]))
      }
      plot(dataByChr[, "Position"], clonalFreq[, 
                                               2], type = "h", col = lohCol[dataByChr[, 
                                                                                      "TITANcall"]], las = 1, xaxt = "n", 
           ylab = "Cellular Prevalence", xlim = xlim, 
           ...)
      
      # plot cluster lines and labels
      if (nrow(clusters) > 0){
        for (j in 1:length(clusters[, 1])) {
          chrLen <- as.numeric(dataByChr[dim(dataByChr)[1], 
                                         "Position"])
          lines(c(1 - chrLen * 0.02, chrLen * 
                    1.02), rep(clusters[j, 2], 2), type = "l", 
                col = "grey", lwd = 3)
          mtext(side = 4, at = clusters[j, 2], 
                text = paste("Z", clusters[j, 1], 
                             "", sep = ""), cex = 1, padj = 0.5, 
                adj = 1, las = 2, outer = FALSE)
          mtext(side = 2, at = clusters[j, 2], 
                text = paste("Z", clusters[j, 1], 
                             "", sep = ""), cex = 1, padj = 0.5, 
                adj = 0, las = 2, outer = FALSE)
        }
      }
      
      if (!is.null(normal)) {
        chrLen <- as.numeric(dataByChr[nrow(dataByChr), 
                                       "Position"])
        lines(c(1 - chrLen * 0.02, chrLen * 
                  1.02), rep((1 - normal), 2), type = "l", 
              col = "#000000", lwd = 3)
        #mtext(side = 4, at = (1 - normal), 
        #text = paste("-T-", sep = ""), padj = 0.5, 
        #adj = 1, cex = 1, las = 2, outer = FALSE)
        #mtext(side = 2, at = (1 - normal), 
        #text = paste("-T-", sep = ""), padj = 0.5, 
        #adj = 0, cex = 1, las = 2, outer = FALSE)
      }
      
      if (!is.null(geneAnnot)) {
        plotGeneAnnotation(geneAnnot, i)
      }
    }
  } else {
    # plot genome-wide
    coord <- getGenomeWidePositions(as.matrix(dataIn[, "Chr"]), 
                                    as.matrix(dataIn[, "Position"]))
    plot(coord$posns, clonalFreq[, 2], type = "h", 
         col = lohCol[as.matrix(dataIn[, "TITANcall"])], pch = 16, 
         xaxt = "n", las = 1, bty = "n", ylab = "Cellular Prevalence",  xlab="Scaffold",
         ...)
    plotChrLines_JYedit(unique(as.matrix(dataIn[, "Chr"])), coord$chrBkpt, 
                        c(-0.1, 1.1))
    
    # plot cluster lines and labels
    for (j in 1:length(clusters[, 1])) {
      chrLen <- as.numeric(coord$posns[length(coord$posns)])
      lines(c(1 - chrLen * 0.02, chrLen * 1.02), 
            rep(clusters[j, 2], 2), type = "l", 
            col = "grey", lwd = 3)
      mtext(side = 4, at = clusters[j, 2], text = paste("Z", 
                                                        clusters[j, 1], "", sep = ""), cex = 1, 
            padj = 0.5, adj = 1, las = 2, outer = FALSE)
      mtext(side = 2, at = clusters[j, 2], text = paste("Z", 
                                                        clusters[j, 1], "", sep = ""), cex = 1, 
            padj = 0.5, adj = 0, las = 2, outer = FALSE)
    }
    if (!is.null(normal)) {
      chrLen <- as.numeric(coord$posns[length(coord$posns)])
      lines(c(1 - chrLen * 0.02, chrLen * 1.02), 
            rep((1 - normal), 2), type = "l", col = "#000000", 
            lwd = 3)
    }
    
  }
  
}


#####################################################################
#####################################################################
#####################################################################


id <- paste0(Isolate)

#infile comes from extractAlleleReadCounts() which outputs a text file
#extdata is from Users/jdyuzon/Library/R/3.3/library/TitanCNA/extdata/
infilename<-paste0(Isolate,"_",Normal,"_tumAlleleCounts.tsv")
infile <- infilename
data <- loadAlleleCounts_JYedit(infile, symmetric=TRUE)
cnData <- correctReadDepth_JYedit(tum_reads.wig,norm_reads.wig,gc.wig,map.wig)
cnData[[3]]<-as.numeric(cnData[[3]])
cnData2chrNo<-as.numeric(gsub("contig_|_alt", "", cnData$chr))
cnData$chrNo<-cnData2chrNo
cnData2chrNo<-as.data.frame(cnData)
cnData2 <- cnData2chrNo[order(cnData2chrNo$chrNo),]
cnData2$chrNo<-NULL

#Assign copy number to each position
#order chromosomes and positions of data and cnData
occurances<-as.data.frame(table(data$chr))
occurances2<-subset(occurances, occurances$Freq > 1)
df<-as.data.frame(data)
data1<-subset(df, chr %in% as.character(occurances2$Var1))
data2<-list()
data2$chr<-as.character(data1$chr)
data2$posn <-data1$posn
data2$ref <-data1$ref
data2$refOriginal <-data1$refOriginal
data2$nonRef <-data1$nonRef
data2$tumDepth <-data1$tumDepth
data2$logR <-data1$logR
#order chromosomes in sequential order
data2chr<-as.character(data2$chr)
data2chrNo<-as.numeric(gsub("contig_|_alt", "", data2chr))
data2$chrNo<-data2chrNo
data2chrNo<-as.data.frame(data2)
data2 <- data2chrNo[order(data2chrNo$chrNo),]
data2$chrNo<-NULL
tmpdata2<-NULL
tmpdata2$chr<-data2$chr
tmpdata2$posn<-data2$posn
tmpdata2$ref<-data2$ref
tmpdata2$refOriginal<-data2$refOriginal
tmpdata2$nonRef<-data2$nonRef
tmpdata2$tumDepth<-data2$tumDepth
tmpdata2$logR<-data2$logR
data2<-tmpdata2
data2$chr<-as.character(data2$chr)



logR <- getPositionOverlap(data2$chr,data2$posn,cnData2)
logR<-as.numeric(logR)
data2$logR <- log(2^logR)
rm(logR,cnData)

#occurances<-as.data.frame(table(data$chr))
#occurances2<-subset(occurances, occurances$Freq > 1)
#df<-as.data.frame(data)
#data1<-subset(df, chr %in% as.character(occurances2$Var1))
#data2<-list()
#data2$chr<-as.character(data1$chr)
#data2$posn <-data1$posn
#data2$ref <-data1$ref
#data2$refOriginal <-data1$refOriginal
#data2$nonRef <-data1$nonRef
#data2$tumDepth <-data1$tumDepth
#data2$logR <-data1$logR

#order chromosomes in sequential order
#data2chr<-as.character(data2$chr)
#data2chrNo<-as.numeric(gsub("contig_", "", data2chr))
#data2$chrNo<-data2chrNo 
#data2chrNo<-as.data.frame(data2)
#data2 <- data2chrNo[order(data2chrNo$chrNo),] 
#data2$chrNo<-NULL
#tmpdata2<-NULL
#tmpdata2$chr<-data2$chr
#tmpdata2$posn<-data2$posn
#tmpdata2$ref<-data2$ref
#tmpdata2$refOriginal<-data2$refOriginal
#tmpdata2$nonRef<-data2$nonRef
#tmpdata2$tumDepth<-data2$tumDepth
#tmpdata2$logR<-data2$logR
#data2<-tmpdata2

data2 <- transform(data2, chr = as.character(chr), ref = as.numeric(ref), refOriginal=as.numeric(refOriginal), nonRef=as.numeric(nonRef), tumDepth=as.numeric(tumDepth), logR=as.numeric(logR))
contigs<-NULL
contigs$chr<-c(as.character(paste0("contig_", 1:302)),as.character(paste0("contig_",1:302,"_alt")))
contigs$num<-as.numeric(gsub("contig_|_alt", "", contigs$chr))
contigs<-as.data.frame(contigs)
contigs2<-contigs[order(contigs$num),]

#Filter out low and high depth positions and positions with NA’s.
#change 1:22,"X" to Pram scaffolds
#Users can also choose to remove positions that are highly mappable.
mScore <- as.data.frame(wigToRangedData(map.wig))
mScore <- getPositionOverlap(data2$chr,data2$posn,mScore[,-4])



##########################################################################################################
############################################## Haplotypes 1 ##############################################
##########################################################################################################
foo <- filterData(data2, c(as.character(paste0("contig_", 1:302))),minDepth=5,maxDepth=200,map=mScore,mapThres=0.8)

#Parameter estimation
convergeParams <- runEMclonalCN(foo, params,maxiter = 20, maxiterUpdate = 1500,
                                txnExpLen = 1e15, txnZstrength = 1e5,  useOutlierState = FALSE,
                                normalEstimateMethod = "map", estimateS = TRUE, estimatePloidy = TRUE)

#Optimal genotype/clonal cluster inference
optimalPath <- viterbiClonalCN(foo,convergeParams)

#output formatted results
TitanCNAOutfile<-paste0("CNV_", Isolate,"_",Normal,"_clust" ,numClusters, "_titan_H1.txt", sep="")
results<-outputTitanResults(foo,convergeParams,optimalPath,filename=TitanCNAOutfile,posteriorProbs=FALSE,subcloneProfiles=TRUE)

#outparam <- paste("CNV_", Isolate,"_",Normal, "_cluster0",numClusters,"_params.txt",sep="")
#outputModelParameters(convergeParams,results$results,outparam)
norm <- convergeParams$n[length(convergeParams$n)]
ploidy <- convergeParams$phi[length(convergeParams$phi)]

outplot<-paste("CNVsc_", Isolate,"_",Normal, "_clust_0", numClusters, "_H1.png", sep = "")
png(outplot,width=1200,height=1000,res=100)
par(mfrow=c(3,1))
plotCNlogRByChr_JYedit(results$results, normal = norm, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-4,6),cex=0.5,main=paste0("\"", Isolate, "\""))
plotAllelicRatio_JYedit(results$results,  geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5)
plotClonalFrequency_JYedit(results$results,  normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5,main=paste0("\"", Isolate, "\""))
dev.off()
##########################################################################################################
############################################## Haplotypes 2 ##############################################
##########################################################################################################
foo <- filterData(data2, c(as.character(paste0("contig_", 1:302,"_alt"))),minDepth=5,maxDepth=200,map=mScore,mapThres=0.8)

#Parameter estimation
convergeParams <- runEMclonalCN(foo, params,maxiter = 20, maxiterUpdate = 1500,
                                txnExpLen = 1e15, txnZstrength = 1e5,  useOutlierState = FALSE,
                                normalEstimateMethod = "map", estimateS = TRUE, estimatePloidy = TRUE)

#Optimal genotype/clonal cluster inference
optimalPath <- viterbiClonalCN(foo,convergeParams)

#output formatted results
TitanCNAOutfile<-paste0("CNV_", Isolate,"_",Normal,"_clust" ,numClusters, "_titan_H2.txt", sep="")
results<-outputTitanResults(foo,convergeParams,optimalPath,filename=TitanCNAOutfile,posteriorProbs=FALSE,subcloneProfiles=TRUE)

#outparam <- paste("CNV_", Isolate,"_",Normal, "_cluster0",numClusters,"_params.txt",sep="")
#outputModelParameters(convergeParams,results$results,outparam)
norm <- convergeParams$n[length(convergeParams$n)]
ploidy <- convergeParams$phi[length(convergeParams$phi)]

outplot<-paste("CNVsc_", Isolate,"_",Normal, "_clust_0", numClusters, "_H2.png", sep = "")
png(outplot,width=1200,height=1000,res=100)
par(mfrow=c(3,1))
plotCNlogRByChr_JYedit(results$results, normal = norm, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-4,6),cex=0.5,main=paste0("\"", Isolate, "\""))
plotAllelicRatio_JYedit(results$results,  geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5)
plotClonalFrequency_JYedit(results$results,  normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5,main=paste0("\"", Isolate, "\""))
dev.off()

##########################################################################################################
###################################### Whole Genome:Both Haplotypes ######################################
##########################################################################################################
foo <- filterData(data2, c(as.character(paste0("contig_", 1:302)),
	as.character(paste0("contig_", 1:302,"_alt"))),minDepth=5,maxDepth=200,map=mScore,mapThres=0.8)

#Parameter estimation
convergeParams <- runEMclonalCN(foo, params,maxiter = 20, maxiterUpdate = 1500, 
				txnExpLen = 1e15, txnZstrength = 1e5,  useOutlierState = FALSE,  
				normalEstimateMethod = "map", estimateS = TRUE, estimatePloidy = TRUE)

#Optimal genotype/clonal cluster inference
optimalPath <- viterbiClonalCN(foo,convergeParams)

#output formatted results
TitanCNAOutfile<-paste0("CNV_", Isolate,"_",Normal,"_clust" ,numClusters, "_titan.txt", sep="")
results<-outputTitanResults(foo,convergeParams,optimalPath,filename=TitanCNAOutfile,posteriorProbs=FALSE,subcloneProfiles=TRUE)

outparam <- paste("CNV_", Isolate,"_",Normal, "_cluster0",numClusters,"_params.txt",sep="")
outputModelParameters(convergeParams,results$results,outparam)

#Plot results:
#CNV: The Y-axis is based on log ratios. Log ratios are computed ratios between normalized tumour and normal read depths. Data points close to 0 represent diploid, above 0 are copy gains, below 0 are deletions. Bright Green - HOMD Green - DLOH Blue - HET, NLOH Dark Red - GAIN Red - ASCNA, UBCNA, BCNA
#Allele: The Y-axis is based on allelic ratios. Allelic ratios are computed as RefCount/Depth. Data points close to 1 represent homozygous reference base, close to 0 represent homozygous non-reference base, and close to 0.5 represent heterozygous. Normal contamination influences the divergence away from 0.5 for LOH events. Grey - HET, BCNA Bright Green - HOMD Green - DLOH, ALOH Blue - NLOH Dark Red - GAIN Red - ASCNA, UBCNA
norm <- convergeParams$n[length(convergeParams$n)]
ploidy <- convergeParams$phi[length(convergeParams$phi)]
#library(SNPchip)  ## use this library to plot chromosome idiogram (optional)


#x-axis are positions along the genome if the genome were one scaffold. Therefore, scaffold_2 positions are added to scaffold_1's last position
#outplot<-paste("CNV_", Isolate,"_",Normal, "_clust_0", numClusters, ".png", sep = "")
#png(outplot,width=1200,height=1000,res=100)
#par(mfrow=c(3,1))
#plotCNlogRByChr(results, normal = norm, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-4,6),cex=0.5,main=paste0("\"", Isolate, "\""))
#plotAllelicRatio(results,  geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5)
#plotClonalFrequency(results,  normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5,main=paste0("\"", Isolate, "\""))
#dev.off()
#x-axis are scaffold numbers
outplot<-paste("CNVsc_", Isolate,"_",Normal, "_clust_0", numClusters, ".png", sep = "")
png(outplot,width=1200,height=1000,res=100)
par(mfrow=c(3,1))
plotCNlogRByChr_JYedit(results$results, normal = norm, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-4,6),cex=0.5,main=paste0("\"", Isolate, "\""))
plotAllelicRatio_JYedit(results$results,  geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5)
plotClonalFrequency_JYedit(results$results,  normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5,main=paste0("\"", Isolate, "\""))
dev.off()

#if (as.numeric(numClusters) <= 1){ 
## NEW IN V1.2.0 ##
## users can choose to plot the subclone copy number profiles for <= 2 clusters
#	plotSubcloneProfiles(results, chr, cex = 2, spacing=6, main=paste0("\"", Isolate, "\""))
#}
#pI <- plotIdiogram(chr,build="hg19",unit="bp",label.y=-4.25,new=FALSE,ylim=c(-2,-1))
#dev.off()

id <- paste0(" -id ",Isolate)
infile<- paste0(" -infile ",TitanCNAOutfile)
outfile<- paste0(" -outfile ", "CNV_", Isolate,"_",Normal, "tmp.bed")
outIGV<- paste0(" -outIGV ", "CNV_tmpIGV.",TitanCNAOutfile)
outfilemask<-paste0("CNV_", Isolate,"_",Normal, ".mask.bed")
createTITANsegmentfiles <- paste("perl", "createTITANsegmentfiles_JYedit.pl", id, infile, outfile, outIGV)
system(createTITANsegmentfiles)

#Get Mask File
#grep -Ev "HET|BCNA|GAIN" outfile | awk '$5!=1' | awk '{print $2 "\t" $3 "\t" $4 }' > outfilebed
outfile.0.txt = read.table(paste0("CNV_", Isolate,"_",Normal,"tmp.bed"))
outfile.1.txt = outfile.0.txt[-1,]

scaffoldNo<-as.character(outfile.1.txt$V2)
scaffoldNo1<-as.numeric(gsub("contig_", "", scaffoldNo))
outfile.1.txt$y<-scaffoldNo1 
outfile.1.txt <- outfile.1.txt[order(outfile.1.txt$y),] 

outfile.2.txt <- outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9!="HET" & outfile.1.txt$V9!="BCNA" & outfile.1.txt$V9!="GAIN" & outfile.1.txt$V9!="ASCNA" & outfile.1.txt$V9!="UBCNA", ]
outfile.3.txt <- subset(outfile.2.txt, select=c(V2, V3, V4 ))

#to remove first base pair of SV for the mask need, subtract 1 from the start position of the SV in the mask file
outfile.3a.txt <-outfile.3.txt
outfile.3a.txt$V3<-as.numeric(as.character(outfile.3.txt$V3))-1

write.table(outfile.3a.txt, outfilemask , quote = FALSE, sep="\t", row.names = FALSE,
            col.names = FALSE)

###
#make file for CNA 
outfileCNVtracks<-paste0("CNV_", Isolate,"_",Normal,"_CNVtracks.bed")
outfile.1.txt$x<-paste0(Isolate, outfile.1.txt$V9,outfile.1.txt$V3,outfile.1.txt$V4 )
outfile.4.txt <-outfile.1.txt[as.numeric(outfile.1.txt$V5) > 2 & outfile.1.txt$V9!="HET" , ]
outfile.5.txt<-subset(outfile.4.txt, select=c(V2,V3,V4,x))


write.table(outfile.5.txt, outfileCNVtracks , quote = FALSE,  sep="\t", row.names = FALSE,
            col.names = FALSE)

#bed file for each CNA type to use in Phylogeny
LOH<-subset(outfile.1.txt[outfile.1.txt$V9=="NLOH"|outfile.1.txt$V9=="ALOH"|outfile.1.txt$V9=="DLOH", ],select=c(V2,V3,V4,x))
TRI<-subset(outfile.1.txt[outfile.1.txt$V9=="GAIN"|outfile.1.txt$V9=="ASCNA"|outfile.1.txt$V9=="BCNA"|outfile.1.txt$V9=="UBCNA", ],select=c(V2,V3,V4,x))
#fltLOH<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & as.numeric(outfile.1.txt$V6)>0.99 & outfile.1.txt$V9=="NLOH"|outfile.1.txt$V9=="ALOH"|outfile.1.txt$V9=="DLOH", ],select=c(V2,V3,V4,x))
#fltTRI<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) >=1000 & as.numeric(outfile.1.txt$V6)>=0.65 & as.numeric(outfile.1.txt$V6)<=0.69 & as.numeric(outfile.1.txt$V10)==3, ],select=c(V2,V3,V4,x))
#ALOH<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="ALOH" , ],select=c(V2,V3,V4,x))
#DLOH<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="DLOH" , ],select=c(V2,V3,V4,x))
#NLOH<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="NLOH" , ],select=c(V2,V3,V4,x))
#ASCNA<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="ASCNA" , ],select=c(V2,V3,V4,x))
#BCNA<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="BCNA" , ],select=c(V2,V3,V4,x))
#UBCNA<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="UBCNA" , ],select=c(V2,V3,V4,x))
#GAIN<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="GAIN" , ],select=c(V2,V3,V4,x))
#HOMD<-subset(outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9=="HOMD" , ],select=c(V2,V3,V4,x))

write.table(LOH, paste0(Isolate,"_",Normal,"_LOH.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
write.table(TRI, paste0(Isolate,"_",Normal,"_TRI.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)

#write.table(fltLOH, paste0("CNV_",Isolate,"_",Normal, "_fltLOH.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(fltTRI, paste0("CNV_",Isolate,"_",Normal, "_fltTRI.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(ALOH, paste0("CNV_",Isolate,"_",Normal, "_ALOH.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(DLOH, paste0("CNV_",Isolate,"_",Normal, "_DLOH.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(NLOH, paste0("CNV_",Isolate,"_",Normal, "_NLOH.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(ASCNA, paste0("CNV_",Isolate,"_",Normal, "_ASCNA.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(BCNA, paste0("CNV_",Isolate,"_",Normal, "_BCNA.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(UBCNA, paste0("CNV_",Isolate,"_",Normal, "_UBCNA.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(GAIN, paste0("CNV_",Isolate,"_",Normal, "_GAIN.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)
#write.table(HOMD,paste0("CNV_",Isolate,"_",Normal, "_HOMD.bed") , quote = FALSE,  sep="\t", row.names = FALSE, col.names = FALSE)

###PHASING
txt<-read.table(TitanCNAOutfile)
tsv<-read.table(infilename)
phase.tsv<-merge(txt,tsv, by = c('V1', 'V2'), sort=FALSE) #merge so have allele counts, alleles, copy number
phase.tsv$ampHaplo<-NA
phase.tsv$altHaplo<-NA

for(i in 1:nrow(phase.tsv)){
    row<-phase.tsv[i,]
    if(row$V8=='2' & row$V10 !="HET" | row$V8!='2'){
    if(as.numeric(levels(row$V3.x))[row$V3.x]>as.numeric(levels(row$V4.x))[row$V4.x]){
      phase.tsv[i,20]<-as.character(row$V3.y)
      phase.tsv[i,21]<-as.character(row$V5.y)
      }
    else if(as.numeric(levels(row$V3.x))[row$V3.x]<as.numeric(levels(row$V4.x))[row$V4.x]){
      phase.tsv[i,20]<-as.character(row$V5.y)
      phase.tsv[i,21]<-as.character(row$V3.y)
      }
    else {
      phase.tsv[i,20]<-"-"
      phase.tsv[i,21]<-"-"
      }
  }
  else {
    phase.tsv[i,20]<-"-"
    phase.tsv[i,21]<-"-"
    }
}
phaseq<-data.frame(phase.tsv$V1, phase.tsv$V2, phase.tsv$ampHaplo, phase.tsv$altHaplo)
write.table(phaseq, paste0("CNV_", Isolate,"_",Normal,"_phaseq") , quote = FALSE,  sep="\t", row.names = FALSE,
            col.names = FALSE)
#have to merge with vcf that has all positions in consecutive order and allelecounts; phase those out as well
#OR use a Trisomy+LOH bed file to subset a vcf file (that has all positions in consecutive order and allelecounts), and output alleles with highest read coverage and lowest read coverage to phase

###Getting Allele Counts to later get population allele frequency
phase.tsv$RefAlleleCount<-NA
phase.tsv$AltAlleleCount<-NA
for(i in 1:nrow(phaseq)){
  row<-phase.tsv[i,]
  if(row$V10 =="HET"){
    phase.tsv[i,22]<-1
    phase.tsv[i,23]<-1
  }
  else if(row$V10 !="HET" & as.numeric(levels(row$V3.x))[row$V3.x]!=0 & as.numeric(levels(row$V4.x))[row$V4.x]!=0){
    if(as.numeric(levels(row$V3.x))[row$V3.x] > as.numeric(levels(row$V4.x))[row$V4.x]){
      phase.tsv[i,22]<-2
      phase.tsv[i,23]<-1
    }
    else if(as.numeric(levels(row$V3.x))[row$V3.x]<as.numeric(levels(row$V4.x))[row$V4.x]){
      phase.tsv[i,22]<-1
      phase.tsv[i,23]<-2
    }
    else{
      phase.tsv[i,22]<-1
      phase.tsv[i,23]<-1
    }
  }
  else if(row$V10 !="HET" & (as.numeric(levels(row$V3.x))[row$V3.x]==0 | as.numeric(levels(row$V4.x))[row$V4.x]==0)){
    if(as.numeric(levels(row$V3.x))[row$V3.x]>as.numeric(levels(row$V4.x))[row$V4.x]){
      phase.tsv[i,22]<-2
      phase.tsv[i,23]<-0
    }
    else if(as.numeric(levels(row$V3.x))[row$V3.x]<as.numeric(levels(row$V4.x))[row$V4.x]){
      phase.tsv[i,22]<-0
      phase.tsv[i,23]<-2
    }
    else{
      phase.tsv[i,22]<-1
      phase.tsv[i,23]<-1 
    }
  }
}
Isolateallelecount<-data.frame(phase.tsv$V1, phase.tsv$V2, phase.tsv$V3.y, phase.tsv$V5.y, phase.tsv$RefAlleleCount, phase.tsv$AltAlleleCount)
write.table(Isolateallelecount, paste0("CNV_", Isolate,"_",Normal,"_allelecounts_for_freq.bed") , quote = FALSE,  sep="\t", row.names = FALSE,
            col.names = FALSE)

###For Phasing Using VCF
outfile.2b.txt <- outfile.1.txt[as.numeric(outfile.1.txt$V5) > 4 & outfile.1.txt$V9!="HET", ]
outfile.3b.txt <- subset(outfile.2b.txt, select=c(V2, V3, V4 ))

#to remove first base pair of SV for the mask need, subtract 1 from the start position of the SV in the mask file
outfile.3b.txt$V3<-as.numeric(as.character(outfile.3b.txt$V3))-1

write.table(outfile.3b.txt, paste0("CNV_", Isolate, "_",Normal,".forphasing.bed") , quote = FALSE, sep="\t", row.names = FALSE,
            col.names = FALSE)

