######### Load packages ##########
library(data.table)
library("regioneR")
library(BSgenome.Btaurus.UCSC.bosTau9)
genomes<-filterChromosomes(getGenome("BSgenome.Btaurus.UCSC.bosTau9"),organism="bosTau")
## This is an example of enrichment analysis
## This code calculayed the enrichment of CSE(cattle-specific enhancer) for copy number variable regions.

####BED file
system("slopBed -i CSE.bed -g Bos_taurus.chromosome.length -l 50000 -r 50000 | | sort -k 1,1 -k2,2g -k3,3g>CSE.bed.50KB")
Bed<-read.table("CSE.bed.50KB") 
Bed<-toGRanges(Bed, header=FALSE)

Triats <- read.table("ARS_CNV.bed",header=FALSE) 
###ARS_CNV.bed is the file of copy number variable regions(CNVR) which were downloaded from Animal Omics Database (AOD)

Triats <-toGRanges(Triats,header=FALSE)

suppressWarnings(pt <-permTest( A=Bed, B=Triats,ntimes=1000, genome=genomes, 
                               randomize.function=randomizeRegions,
                               evaluate.function=numOverlaps,
                               force.parallel=TRUE,count.once=TRUE))
Observe<-pt[["numOverlaps"]]$observed
Expect<-mean(pt[["numOverlaps"]]$permuted)
write.table(Observe,file=paste("CSE_CNV_observe.txt"))
write.table(Expect,file=paste(CSE_CNV__expect.txt"))
