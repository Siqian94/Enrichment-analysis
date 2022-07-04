######### Load packages ##########
library(data.table)
library("regioneR")
library(BSgenome.Btaurus.UCSC.bosTau9)
genomes<-filterChromosomes(getGenome("BSgenome.Btaurus.UCSC.bosTau9"),organism="bosTau")

####BED file
Bed<-read.table("CSE.bed.50KB")
Bed<-toGRanges(Bed, header=FALSE)

Triats <- read.table("ARS_CNV.bed",header=FALSE)
Triats <-toGRanges(Triats,header=FALSE)

suppressWarnings(pt <-permTest( A=Bed, B=Triats,ntimes=1000, genome=genomes, 
                               randomize.function=randomizeRegions,
                               evaluate.function=numOverlaps,
                               force.parallel=TRUE,count.once=TRUE))
Observe<-pt[["numOverlaps"]]$observed
Expect<-mean(pt[["numOverlaps"]]$permuted)
write.table(Observe,file=paste("CSE_CNV_observe.txt"))
write.table(Expect,file=paste(CSE_CNV__expect.txt"))
