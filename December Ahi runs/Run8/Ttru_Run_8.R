rm(list = ls())
library("gpclib", lib.loc="../../R_Library")
library("swfscMisc", lib.loc="../../R_Library")
library("strataG", lib.loc="../../R_Library")
source(file="../CombineTtruLhosData.R")

#strata names: 2-All,3-Fine, 4-Broad, 5-HI.Fine, 6-CNMI_HIArch, 7-CNMI_Other, 8-CNMI_NoHybrids, 9-STRUCTURE_strata
strat.num <- 2
run.label <- "Ttru_Run8"

Ttru.msats <- read.gen.data("../AS183 msat.csv")
Ttru.strata <- read.gen.data("../AS183 strata.csv")
Lhos.msats <- read.gen.data("../AS179 msat.csv")
Lhos.strata <- read.gen.data("../AS179 strata.csv")

comb.data <- combine.Ttru.Lhos.data(Ttru.msats,Ttru.strata,Lhos.msats,Lhos.strata)
msats <- comb.data$msats
strata <- comb.data$strata
loc.col <- ncol(strata)+1

msat.merge <- merge(strata, msats, by = "LabID", all.y = TRUE)
msat.merge <- msat.merge[order(as.numeric(msat.merge$LabID)), ]
msat.gtypes <- gtypes(msat.merge, id.col = 1, strata.col = strat.num, locus.col = loc.col, description = run.label)

#Don't use prior info for CNMI samples
x <- rownames(msat.gtypes$genotypes)
orig.pop <- subset(strata,LabID %in% x)
orig.pop <- orig.pop[order(as.numeric(orig.pop$LabID)),9]
popflag <- rep(1,length(orig.pop))
popflag[which(orig.pop=="CNMI")] <- 0

# Run STRUCTURE
sr <- structure.run(msat.gtypes, k.range = 2:2, num.k.rep = 10, label = run.label, delete.files = FALSE, num.cores = 2, 
                    burnin = 100000, numreps = 500000, noadmix = FALSE, freqscorr = FALSE, 
                    pop.prior = "usepopinfo", popflag=popflag)

save.image(file=paste(run.label,"_sr.rdata",sep=""))

# Calculate Evanno metrics
#evno <- structure.evanno(sr)
#print(evno)

# Run CLUMPP to combine runs for K = 2
clumpp2 <- clumpp.run(sr, k = 2)
#clumpp2$orig.pop <- strata.num
#print(clumpp2)
# Plot CLUMPP results
x <- rownames(msat.gtypes$genotypes)
orig.pop <- subset(strata,LabID %in% x)
clumpp2$orig.pop <- orig.pop[order(as.numeric(orig.pop$LabID)),9]
orig.pop <- sapply(1:nrow(clumpp2), function(x){
  switch(clumpp2[x,3],Lhos=1,CNMI=2,WPac=3,Pelagic=4,HIArch=5)
})
clumpp2$orig.pop <- orig.pop
structure.plot(clumpp2,sort.probs=FALSE)

save.image(file=paste(run.label,"_all.rdata",sep=""))