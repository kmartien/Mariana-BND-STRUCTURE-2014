rm(list = ls())
library("gpclib", lib.loc="../../R_Library")
library("swfscMisc", lib.loc="../../R_Library")
library("strataG", lib.loc="../../R_Library")
source(file="../CombineTtruLhosData.R")

#strata names: All,Fine, Broad, HI.Fine, CNMI_HIArch, CNMI_Other, CNMI_NoHybrids, STRUCTURE_strata
strat.num <- 9
run.label <- "Ttru_Run1"

Ttru.msats <- read.gen.data("../AS183 msat.csv")
Ttru.strata <- read.gen.data("../AS183 strata.csv")
Lhos.msats <- read.gen.data("../AS179 msat.csv")
Lhos.strata <- read.gen.data("../AS179 strata.csv")
#Ttru.msats <- read.gen.data("AS183 msat.csv")
#Ttru.strata <- read.gen.data("AS183 strata.csv")
#Lhos.msats <- read.gen.data("AS179 msat.csv")
#Lhos.strata <- read.gen.data("AS179 strata.csv")

comb.data <- combine.Ttru.Lhos.data(Ttru.msats,Ttru.strata,Lhos.msats,Lhos.strata)
msats <- comb.data$msats
strata <- comb.data$strata
loc.col <- ncol(strata)+1

msat.merge <- merge(strata, msats, by = "LabID", all.y = TRUE)
msat.merge <- msat.merge[order(as.numeric(msat.merge$LabID)), ]
msat.gtypes <- gtypes(msat.merge, id.col = 1, strata.col = strat.num, locus.col = loc.col, description = run.label)

#for Run1, remove Lhos
strata.names <- attr(msat.gtypes,"strata.names")
strata.names <- strata.names[!is.na(strata.names)]
strata.to.keep <- as.vector(strata.names[-which(strata.names=="Lhos")])
msat.gtypes <- subset(msat.gtypes, strata=strata.to.keep)

# Run STRUCTURE
sr <- structure.run(msat.gtypes, k.range = 1:6, num.k.rep = 10, label = run.label, delete.files = FALSE, num.cores = 2, 
                    burnin = 100000, numreps = 500000, noadmix = FALSE, freqscorr = TRUE, 
                    pop.prior = NULL)

save.image(file=paste(run.label,"_sr.rdata",sep=""))
# Calculate Evanno metrics
#evno <- structure.evanno(sr)
#print(evno)

#strata.num <- na.omit(strata[strat.num])
#strata.num <- strata.num[1:166] + 5

# Run CLUMPP to combine runs for K = 2
clumpp2 <- clumpp.run(sr, k = 2)
#clumpp2$orig.pop <- strata.num
#print(clumpp2)
# Plot CLUMPP results
orig.pop <- sapply(1:nrow(clumpp2), function(x){
  switch(clumpp2[x,3],CNMI=1,WPac=2,Pelagic=3,HIArch=4)
})
clumpp2$orig.pop <- orig.pop
structure.plot(clumpp2,sort.probs=FALSE)

# Run CLUMPP to combine runs for K = 3
clumpp3 <- clumpp.run(sr, k = 3)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp3)

# Run CLUMPP to combine runs for K = 4
clumpp4 <- clumpp.run(sr, k = 4)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp4)

# Run CLUMPP to combine runs for K = 5
clumpp5 <- clumpp.run(sr, k = 5)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp5)

# Run CLUMPP to combine runs for K = 6
clumpp6 <- clumpp.run(sr, k = 6)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp6)

save.image(file=paste(run.label,"_all.rdata",sep=""))