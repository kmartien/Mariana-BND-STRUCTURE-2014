rm(list = ls())
library("gpclib", lib.loc="../../R_Library")
library("swfscMisc", lib.loc="../../R_Library")
library("strataG", lib.loc="../../R_Library")

dolph.msats <- read.gen.data("../Dolphin Structure msat.csv")
dolph.strata <- read.gen.data("../Dolphin Structure strata.csv")

msat.merge <- merge(dolph.strata, dolph.msats, by = "LabID", all.y = TRUE)
msat.merge <- msat.merge[order(as.numeric(msat.merge$LabID)), ]
msats <- gtypes(msat.merge, id.col = 1, strata.col = 4, locus.col = 10, description = "Ttru Run 6")

# Run STRUCTURE
sr <- structure.run(msats, k.range = 2:2, num.k.rep = 10, label = "Run 6", delete.files = FALSE, num.cores = 2, 
                    burnin = 100000, numreps = 500000, noadmix = TRUE, freqscorr = FALSE, 
                    pop.prior = "usepopinfo", popflag = dolph.strata$Pop_Flag_CNMI)

save.image("Ttru_Run6_sr.rdata")


# Run CLUMPP to combine runs for K = 2 
clumpp2 <- clumpp.run(sr, k = 2)
#set stratification for plots
clumpp2$orig.pop <- as.numeric(dolph.strata$All_Fine)
# Plot CLUMPP results
structure.plot(clumpp2)

save.image("Ttru_Run6_all.rdata")
