rm(list = ls())
library("gpclib", lib.loc="../../R_Library")
library("swfscMisc", lib.loc="../../R_Library")
library("strataG", lib.loc="../../R_Library")

dolph.msats <- read.gen.data("../Dolphin Structure msat.csv")
dolph.strata <- read.gen.data("../Dolphin Structure strata.csv")

msat.merge <- merge(dolph.strata, dolph.msats, by = "LabID", all.y = TRUE)
msat.merge <- msat.merge[order(as.numeric(msat.merge$LabID)), ]
msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 10, description = "Ttru Run 4k7")

# Run STRUCTURE
sr <- structure.run(msats, k.range = 1:7, num.k.rep = 10, label = "Run 4k7", delete.files = TRUE, num.cores = 8, 
                    burnin = 100000, numreps = 500000, noadmix = TRUE, freqscorr = FALSE, 
                    pop.prior = NULL)

save.image("Ttru_Run4_srk7.rdata")

# Calculate Evanno metrics
evno <- structure.evanno(sr)
print(evno)

# Run CLUMPP to combine runs for K = 2
clumpp2 <- clumpp.run(sr, k = 2)
strata.num <- as.numeric(dolph.strata$All_Fine)
clumpp2$orig.pop <- strata.num
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp2)

# Run CLUMPP to combine runs for K = 3
clumpp3 <- clumpp.run(sr, k = 3)
strata.num <- as.numeric(dolph.strata$All_Fine)
clumpp3$orig.pop <- strata.num
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp3)

# Run CLUMPP to combine runs for K = 4
clumpp4 <- clumpp.run(sr, k = 4)
strata.num <- as.numeric(dolph.strata$All_Fine)
clumpp4$orig.pop <- strata.num
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp4)

save.image("Ttru_Run4_allk7.rdata")