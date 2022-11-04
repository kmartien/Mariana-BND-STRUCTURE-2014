rm(list = ls())
library("gpclib", lib.loc="../../R_Library")
library("swfscMisc", lib.loc="../../R_Library")
library("strataG", lib.loc="../../R_Library")

dolph.msats <- read.gen.data("../Dolphin Structure msat.csv")
dolph.strata <- read.gen.data("../Dolphin Structure strata.csv")


msat.merge <- merge(dolph.strata, dolph.msats, by = "LabID", all.y = TRUE)
msat.merge <- msat.merge[order(as.numeric(msat.merge$LabID)), ]
msats <- gtypes(msat.merge, id.col = 1, strata.col = 7, locus.col = 10, description = "Ttru Run 1")

# Run STRUCTURE
sr <- structure.run(msats, k.range = 1:6, num.k.rep = 10, label = "Run 1", delete.files = FALSE, num.cores = 2, 
                    burnin = 100000, numreps = 500000, noadmix = FALSE, freqscorr = TRUE, 
                    pop.prior = NULL)

save.image("Ttru_Run1_sr.rdata")

# Calculate Evanno metrics
evno <- structure.evanno(sr)
print(evno)

strata.num <- na.omit(as.numeric(dolph.strata$AS183_Ttru))
strata.num <- strata.num[1:166] + 5

# Run CLUMPP to combine runs for K = 2
clumpp2 <- clumpp.run(sr, k = 2)
clumpp2$orig.pop <- strata.num
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp2)

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

save.image("Ttru_Run1_all.rdata")