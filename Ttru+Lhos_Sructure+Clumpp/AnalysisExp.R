rm(list = ls())
library("eiaMisc", lib.loc="~/R/win-library/3.1")
library("strataG", lib.loc="~/R/win-library/3.1")

dolph.msats <- read.csv("Dolphin Structure msat.csv")
dolph.strata <- read.csv("Dolphin Structure strata.csv")

msat.merge <- merge(dolph.strata, dolph.msats, by = "LabID", all.y = TRUE)
msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 8, description = "Dolphin msat example")

# Run STRUCTURE
sr <- structure.run(msats, k.range = 1:2, num.k.rep = 1, delete.files = TRUE)

# Calculate Evanno metrics
evno <- structure.evanno(sr)
print(evno)

# Run CLUMPP to combine runs for K = 2
clumpp <- clumpp.run(sr, k = 2)
print(clumpp)

# Plot CLUMPP results
structure.plot(clumpp)

## End(Not run)