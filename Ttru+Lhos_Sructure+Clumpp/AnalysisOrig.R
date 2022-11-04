rm(list = ls())
library(strataG)

data(dolph.strata)
data(dolph.msats)
msat.merge <- merge(dolph.strata, dolph.msats, by = "ids", all.y = TRUE)
msats <- gtypes(msat.merge, id.col = 1, strata.col = 3, locus.col = 5, description = "Dolphin msat example")

# Run STRUCTURE
sr <- structure.run(msats, k.range = 1:6, num.k.rep = 10, delete.files = TRUE)

# Calculate Evanno metrics
evno <- structure.evanno(sr)
print(evno)

# Run CLUMPP to combine runs for K = 2
clumpp <- clumpp.run(sr, k = 2)
print(clumpp)

# Plot CLUMPP results
structure.plot(clumpp)

## End(Not run)