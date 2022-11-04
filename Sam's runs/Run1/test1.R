rm(list = ls())
library("eiaMisc", lib.loc="../../R_Library")
library("strataG", lib.loc="../../R_Library")

load("Ttru_Run3_sr.rdata")

# Run CLUMPP to combine runs for K = 2
clumpp2 <- clumpp.run(sr, k = 2)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp2)

# Run CLUMPP to combine runs for K = 3
clumpp3 <- clumpp.run(sr, k = 3)
#print(clumpp)
# Plot CLUMPP results
structure.plot(clumpp3)

save.image("Ttru_Run3_test.rdata")