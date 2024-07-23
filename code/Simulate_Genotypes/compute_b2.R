# Compute correlation between population and sample PCs
args=commandArgs(TRUE)

library(data.table)
library(dplyr)

# Assign all arguments
evec_file = args[1]
pop_file = args[2]
outfile = args[3]

# Read in PCs and pop info and joing
dfE <- fread(evec_file, header = T)
dfPop <- fread(pop_file, header = T)
colnames(dfPop)[1] <- "#FID"
df <- inner_join(dfE, dfPop)

# Manually create population PCs
df <- df %>% mutate(popPC1 = case_when(POP == "A" ~ 1, POP == "B" ~ 0, POP == "C" ~ 0),
                    popPC2 = case_when(POP == "A" ~ 0, POP == "B" ~ 1, POP == "C" ~ -1))

# Calculate b2
b2_PC1 <- cor(df$PC1, df$popPC1)^2
b2_PC2 <- cor(df$PC2, df$popPC2)^2

# Set up output table
dfOut <- as.data.frame(matrix(NA, nrow = 1, ncol = 2))
dfOut[1, ] <- c(b2_PC1, b2_PC2)
colnames(dfOut) <- c("PC1", "PC2")

# Save output
write.table(dfOut, outfile, quote = F, col.names = T, row.names = F)
