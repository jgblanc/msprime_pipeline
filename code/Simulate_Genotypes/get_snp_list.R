# Sample L total SNPs
args=commandArgs(TRUE)

library(data.table)
library(dplyr)

# Assign all arguments
pvar = args[1]
outfile = args[2]
nSNP = as.numeric(args[3])

# Test panel.freq file
df <- fread(pvar, header = T)
print(paste0("The numer of total SNPs is ", nrow(df)))

# Randomly sample only nSNP number of SNPs
df <- df %>% sample_n(nSNP) %>% select(ID)

write.table(df, outfile, quote = F, col.names = F, row.names = F)
