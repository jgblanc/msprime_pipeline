# Concatenate all results
args=commandArgs(TRUE)

library(data.table)
library(dplyr)

# Assign all arguments
outfile = args[1]

# Loop through files to get results
dfOut <- matrix(NA, nrow = 1, ncol = 9)
for (i in 2:length(args)) {
  print(i)
  filename <- args[i]

  # Extract parameters
  rep <- strsplit(filename, "/")[[1]][3]
  M <- as.numeric(strsplit(strsplit(filename, "/")[[1]][4], "-")[[1]][2])
  theta <- as.numeric(strsplit(strsplit(filename, "/")[[1]][5], "-")[[1]][2])
  L <- as.numeric(strsplit(strsplit(filename, "/")[[1]][6], "-")[[1]][2])

  # Calculate dimensions
  D <- M * L

  # add PC stats
  tmp <- fread(filename)
  b2PC1 <- tmp[1,1]
  b2PC2 <- tmp[1,2]

  # Read in Fst file and extract info
  fstFile <- gsub("b2.txt", "genos.fst.summary", filename)
  dfFst <- as.data.frame(fread(fstFile))
  Fst1 <- mean(dfFst[1:2,3]) # Take average of AvB and AvC
  Fst2 <- dfFst[3,3]

  # Create row to add to output
  outRow <- c(rep, M, theta, L, D, b2PC1, b2PC2, Fst1, Fst2)
  dfOut <- rbind(dfOut, outRow)

}

# Save output
colnames(dfOut) <- c("Rep", "M", "theta", "L", "D", "b2PC1", "b2PC2", "Fst1", "Fst2")
dfOut <- dfOut[2:nrow(dfOut),]
write.table(dfOut, outfile, quote = F, col.names = T, row.names = F)
