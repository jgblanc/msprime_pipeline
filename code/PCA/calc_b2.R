## Compute b2
args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_b2.R")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

infile = args[1]
outfile = args[2]
M = as.numeric(args[3])
L = as.numeric(args[4])
theta = as.numeric(args[5])
fst1 = as.numeric(args[6])
fst2 = as.numeric(args[7])

# Define population sizes with rounding
size_B <- round(theta * M * 0.5)
size_C <- round(theta * M * 0.5)
size_A <- M - size_B - size_C

# Construct population eigenvectors
pop_e1 <- c(rep(-theta, size_A), rep(1 - theta, size_B + size_C)) / (theta * sqrt(((theta^2 + 1)/theta) - 1))
pop_e2 <- c(rep(0, size_A), rep(1, size_B), rep(-1, size_C)) / sqrt(theta)


# Read in sample eigenvectors
df <- fread(infile, header=FALSE)
sample_e1 <- df[,1]
sample_e2 <- df[,2]

# Compute b2
b2_e1 <- cor(sample_e1, pop_e1)^2
b2_e2 <- cor(sample_e2, pop_e2)^2

# Format and save output
dfOut <- data.frame(b2_e1 = b2_e1, b2_e2 = b2_e2)
fwrite(dfOut, outfile, quote = F, row.names = F, sep = "\t")
