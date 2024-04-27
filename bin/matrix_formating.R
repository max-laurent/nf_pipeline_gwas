#!/usr/bin/env Rscript

# Import required packages
library("tidyverse")

# Add command line argument
args <- commandArgs(trailingOnly = TRUE)

# Load the phenotypic data and format 
pheno_file <- read.table(file = args[1], header = TRUE, sep = " ")
pheno_file <- pheno_file %>% rename(FID = genotype)

#Load the genotypic data
geno_file <- read.table(args[2], stringsAsFactors = FALSE, header = FALSE)

# Create a phenotypic file designed to be used by FASTLMM
frame <- geno_file[, c(1:2)]
frame <- left_join(frame, pheno_file, by = join_by("V2" == "FID"))
frame[is.na(frame)] <- -9
colnames(frame)[1] <- "FID"
colnames(frame)[2] <- "IID"

#### Save output ####
write.table(frame, file = paste0("FORMATED_", args[1]), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)