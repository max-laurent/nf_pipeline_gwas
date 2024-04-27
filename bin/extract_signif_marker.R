#!/usr/bin/env Rscript

# Import required packages
library(tidyverse)

# Add command line argument
args <- commandArgs(TRUE)

# Load the data
output_gwas <- read.table(file = args[1], header = TRUE)

# Extract significant markers
output_gwas <- output_gwas %>% filter(LogPvalue > 5)

args[1] <- gsub(pattern = "Formated_GWAS_", replacement = "", x = args[1])

#### Save output ####
write.table(output_gwas, file = paste0("Signif_markers_", args[1]), quote = FALSE, row.names = FALSE)