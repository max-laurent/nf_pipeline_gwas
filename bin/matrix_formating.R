#!/usr/bin/env Rscript
library("tidyverse")

args <- commandArgs(trailingOnly = TRUE)

pheno_file <- read.table(file = args[1], header = TRUE, sep = " ")

pheno_file <- pheno_file %>% rename(FID = genotype)

geno_file <- read.table(args[2], stringsAsFactors = FALSE, header = FALSE)

frame <- geno_file[, c(1:2)]
frame <- left_join(frame, pheno_file, by = join_by("V2" == "FID"))
frame[is.na(frame)] <- -9

colnames(frame)[1] <- "FID"
colnames(frame)[2] <- "IID"

write.table(frame, file = paste0("FORMATED_", args[1]), sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)