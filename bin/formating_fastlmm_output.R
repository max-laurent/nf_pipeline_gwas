#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(TRUE)

output_gwas <- read.table(file = args[1], header = TRUE)

#### Add relevant columns to the df ####
output_gwas <- output_gwas %>% 
  mutate(R2_LR = 1 - exp(-2*(AltLogLike - NullLogLike)/N), 
         LogPvalue = -log10(Pvalue),
         LogQvalue = - log10(Qvalue)) %>% 
  select(-Pvalue, -Qvalue, -AltLogLike, -NullLogLike)

#### Add the MAF to the GWAS output and filter markers with MAF less than 0.05 ####
MAF <- read.table(file = args[2], header = TRUE)
MAF <- MAF %>% filter(MAF > 0.05 & SNP %in% output_gwas$SNP)
output_gwas <- output_gwas %>% filter(SNP %in% MAF$SNP)
output_gwas <- output_gwas %>% left_join(., MAF[,c("SNP", "MAF", "NCHROBS")], by = join_by("SNP"))
rm(MAF)

#### Add the information about the reference alleles ####
geno_file <- read.table(file = args[3], header = TRUE)
colnames(geno_file) <- c("Chromosome", "SNP", "GeneticDistance", "Position", "Allele1", "Allele2", "Allele1_Dose", "Allele2_Dose")
geno_file <- geno_file %>% filter(SNP %in% output_gwas$SNP)

markers_file <- readRDS(file = args[4])
markers_file <- t(markers_file)
markers_file <- markers_file[output_gwas$SNP, "B73"] %>% as.data.frame() %>% rownames_to_column(., var = "SNP")
colnames(markers_file)[2] <- "RefAllele_dose"

geno_file <- geno_file %>% left_join(., markers_file, by = join_by("SNP")) %>% select(-Chromosome,-GeneticDistance,-Position)
rm(markers_file)
output_gwas <- output_gwas %>% left_join(., geno_file, by = join_by("SNP"))

output_gwas <- output_gwas %>% mutate(SNP_Weight_RefAllele = ifelse(RefAllele_dose == Allele1_Dose, SNPWeight, -SNPWeight),
                                      RefAllele_letter = ifelse(RefAllele_dose == Allele1_Dose, Allele1,
                                                                ifelse(RefAllele_dose == Allele2_Dose, Allele2, paste(Allele1, Allele2, sep = "/"))))
rm(geno_file)

#### Convert position from v2 to v4 ####
info <- read.table(file = args[5])
info <- info %>% filter(snp.name %in% output_gwas$SNP) %>% select(snp.name, Chromosome, PositionPhys)
output_gwas <- output_gwas %>% 
  select(-Chromosome, -Position) %>% 
  inner_join(., info, by = join_by("SNP" == "snp.name")) %>%
  rename(Position = PositionPhys)

rm(info)

#### Save output ####
write.table(output_gwas, file = paste0("Formated_", args[1]), quote = FALSE, row.names = FALSE)