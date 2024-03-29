library(tidyverse)

args <- commandArgs(TRUE)

output_gwas <- read.table(file = args[1], header = TRUE)
output_gwas <- output_gwas %>% filter(LogPvalue > 5)

args[1] <- gsub(pattern = "Formated_GWAS_", replacement = "", x = args[1])

#### Save output ####
write.table(output_gwas, file = paste0("Signif_markers_", args[1]), quote = FALSE, row.names = FALSE)