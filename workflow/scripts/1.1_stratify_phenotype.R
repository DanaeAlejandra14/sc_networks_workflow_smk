log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)

clinical_csv <- snakemake@config[["phenotype"]][["clinical_csv"]]
th <- snakemake@config[["phenotype"]][["thresholds"]]

message("Leyendo datos clinicos: ", clinical_csv)
clinical <- read_csv(clinical_csv, show_col_types = FALSE)

clinical <- clinical %>%
  mutate(phenotype_label = case_when(
    cogdx %in% th$cogdx_control & braaksc >= th$braaksc_ad_min & ceradsc %in% th$ceradsc_ad ~ "AD-NC_ASYM",
    cogdx %in% th$cogdx_control & ceradsc %in% th$ceradsc_control ~ "control",
    cogdx %in% th$cogdx_ad_symptomatic & braaksc >= th$braaksc_ad_min & ceradsc %in% th$ceradsc_ad ~ "AD-NC_SYM",
    cogdx %in% th$cogdx_mci ~ "MCI",
    TRUE ~ NA_character_
  ))

out <- clinical %>% select(individualID, phenotype_label)

message("Distribucion de fenotipos:")
message(paste(capture.output(print(table(out$phenotype_label, useNA = "ifany"))), collapse = "\n"))

write_csv(out, snakemake@output[[1]])
message("Tabla de fenotipos guardada: ", nrow(out), " individuos")

sink(type = "message"); sink()