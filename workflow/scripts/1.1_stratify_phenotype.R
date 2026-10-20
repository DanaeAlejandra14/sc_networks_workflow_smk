log <- file(snakemake@log[[1]], open = "wt")
sink(log); sink(log, type = "message")

library(tidyverse)

clinical_csv     <- snakemake@config[["phenotype"]][["clinical_csv"]]
biospecimen_csv  <- snakemake@input[["biospecimen"]]
th <- snakemake@config[["phenotype"]][["thresholds"]]

message("Reading clinical data: ", clinical_csv)
clinical <- read_csv(clinical_csv, show_col_types = FALSE)

message("Reading biospecimen data: ", biospecimen_csv)
biospecimen <- read_csv(biospecimen_csv, show_col_types = FALSE)

# Confirm that there are no duplicated individualIDs in the clinical data and that the AD and control CERAD scores do not overlap
stopifnot(!anyDuplicated(clinical$individualID))
stopifnot(length(intersect(unlist(th$ceradsc_ad), unlist(th$ceradsc_control))) == 0)


#Stratify phenotype based on cogdx, braaksc and ceradsc, using the thresholds defined in the config file
clinical <- clinical %>%
  mutate(phenotype_label = case_when(
    cogdx %in% th$cogdx_control & braaksc >= th$braaksc_ad_min & ceradsc %in% th$ceradsc_ad ~ "AD-NC_ASYM",
    cogdx %in% th$cogdx_control & ceradsc %in% th$ceradsc_control ~ "control",
    cogdx %in% th$cogdx_ad_symptomatic & braaksc >= th$braaksc_ad_min & ceradsc %in% th$ceradsc_ad ~ "AD-NC_SYM",
    cogdx %in% th$cogdx_mci ~ "MCI",
    TRUE ~ NA_character_

  ))

sc_individuals <- biospecimen %>%
  filter(assay == "scrnaSeq", is.na(exclude) | exclude != TRUE) %>%
  distinct(individualID) %>%
  pull(individualID)

out <- clinical %>%
  filter(individualID %in% sc_individuals, !is.na(phenotype_label)) %>%
  select(individualID, phenotype_label, msex, educ, pmi, braaksc,
         apoe_genotype, ceradsc, cogdx, dcfdx_lv)

message("Phenotype distribution (individuals with single-nucleus data only):")
message(paste(capture.output(print(table(out$phenotype_label, useNA = "ifany"))), collapse = "\n"))

write_csv(out, snakemake@output[[1]])
message("Final table saved: ", nrow(out), " individuals")

sink(type = "message"); sink()