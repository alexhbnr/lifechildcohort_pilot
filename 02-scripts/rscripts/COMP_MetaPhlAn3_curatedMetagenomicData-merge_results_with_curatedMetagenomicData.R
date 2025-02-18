library(mia)
library(curatedMetagenomicData)
library(phyloseq)
library(tidyverse)

# Load HMP 2012 tissue info
hmp2012_tissues <- read_tsv(snakemake@params[["hmp_tissues"]])
# Load additional profiles
relab_tbl <- read_tsv(snakemake@input[[1]])

# Create a phyloSeq object from these profiles
## OTU table
otu_mat <- relab_tbl |>
  pivot_wider(names_from = "sample",
              values_from = "relative_abundance",
              values_fill = 0) |>
  mutate(taxon = str_extract(clade_name, "s__.+$")) |>
  column_to_rownames(var = "taxon") |>
  select(-clade_name) |>
  as.matrix()
## Tax table
tax_mat <- relab_tbl %>%
  select(clade_name) %>%
  distinct() %>%
  mutate(clade_name = str_replace_all(clade_name, "[a-z]__", "")) %>%
  pull(clade_name) %>%
  str_split_fixed("\\|", n=7)
rownames(tax_mat) <- rownames(otu_mat)
colnames(tax_mat) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tax_mat[, "Species"] = str_replace_all(tax_mat[, "Species"], "_", " ")
## Sample table
sample_info <- relab_tbl |>
  select(sample) |>
  distinct() |>
  mutate(sample_type = case_when(
      str_sub(sample, 1, 3) == "EXB" ~ "extraction blank",
      str_sub(sample, 1, 3) == "LIB" ~ "library blank",
      str_sub(sample, 1, 3) == "LLP" ~ "LIFE Child cohort plaque",
      str_sub(sample, 1, 3) == "SRS" ~ "HMP subgingival plaque",
    ),
    sample_id = sample) |>
  column_to_rownames(var = "sample_id")
## Create the phyloSeq object
mp3_phy <- phyloseq(otu_table(otu_mat, taxa_are_rows = T),
                    tax_table(tax_mat),
                    sample_data(sample_info))

# Prepare the curatedMetagenomicData HMP samples
## Select samples
hmp_oral_samples <- filter(sampleMetadata, body_site == "oralcavity",
                           study_name == "HMP_2012",
                           number_reads >= 1e7)
## Pull the profiles and convert to phyloSeq
cmd_profiles <- returnSamples(hmp_oral_samples, "relative_abundance",
                              counts = FALSE, rownames = "long")
cmd_profiles_ps <- convertToPhyloseq(cmd_profiles, assay.type = "relative_abundance")
## Make the phyloSeq objects compatible prior to merging
hmp_otu_table <- otu_table(cmd_profiles_ps)
rownames(hmp_otu_table) <- str_extract(rownames(hmp_otu_table), "s__.+$")
hmp_taxa_table <- tax_table(cmd_profiles_ps)
rownames(hmp_taxa_table) <- str_extract(rownames(hmp_taxa_table), "s__.+$")
colnames(hmp_taxa_table) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
hmp_sample_info <- filter(hmp2012_tissues,
                          sample_id %in% hmp_oral_samples$sample_id) |>
  mutate(sample_type = case_when(
          tissue_type == "G_DNA_Supragingival plaque" ~ "HMP supragingival plaque",
          tissue_type == "G_DNA_Buccal mucosa" ~ "HMP buccal mucosa",
          tissue_type == "G_DNA_Tongue dorsum" ~ "HMP tongue dorsum",
          tissue_type == "G_DNA_Subgingival plaque" ~ "HMP subgingival plaque"
      ),
      sample = sample_id) |>
  select(-tissue_type) |>
  column_to_rownames(var = "sample_id")
hmp_phy <- phyloseq(hmp_otu_table,
                    hmp_taxa_table,
                    sample_data(hmp_sample_info))
hmp_phy <- prune_samples(as.logical(!is.na(sample_data(hmp_phy)$sample_type)),
                         hmp_phy)

phy <- merge_phyloseq(mp3_phy, hmp_phy)

save(phy, file = snakemake@output[[1]])
