## Overview of the folder `02-scripts`

This folder contains all the workflows and scripts that are necessary to conduct the experiments of
this project.

### `PREP`: data preparation

- `PREP_link_sequencing_data.Snakefile`: link the de-multiplexed FastQ files into the current file
  structure

### `COMP`: compositional analysis

- `COMP_taxprofiling.Snakefile`: taxonomically profile the plaque samples and blanks with Kraken2
  and MetaPhlAn4 using nf-core/taxprofiler
- `COMP_MetaPhlAn3_curatedMetagenomicData.Snakefile`: merge the MetaPhlAn3 taxonomic profiles of the
  pilot samples with the HMP oral cavity samples from the curatedMetagenomicData R package

### `STRA`: strain analysis

- `STRA_pilot_top3_strains.Snakefile`: strain analysis of the top 3 most abundant species
