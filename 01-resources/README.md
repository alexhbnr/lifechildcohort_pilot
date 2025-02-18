## Overview of the content of the folder `01-resources`

This folder contains files that are unchangeable resources for the analyses of this project.

- `databases.csv`: database input table for running Kraken2 and MetaPhlAn4 using nf-core/taxprofiler
- `hmp2012_curatedMetagenomicData_tissues.tsv`: the information which tissue was samples for each
  entry in the curatedMetagenomicData dataset of the HMP samples
- `hmp2012_subgingival_plaque.tsv`: the list of subgingival plaque samples from the HMP 2012 that
  Irina downloaded
- `kraken_mem.config`: Nextflow input file adapted for running Kraken2 locally with a RAM disk
- `overview_pilot_batch.tsv`: the overview of the plaque samples and the extraction and library
  blanks that are part of the pilot study
