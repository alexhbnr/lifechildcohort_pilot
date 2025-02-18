## Overview of the folder `02-scripts`

This folder contains all the workflows and scripts that are necessary to conduct the experiments of
this project.

### `PREP`: data preparation

- `PREP_link_sequencing_data.Snakefile`: link the de-multiplexed FastQ files into the current file
  structure

### `COMP`: compositional analysis

- `COMP_taxprofiling.Snakefile`: taxonomically profile the plaque samples and blanks with Kraken2
  and MetaPhlAn4 using nf-core/taxprofiler
