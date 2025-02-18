## The overview of the folder `05-results`

This folder contains the results of the analyses of this project.

### `PREP`: results of analyses related to the data preparation

  - `PREP_seqdata_processing_nReads_pilot.tsv`: overview of the number of sequences that were
    retained throughout the different pre-processing steps
  - `PREP_seqdata_processing_insertsizes_pilot.tsv`: overview of the insert size distribution of the
    samples

### `COMP`: results of analyses related to the taxonomic profiling

  - `COMP_TaxProfiler_Kraken2_NCBIRefSeq_universal_240912_pilot.tsv`: the taxonomic profiles of the
    pilot samples obtained from Kraken2 with a NCBI RefSeq databases containing microbes
  - `COMP_TaxProfiler_MetaPhlAnv4_Jun23_pilot.tsv`: the taxonomic profiles of the pilot samples
    obtained from MetaPhlAn with its database v4 from June 2023
  - `COMP_MetaPhlAn3_pilot_HMP2012subgplaque.tsv`: the taxonomic profiles of the pilot samples and
    the subgingival plaque samples obtained from MetaPhlAn with its database v3 from January 2019
  - `COMP_MetaPhlAn3_pilot_HMP2012_phyloSeq.RData`: a phyloSeq object of the pilot samples and the
    HMP samples from the oral sites subgingival plaque, supragingival plaque, buccal mucosa, and
    tongue dorsum obtained from MetaPhlAn with its database v3 from January 2019
