################################################################################
# Project: LIFE Child cohort - plaque samples
# Part: Taxonomic profiling
# Step: Profile samples with MetaPhlAn3 to compare against the HMP samples
#
# Dependent on:
#   - COMP_taxprofiling.Snakefile
#
# Alex Huebner, 12/02/25
################################################################################

from glob import glob
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES = {}
# Pilot
pilot = pd.read_csv("04-analysis/taxprofiler/pilot.samplesheet.csv", sep=",")
for s in pilot['sample'].unique():
    SAMPLES[s] = ('pilot', f'/mnt/archgen/microbiome_calculus/lifeplaque/04-analysis/taxprofiler/results_pilot/run_merging/{s}')
# Subgingival plaque
subplq = pd.read_csv("01-resources/hmp2012_subgingival_plaque.tsv", sep="\t")
for s in subplq['sample_id'].unique():
    SAMPLES[s] = ('subplq', f'/mnt/archgen/microbiome_calculus/Cameroon_plaque/03-preprocessing/all_data_combined/hmp_pq_eager2/samtools/filter/{s}')
################################################################################

#### Auxilliary functions ######################################################

def return_reads(wildcards):
    if SAMPLES[wildcards.sample][0] == "pilot":
        p = ",".join([f"{SAMPLES[wildcards.sample][1]}_{i}.merged.fastq.gz" for i in [1, 2]])
    else:
        p = f"{SAMPLES[wildcards.sample][1]}.unmapped.fastq.gz"
    return p

################################################################################

rule all:
    input:
        "05-results/COMP_MetaPhlAn3_pilot_HMP2012subgplaque.tsv",
        "05-results/COMP_MetaPhlAn3_pilot_HMP2012_phyloSeq.RData"

#### MetaPhlAn3 ################################################################

rule metaphlan3:
    output:
        profile = "04-analysis/metaphlan3/{sample}.metaphlan.profile.txt",
        sam = "04-analysis/metaphlan3/{sample}.metaphlan.sam",
        bowtie2out = temp("04-analysis/metaphlan3/{sample}.metaphlan.bowtie2out")
    message: "Run MetaPhlAn3 with default settings for sample {wildcards.sample}"
    container: "https://depot.galaxyproject.org/singularity/metaphlan%3A3.1.0--pyhb7b1952_0"
    resources:
        mem = 16
    params:
        db = "/mnt/archgen/users/huebner/refdbs/metaphlan/mpa_v31_CHOCOPhlAn_201901",
        fastq = lambda wildcards: return_reads(wildcards)
    threads: 8
    shell:
        """
        metaphlan \
            {params.fastq} \
            --input_type fastq \
            --force \
            --index mpa_v31_CHOCOPhlAn_201901 \
            --bowtie2db {params.db} \
            --bowtie2out {output.bowtie2out} \
            --ignore_eukaryotes \
            -t rel_ab_w_read_stats \
            --sample_id {wildcards.sample} \
            -o {output.profile} \
            -s {output.sam} \
            --nproc {threads}
        """

rule summarise_metaphlan3:
    input:
        expand("04-analysis/metaphlan3/{sample}.metaphlan.profile.txt", sample=SAMPLES)
    output:
        "05-results/COMP_MetaPhlAn3_pilot_HMP2012subgplaque.tsv"
    message: "Summarise the MetaPhlAn3 results into a single table"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    run:
        mp3_profiles = []
        for fn in input:
            sample = os.path.basename(fn).replace(".metaphlan.profile.txt", "")
            tbl = (pd.read_csv(fn, sep="\t", skiprows=4, usecols=[0, 2])
                .rename({'#clade_name': 'clade_name'}, axis=1)
                .query('clade_name.str.contains("s__")')
                .assign(sample=sample)
            )
            mp3_profiles.append(tbl)

        (pd.concat(mp3_profiles)[['sample', 'clade_name', 'relative_abundance']]
         .sort_values(['sample', 'relative_abundance'], ascending=[True, False])
         .to_csv(output[0], sep="\t", index=False)
        )

rule merge_results_with_curatedMetagenomicData:
    input:
        "05-results/COMP_MetaPhlAn3_pilot_HMP2012subgplaque.tsv"
    output:
        "05-results/COMP_MetaPhlAn3_pilot_HMP2012_phyloSeq.RData"
    message: "Merge the results with the HMP 2012 oral samples from the curatedMetagenomeData package"
    conda: "ENVS_Renv.yaml"
    resources:
        mem = 4,
        cores = 1
    params:
        hmp_tissues = "01-resources/hmp2012_curatedMetagenomicData_tissues.tsv"
    threads: 1
    script:
        "rscripts/COMP_MetaPhlAn3_curatedMetagenomicData-merge_results_with_curatedMetagenomicData.R"
