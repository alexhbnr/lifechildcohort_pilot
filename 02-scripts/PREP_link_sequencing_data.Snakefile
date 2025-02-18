################################################################################
# Project: LIFE Child cohort - plaque samples
# Part: Data preparation
# Step: Link the sequencing data into the folder structure
#
# Dependent on:
#
# Alex Huebner, 29/01/25
################################################################################

from glob import glob
import os

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SEQRUNS = {'pilot': '01-resources/overview_pilot_batch.tsv'}
################################################################################


rule all:
    input:
        expand("03-data/rawdata/{seqrun}.done", seqrun=SEQRUNS)

rule link_fastqs:
    output:
        touch("03-data/rawdata/{seqrun}.done")
    message: "Link the FastQ files of the sequencing run: {wildcards.seqrun}"
    resources:
        mem = 2,
        cores = 1
    params:
        seq_tbl = lambda wildcards: SEQRUNS[wildcards.seqrun],
        dir = "/mnt/archgen/data_releases",
        outdir = "03-data/raw_data/{seqrun}"
    threads: 1
    run:
        os.makedirs(params.outdir, exist_ok=True)
        tbl = pd.read_csv(params.seq_tbl, sep="\t")
        for s in tbl.itertuples():
            year = f"20{s.seqrun_id[:2]}"
            fastqs = glob(f"{params.dir}/{year}/{s.seqrun_id}/{s.sample_id}*/*.fastq.gz")
            for fastq in fastqs:
                if not os.path.islink(f"{params.outdir}/{os.path.basename(fastq)}"):
                    os.symlink(fastq, f"{params.outdir}/{os.path.basename(fastq)}")
