################################################################################
# Project: LIFE Child cohort - plaque samples
# Part: Taxonomic profiling
# Step: Taxonomic profiling with nf-core/taxprofiler
#
# Using nf-core/taxprofiler, I will initally process the metagenomic sequencing
# data and subsequently profile them taxonomically using Kraken2 and
# MetaPhlAn4.
#
# Dependent on:
#   - PREP_link_sequencing_data.Snakefile
#
# Alex Huebner, 29/01/25
################################################################################

from glob import glob
import json
from math import ceil
import os
import re

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
        expand("04-analysis/taxprofiler/results_{seqrun}/multiqc/multiqc_report.html", seqrun=SEQRUNS),
        expand("05-results/COMP_TaxProfiler_MetaPhlAnv4_Jun23_{seqrun}.tsv", seqrun=SEQRUNS),
        expand("05-results/COMP_TaxProfiler_Kraken2_NCBIRefSeq_universal_240912_{seqrun}.tsv", seqrun=SEQRUNS),
        expand("05-results/PREP_seqdata_processing_nReads_{seqrun}.tsv", seqrun=SEQRUNS),
        expand("05-results/PREP_seqdata_processing_insertsizes_{seqrun}.tsv", seqrun=SEQRUNS)

#### Prepare auxilliary files for assigning taxonomy from ids ##################

rule download_taxdump:
    output:
        "tmp/taxdump.tar.gz"
    message: "Download NCBI taxonomy information"
    params:
        url = "https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    shell:
        """
        wget -O {output} {params.url}
        """

rule extract_ncbi_tax_dumps:
    input:
        "tmp/taxdump.tar.gz"
    output:
        "tmp/taxdump/merged.dmp"
    message: "Extract taxdump files required for taxpasta"
    params:
        dir = "tmp/taxdump"
    shell:
        """
        mkdir -p {params.dir}
        for f in names nodes merged; do
            tar -C {params.dir} -xzf {input} ${{f}}.dmp
        done
        """

################################################################################

#### nf-core/taxprofiler #######################################################

rule prepare_samplesheet:
    output:
        "04-analysis/taxprofiler/{seqrun}.samplesheet.csv"
    message: "Prepare the samplesheet for nf-core/taxprofiler: {wildcards.seqrun}"
    resources:
        mem = 2,
        cores = 1
    params:
        dir = "03-data/raw_data/{seqrun}"
    threads: 1
    run:
        fastqs = glob(f"{params.dir}/*_R1_001.fastq.gz")
        runs = []
        for fq in fastqs:
            sample, lane = re.search(r'([A-Z]+[0-9]+\.[A-Z]+[0-9]+)\.SG[0-9]\.[0-9]_S[0-9]_(L[0-9]+)_R1_001.fastq.gz',
                                     os.path.basename(fq)).groups()
            runs.append((sample, lane, "ILLUMINA", fq, fq.replace("_R1_", "_R2_")))
        samplesheet = pd.DataFrame(runs, columns=['sample', 'run_accession',
                                                  'instrument_platform',
                                                  'fastq_1', 'fastq_2'])
        samplesheet['fasta'] = ""
        samplesheet.to_csv(output[0], sep=",", index=False)

rule run_taxprofiler:
    input:
        tsv = "04-analysis/taxprofiler/{seqrun}.samplesheet.csv",
        dmp = "tmp/taxdump/merged.dmp"
    output:
        "04-analysis/taxprofiler/results_{seqrun}/multiqc/multiqc_report.html"
    message: "Run nf-core/taxprofiler: {wildcards.seqrun}"
    params:
        db_conf = "01-resources/databases.csv",
        outdir = "04-analysis/taxprofiler/results_{seqrun}",
        hostref = "03-data/refgenomes/chm13.draft_v2.0.fasta",
        hostref_dir = "03-data/refgenomes",
        taxdump_dir = "tmp/taxdump",
        custom_config = "01-resources/kraken_mem.config",
        workdir = "tmp/taxprofiler_{seqrun}"
    shell:
        """
        nextflow run nf-core/taxprofiler -r 1.2.2 -profile eva,archgen,custom \
            --input {input.tsv} \
            --databases {params.db_conf} \
            --outdir {params.outdir} \
            --perform_shortread_qc \
            --perform_shortread_hostremoval \
            --hostremoval_reference {params.hostref} \
            --shortread_hostremoval_index {params.hostref_dir} \
            --perform_runmerging \
            --save_runmerged_reads \
            --run_kraken2 \
            --run_metaphlan \
            -c {params.custom_config} \
            -w {params.workdir}
        """

################################################################################

#### Post-processing ###########################################################

rule summarise_metaphlan:
    input:
        "04-analysis/taxprofiler/results_{seqrun}/multiqc/multiqc_report.html"
    output:
        "05-results/COMP_TaxProfiler_MetaPhlAnv4_Jun23_{seqrun}.tsv"
    message: "Summarise the taxonomic profiles inferred by MetaPhlAn4: {wildcards.seqrun}"
    resources:
        mem = 2,
        cores = 1
    params:
        tbl = "04-analysis/taxprofiler/results_{seqrun}/metaphlan/metaphlan_metaphlanv4_Oct22_combined_reports.txt"
    threads: 1
    run:
        mp4 = (pd.read_csv(params.tbl, sep="\t", skiprows=1)
               .query('clade_name.str.contains("s__") and ~clade_name.str.contains("t__")'))
        mp4.columns = [
            c.replace("_metaphlanv4_Oct22.metaphlan", "") for c in mp4.columns
        ]
        col_order = mp4.columns.sort_values().tolist()[:-1]

        (mp4.loc[:, ['clade_name'] + col_order]
         .sort_values(col_order[0], ascending=False)
          .to_csv(output[0], sep="\t", index=False, float_format="%.6f")
        )

rule summarise_kraken2:
    input:
        multiqc = "04-analysis/taxprofiler/results_{seqrun}/multiqc/multiqc_report.html",
        taxdump = "tmp/taxdump/merged.dmp"
    output:
        "05-results/COMP_TaxProfiler_Kraken2_NCBIRefSeq_universal_240912_{seqrun}.tsv"
    message: "Summarise the taxonomic profiles inferred by Kraken2: {wildcards.seqrun}"
    container: "https://depot.galaxyproject.org/singularity/taxpasta%3A0.7.0--pyhdfd78af_0"
    resources:
        mem = 8,
        cores = 1
    params:
        dmpdir = "tmp/taxdump",
        k2dir = "04-analysis/taxprofiler/results_{seqrun}/kraken2/k2_ncbirefseq_universal_240912"
    threads: 1
    shell:
        """
        taxpasta merge -p kraken2 -o {output} \
            --output-format TSV \
            --wide \
            --add-name --add-rank \
            --add-lineage --add-rank-lineage \
            --taxonomy {params.dmpdir} \
            $(find {params.k2dir} -type f -name "*.report.txt")
        sed -i "s/_k2_ncbirefseq_universal_240912.kraken2.kraken2.report//g" {output}
        """

rule summarise_readcounts:
    input:
        "04-analysis/taxprofiler/results_{seqrun}/multiqc/multiqc_report.html"
    output:
        nreads = "05-results/PREP_seqdata_processing_nReads_{seqrun}.tsv",
        insertsizes = "05-results/PREP_seqdata_processing_insertsizes_{seqrun}.tsv",
    message: "Summarise the statistics of the read processing issued by FastP and BowTie2: {wildcards.seqrun}"
    resources:
        mem = 4,
        cores = 1
    params:
        fastp = "04-analysis/taxprofiler/results_{seqrun}/fastp",
        bowtie2 = "04-analysis/taxprofiler/results_{seqrun}/multiqc/multiqc_data/bowtie2_pe_plot.txt"
    threads: 1
    run:
        # Read the FastP Json files
        fastp_reads = []
        fastp_duplication = []
        fastp_insertsizes = []
        for fn in glob(f"{params.fastp}/*.json"):
            sample = os.path.basename(fn).replace(".fastp.json", "")
            fastp = json.load(open(fn, "rt"))
            fastp_reads.append(pd.DataFrame.from_dict(fastp['filtering_result'], orient="index")
             .transpose()
             .assign(sample=sample)
             .assign(raw_reads=fastp['summary']['before_filtering']['total_reads'])
             .assign(final_reads=fastp['summary']['after_filtering']['total_reads'])
             )
            fastp_duplication.append((sample, fastp['duplication']['rate']))
            fastp_histogram = pd.DataFrame.from_dict({'length': list(range(len(fastp['insert_size']['histogram']))),
                                                      'count': fastp['insert_size']['histogram']})
            fastp_histogram = pd.concat([fastp_histogram,
                       pd.DataFrame.from_dict({'length': [len(fastp['insert_size']['histogram'])],
                                               'count': [fastp['insert_size']['unknown']]})
            ])
            fastp_insertsizes.append(fastp_histogram.assign(sample=sample))
        
        # Read the BowTie2 summary
        bt2 = pd.read_csv(params.bowtie2, sep="\t").rename({'Sample': 'sample'}, axis=1)
        bt2['non_host_reads'] = (bt2['PE neither mate aligned'] * 2).astype(int)

        # Summarise the number of reads
        nreads = (pd.concat(fastp_reads)
         .merge(pd.DataFrame(fastp_duplication, columns=['sample', 'duplication_rate']),
                             how="left", on="sample")
         .merge(bt2[['sample', 'non_host_reads']], how="left", on="sample")
         )
        nreads['sample'] = nreads['sample'].str.split("_L").str[0]
        (nreads.groupby(['sample'], as_index=False)
         .agg({
            'passed_filter_reads': 'sum',
            'low_quality_reads': 'sum',
            'too_many_N_reads': 'sum',
            'low_complexity_reads': 'sum',
            'too_short_reads': 'sum',
            'too_long_reads': 'sum',
            'raw_reads': 'sum',
            'final_reads': 'sum',
            'duplication_rate': 'mean',
            'non_host_reads': 'sum',
        })
         .iloc[:, [0, 7, 8, 9, 1, 10, 2, 3, 4, 5, 6]]
        .to_csv(output['nreads'], sep="\t", index=False, float_format="%.2f")
        )

        # Summarise the insert sizes
        insert_size = pd.concat(fastp_insertsizes)
        insert_size['sample'] = insert_size['sample'].str.split("_L").str[0]
        insert_size = (insert_size.groupby(['sample', 'length'],
                                          as_index=False)['count'].agg('sum'))
        insert_size_wide = (pd.pivot(insert_size,
                                index="sample",
                                columns="length",
                                values="count")
                       .reset_index()
        )
        insert_size_wide.to_csv(output['insertsizes'], sep="\t", index=False)

################################################################################
