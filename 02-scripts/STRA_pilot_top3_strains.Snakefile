################################################################################
# Project: LIFE Child cohort - plaque samples
# Part: Strain analysis
# Step: Strain analysis of the 3 most abundant strains in the pilot samples
#
# Dependent on:
#   - COMP_taxprofiling.Snakefile
#
# Alex Huebner, 13/02/25
################################################################################

import os
from collections import Counter

import pandas as pd

if not os.path.isdir("snakemake_tmp"):
    os.makedirs("snakemake_tmp")

os.environ["OPENBLAS_NUM_THREADS"] = '1'
os.environ["OMP_NUM_THREADS"] = '1'

#### SAMPLES ###################################################################
SAMPLES = ['LLP001.A0101', 'LLP002.A0101', 'LLP003.A0101']
# Identify top species
# metaphlan4 = pd.read_csv("05-results/COMP_TaxProfiler_MetaPhlAnv4_Jun23_pilot.tsv",
#                          sep="\t",
#                          usecols=['clade_name', 'LLP001.A0101', 'LLP002.A0101', 'LLP003.A0101'])
# median_relab = (pd.melt(metaphlan4, id_vars='clade_name', var_name="sample",
#                        value_name="relAb")
#     .groupby(['clade_name'], as_index=False)['relAb'].median()
#     .sort_values(['relAb'], ascending=[False])
# )
# top_species = median_relab.iloc[:3,:]['clade_name'].values
SPECIES = {
    'Nsicca': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/753/665/GCF_017753665.1_ASM1775366v1/GCF_017753665.1_ASM1775366v1_genomic.fna.gz',
    'Lmirabilis': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/637/555/GCF_900637555.1_51726_E01/GCF_900637555.1_51726_E01_genomic.fna.gz',
    'Soralis': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/637/025/GCF_900637025.1_46338_H01/GCF_900637025.1_46338_H01_genomic.fna.gz',
}
################################################################################

rule all:
    input:
        expand("03-data/refgenomes/{species}.fna.gz", species=SPECIES),
        expand("tmp/strain_analysis_pilot/{sample}.{species}.calmd.sorted.bam.bai", sample=SAMPLES, species=SPECIES),
        expand("04-analysis/strain_analysis/pilot_samples.{species}.vcf", species=SPECIES),
        "05-results/STRA_pilot_strains_sitefreqspectrum.tsv"

#### Prepare the referene genomes ##############################################

rule download_genome:
    output:
        "03-data/refgenomes/{species}.fna.gz"
    message: "Download the reference genome: {wildcards.species}"
    container: "https://depot.galaxyproject.org/singularity/samtools%3A1.21--h96c455f_1"
    resources:
        mem = 2,
        cores = 1
    params:
        url = lambda wildcards: SPECIES[wildcards.species],
        fasta = "03-data/refgenomes/{species}.fna"
    threads: 1
    shell:
        """
        wget -O {output} {params.url}
        gunzip {output}
        bgzip {params.fasta}
        """

################################################################################

#### Align the sequencing data against the genomes #############################

rule build_BowTie2_index:
    input:
        "03-data/refgenomes/{species}.fna.gz"
    output:
        bt1 = temp("tmp/strain_analysis_pilot/refgenomes/{species}.1.bt2"),
        bt2 = temp("tmp/strain_analysis_pilot/refgenomes/{species}.2.bt2"),
        bt3 = temp("tmp/strain_analysis_pilot/refgenomes/{species}.3.bt2"),
        bt4 = temp("tmp/strain_analysis_pilot/refgenomes/{species}.4.bt2"),
        revbt1 = temp("tmp/strain_analysis_pilot/refgenomes/{species}.rev.1.bt2"),
        revbt2 = temp("tmp/strain_analysis_pilot/refgenomes/{species}.rev.2.bt2")
    message: "Index the contigs for alignment using BowTie2: {wildcards.species}"
    container: "docker://quay.io/biocontainers/bowtie2:2.5.1--py310ha0a81b8_2"
    resources:
        mem = 8,
        cores = 4
    params:
        index = "tmp/strain_analysis_pilot/refgenomes/{species}"
    threads: 4
    shell:
        """
        bowtie2-build -f --threads {threads} \
            {input} \
            {params.index}
        """

rule BowTie2_alignment:
    input:
        bt1 = "tmp/strain_analysis_pilot/refgenomes/{species}.1.bt2",
        bt2 = "tmp/strain_analysis_pilot/refgenomes/{species}.2.bt2",
        bt3 = "tmp/strain_analysis_pilot/refgenomes/{species}.3.bt2",
        bt4 = "tmp/strain_analysis_pilot/refgenomes/{species}.4.bt2",
        revbt1 = "tmp/strain_analysis_pilot/refgenomes/{species}.rev.1.bt2",
        revbt2 = "tmp/strain_analysis_pilot/refgenomes/{species}.rev.2.bt2"
    output:
        pipe("tmp/strain_analysis_pilot/{sample}.{species}.sam")
    message: "Align the reads using BowTie2's very-sensitive setting: {wildcards.sample} against the species {wildcards.species}"
    container: "docker://quay.io/biocontainers/bowtie2:2.5.1--py310ha0a81b8_2"
    group: "ref_alignment"
    resources:
        mem = 14,
        cores = 16
    params:
        index = "tmp/strain_analysis_pilot/refgenomes/{species}",
        pe1 = "04-analysis/taxprofiler/results_pilot/run_merging/{sample}_1.merged.fastq.gz",
        pe2 = "04-analysis/taxprofiler/results_pilot/run_merging/{sample}_2.merged.fastq.gz"
    threads: 14
    shell:
        """
        bowtie2 -p {threads} --very-sensitive -x {params.index} \
                -1 {params.pe1} -2 {params.pe2} -S {output}
        """

rule samtools_sort:
    input:
        bam = "tmp/strain_analysis_pilot/{sample}.{species}.sam",
        fasta = "03-data/refgenomes/{species}.fna.gz"
    output:
        "tmp/strain_analysis_pilot/{sample}.{species}.calmd.sorted.bam",
    message: "Sort the sequencing data by coordinate and calculate the MD field: {wildcards.sample} against the species {wildcards.species}"
    container: "https://depot.galaxyproject.org/singularity/samtools%3A1.21--h96c455f_1"
    group: "ref_alignment"
    resources:
        mem = 8,
        cores = 2
    threads: 2
    shell:
        """
        samtools view -Sb {input.bam} | \
        samtools calmd -u /dev/stdin {input.fasta} | \
        samtools sort -l 4 -o {output} -
        """

rule samtools_addrg:
    input:
        "tmp/strain_analysis_pilot/{sample}.{species}.calmd.sorted.bam",
    output:
        "04-analysis/strain_analysis/{sample}.{species}.calmd.sorted.bam"
    message: "Add RG tag: {wildcards.sample} against the species {wildcards.species}"
    container: "https://depot.galaxyproject.org/singularity/samtools%3A1.21--h96c455f_1"
    resources:
        mem = 4,
        cores = 1
    params:
        rgtag = lambda wildcards: f"-r ID:{wildcards.sample} -r SM:{wildcards.sample}"
    threads: 1
    shell:
        "samtools addreplacerg {params.rgtag} -o {output} {input}"

rule samtools_index:
    input:
        "04-analysis/strain_analysis/{sample}.{species}.calmd.sorted.bam"
    output:
        "04-analysis/strain_analysis/{sample}.{species}.calmd.sorted.bam.bai"
    message: "Index the BAM file: {wildcards.sample} against the species {wildcards.species}"
    container: "https://depot.galaxyproject.org/singularity/samtools%3A1.21--h96c455f_1"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    shell:
        """
        samtools index {input}
        """

################################################################################

#### Infer the genotypes #######################################################

rule samtools_depth:
    input:
        bam = "04-analysis/strain_analysis/{sample}.{species}.calmd.sorted.bam",
        bai = "04-analysis/strain_analysis/{sample}.{species}.calmd.sorted.bam.bai"
    output:
        "tmp/strain_analysis_pilot/{sample}.{species}.depth.tsv"
    message: "Determine the depth along the genome: {wildcards.sample} against the species {wildcards.species}"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    shell:
        "samtools depth -a {input.bam} > {output}"

rule decompress_fasta:
    input:
        "03-data/refgenomes/{species}.fna.gz"
    output:
        temp("tmp/strain_analysis_pilot/refgenomes/{species}.fna")
    message: "Decompress the FastA file: {wildcards.species}"
    resources:
        mem = 2,
        cores = 1
    threads: 1
    shell:
        "gunzip -c {input} > {output}"

rule faidx_reffasta:
    input:
        "tmp/strain_analysis_pilot/refgenomes/{species}.fna"
    output:
        temp("tmp/strain_analysis_pilot/refgenomes/{species}.fna.fai")
    message: "Generate FastA index: {wildcards.species}"
    container: "https://depot.galaxyproject.org/singularity/samtools%3A1.21--h96c455f_1"
    resources:
        mem = 2,
        cores = 1
    shell:
        """
        samtools faidx {input}
        """

rule determine_chunks_freebayes:
    input:
        fai = "tmp/strain_analysis_pilot/refgenomes/{species}.fna.fai",
        depth = "tmp/strain_analysis_pilot/{sample}.{species}.depth.tsv" 
    output:
        temp("tmp/strain_analysis_pilot/{sample}.{species}.chunks")
    message: "Determine regions with approx. equal coverage to have 100 chuncks: {wildcards.sample}"
    container: "docker://quay.io/biocontainers/freebayes:1.3.7--h1870644_0"
    resources:
        mem = lambda wildcards, attempt: 16 + attempt * 16,
        cores = 1
    threads: 1
    shell:
        """
        cat {input.depth} | \
        coverage_to_regions.py {input.fai} 100 > {output}
        """

rule freebayes:
    input:
        bam = lambda wildcards: [f"04-analysis/strain_analysis/{sample}.{wildcards.species}.calmd.sorted.bam" for sample in SAMPLES],
        bai = lambda wildcards: [f"04-analysis/strain_analysis/{sample}.{wildcards.species}.calmd.sorted.bam.bai" for sample in SAMPLES],
        fasta = "tmp/strain_analysis_pilot/refgenomes/{species}.fna",
        fai = "tmp/strain_analysis_pilot/refgenomes/{species}.fna.fai",
        chunks = lambda wildcards: [f"tmp/strain_analysis_pilot/{sample}.{wildcards.species}.chunks" for sample in SAMPLES]
    output:
        "04-analysis/strain_analysis/pilot_samples.{species}.vcf"
    message: "Call genotypes across the pilot samples: {wildcards.species}"
    conda: "ENVS_freebayes.yaml"
    resources:
        mem = 16,
        cores = 8
    params:
        chunks = lambda wildcards: [f"tmp/strain_analysis_pilot/{sample}.{wildcards.species}.chunks" for sample in SAMPLES[:1]]
    threads: 8
    shell:
        """
        freebayes-parallel {params.chunks} \
                {threads} -f {input.fasta} -p 1 -q 30 {input.bam} > {output}
        """

################################################################################

#### Filter VCF file ###########################################################

rule filter_vcf:
    input:
        "04-analysis/strain_analysis/pilot_samples.{species}.vcf"
    output:
        "tmp/strain_analysis_pilot/pilot_samples.{species}.biallelic.vcf"
    message: "Filter for high-quality bi-allelic sites: {wildcards.species}"
    container: "https://depot.galaxyproject.org/singularity/bcftools%3A1.21--h3a4d415_1"
    resources:
        mem = 4,
        cores = 1
    threads: 1
    shell:
        """
        bcftools view -v snps -i 'QUAL >= 30 && INFO/AN == 3' {input} | \
        bcftools annotate -x ^FORMAT/GT,INFO - > {output}
        """

rule sitefreq_spectrum:
    input:
        "tmp/strain_analysis_pilot/pilot_samples.{species}.biallelic.vcf"
    output:
        "04-analysis/strain_analysis/pilot_samples.{species}.sitespectrum.txt"
    message: "Identify the bi-allelic sites: {wildcards.species}"
    run:
        site_freq = []
        with open(input[0], "rt") as vcffile:
            for line in vcffile:
                if line.startswith("#CHROM"):
                    samples = line.rstrip().split("\t")[9:]
                elif not line.startswith("#"):
                    fields = line.rstrip().split("\t")
                    if '2' not in fields[9:] and '3' not in fields[9:]:
                        site_freq.append("".join(fields[9:]))
        counts = Counter(site_freq)
        with open(output[0], "wt") as outfile:
            for p, c in counts.items():
                outfile.write(f"{p}\t{c}\n")

rule summarise_sfs:
    input:
        expand("04-analysis/strain_analysis/pilot_samples.{species}.sitespectrum.txt", species=SPECIES)
    output:
        "05-results/STRA_pilot_strains_sitefreqspectrum.tsv"
    message: "Summarise the site frequency spectra across the species"
    run:
        sfspectra = []
        for fn in input:
            sfs = (pd.read_csv(fn, sep="\t", header=None, names=['pattern', 'count'],
                               dtype={'pattern': str,
                                      'count': int})
                   .assign(filename=fn))
            sfspectra.append(sfs)

        sfs_df = pd.concat(sfspectra)
        sfs_df['species'] = sfs_df['filename'].str.extract(r'04-analysis/strain_analysis/pilot_samples.([A-Za-za]+).sitespectrum.txt')

        (sfs_df.drop(['filename'], axis=1)[['species', 'pattern', 'count']]
         .to_csv(output[0], sep="\t", index=False))

################################################################################
