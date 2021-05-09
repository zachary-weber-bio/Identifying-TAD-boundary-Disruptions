# File Name: Snakefile
# Created By: ZW
# Created On: 2021-05-04
# description- runs MANTA structural variant calling on WGS from 1000 genomes samples
# and looks for intersections with  approximate TAD boundaries

# module imports
import glob
import os.path
import pandas as pd

# read in analysis configuration details
configfile: "configs/config.yaml"
workdir: config["working_dir"]
tmpdir = config["tmp"]
anaconda_env = config["conda"]
ref_fasta = config["reference_fasta"]
call_regions = config["call_regions"]
samples = pd.read_csv(config["sample_file"],sep="\t")["SampleID"].tolist()

#file format strings
seqinfo = "alt_bwamem_GRCh38DH.20150718"
popinfo = "YRI"

# define pipeline output files
rule all:
    input:
        expand("data/{sample}.{seqinfo}.{popinfo}.low_coverage.cram.crai",
                sample=samples,seqinfo=seqinfo, popinfo=popinfo)

# 1.) index our alignment files with samtools
rule samtools_index_cram:
    input:  "data/{sample}.{seqinfo}.{popinfo}.low_coverage.cram"
    output: "data/{sample}.{seqinfo}.{popinfo}.low_coverage.cram.crai"
    conda: anaconda_env["samtools"]
    shell:
        "samtools index {input}"

# 2.) run Manta SV calling on each sample to produce a
# VCF of structural variants found in WGS
rule manta_sv_calling:
    input:
        sample="{sample}"
        ref_genome="{ref_fasts}"
        call_regions="{callregions}"
        bam="data/{sample}.{seqinfo}.{popinfo}.low_coverage.cram"
        bai="data/{sample}.{seqinfo}.{popinfo}.low_coverage.cram.crai"
        tmpdir="{tmpdir}"
        outdir="results/"
        jobs=config["rule_params"]["manta"]["jobs"]
        memGB=config["rule_params"]["manta"]["memGB"]
    output:
        "results/{sample}_diploidSV.vcf.gz"
        "results/{sample}_diploidSV.vcf.gz.tbi"
    conda:
        anaconda_env["manta"]
    shell:
        "./scripts/runManta.sh -b {input.bam} -i {input.bai}"
            " -c {input.call_regions} -r {input.ref_genome}"
            " -s {input.sample} -t {input.tmpdir} -o {input.outdir}"
            " -j {input.jobs} -g {input.memGB}"
