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
datadir = config["data_dir"]
resultsdir = config["results_dir"] 
anaconda_env = config["conda"]
ref_fasta = config["reference_fasta"]
call_regions = config["call_regions"]
samples = pd.read_csv(config["sample_file"],sep="\t")["SampleID"].tolist()

# define pipeline output files
rule all:
    input:
        expand(datadir + "{sample}.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram.crai",sample=samples),
        expand(resultsdir + "{sample}_diploidSV.vcf.gz",sample=samples),
        expand(resultsdir + "{sample}_diploidSV.vcf.gz.tbi",sample=samples)

# 1.) index our alignment files with samtools
rule samtools_index:
    input: cram = datadir + "{sample}.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram",
    params: threads = str(config["rule_params"]["samtools"]["threads"])
    output: datadir + "{sample}.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram.crai"
    conda: anaconda_env["samtools"]
    shell:
        "samtools index -@ {params.threads} {input.cram}"

# 2.) run Manta SV calling on each sample to produce a
# VCF of structural variants found in WGS
print(samples)
rule manta_sv_calling:
    input:
        ref_genome = ref_fasta,
	call_regions = call_regions,
        cram = datadir + "{sample}.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram",
        crai = datadir + "{sample}.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram.crai",
    params:
        sample = "{sample}",
	tmpdir = tmpdir,
	outdir = resultsdir,
	jobs = str(config["rule_params"]["manta"]["jobs"]),
        memGB = str(config["rule_params"]["manta"]["memGB"])
    output:
        resultsdir + "{sample}_diploidSV.vcf.gz",
        resultsdir + "{sample}_diploidSV.vcf.gz.tbi"
    conda:
        anaconda_env["manta"]
    shell:
        "./scripts/runManta.sh -b {input.cram} -i {input.crai}"
        " -s {params.sample} -r {input.ref_genome} -c {input.call_regions}"
        " -s {params.sample} -t {params.tmpdir} -o {params.outdir}"
        " -j {params.jobs} -g {params.memGB}"
