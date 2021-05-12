# File Name: Snakefile
# Created By: ZW
# Created On: 2021-05-04
# description- runs MANTA structural variant calling on WGS from 1000 genomes samples
# and looks for intersections with  approximate TAD boundaries

# module imports
import os.path
import pandas as pd

# read in analysis configuration details
configfile: "configs/config.yaml"
workdir: config["working_dir"]
tmpdir = config["tmp"]
datadir = config["data_dir"]
resultsdir = config["results_dir"] 
anaconda_env = config["conda"]
samples = pd.read_csv(config["sample_file"],sep="\t")["SampleID"].tolist()

#-------------------
#---------
#--
#-
#

# define pipeline output files
rule all:
    input:
        expand(datadir + "{sample}.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram.crai",sample=samples),
        expand(resultsdir + "{sample}_diploidSV.vcf.gz",sample=samples),
        expand(resultsdir + "{sample}_diploidSV.vcf.gz.tbi",sample=samples),
	resultsdir + "SV_TAD-boundary_intersects.sorted.bed",
	expand(resultsdir + "evidence/{sample}.SVevidence.bam",sample=samples)

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
rule manta_sv_calling:
    input:
        ref_genome = config["reference_fasta"],
	call_regions = config["call_regions"],
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

# 3.) find and enumerate structural variants which are thought to intersect TAD
# boundary regions in one or more samples
rule tad_boundary_intersect:
    input:
        tad_boundaries = config["tad_boundaries"],
        sv_calls_vcfs = expand(resultsdir + "{sample}_diploidSV.vcf.gz",sample=samples)
    params: samples = expand("{sample}",sample=samples)
    output: sorted_intersects = resultsdir + "SV_TAD-boundary_intersects.sorted.bed"
    conda: anaconda_env["samtools"]
    shell:
        "bedtools intersect -wo -a {input.tad_boundaries} -b {input.sv_calls_vcfs}"
        " -names {params.samples} | sort -k1,1V -k2,2n > {output.sorted_intersects}"

# 4.) generate "evidence bam" containing reads in the areas where called
# structural variants have been found.
rule generate_evidence_bam:
    input:
        vcf = resultsdir + "{sample}_diploidSV.vcf.gz",
        cram = datadir + "{sample}.alt_bwamem_GRCh38DH.20150718.YRI.low_coverage.cram",
        genome_file = config["genome_file"]
    params:
        sample = "{sample}",
        margin = str(config["rule_params"]["bedtools"]["evidence_roi_margin"])
    output: resultsdir + "evidence/{sample}.SVevidence.bam"
    conda: anaconda_env["samtools"]
    shell:
        "./scripts/generate_evidence_BAM.sh -s {params.sample}"
        " -v {input.vcf} -b {input.cram} -m {params.margin}"
        " -g {input.genome_file} -o {output}"








	
