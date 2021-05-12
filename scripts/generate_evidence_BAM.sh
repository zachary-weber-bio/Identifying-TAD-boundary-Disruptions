#!/bin/bash
# File Name: generate_evidence_BAM.sh
# Created By: ZW
# Created On: 2021-05-11
# Purpose: uses bcftools to extract VCF record info from SV calls into
# into bed format, then uses bedtools to add "slop" around the locus, before
# finally using samtools generate a small BAM file of the reads around the variant

# Usage message for the generate_evidence_BAM.sh script
usage () {
    echo "Usage: generate_evidence_BAM.sh [-s <SAMPLE ID>]
            [-v <VCF file>][-b <CRAM/BAM to subset>]
            [-m <MARGIN (bp to extend interval by)>]
            [-g <GENOME file (as described by bedtools)>]
            [-o <OUTPUT bam>]" 1>&2
    exit 1;
}

# getopt setup for named arguments
while getopts ":s:v:b:m:g:o:" x; do
    case "${x}" in
        s) s=${OPTARG} ;;
        v) v=${OPTARG} ;;
        b) b=${OPTARG} ;;
        m) m=${OPTARG} ;;
        g) g=${OPTARG} ;;
        o) o=${OPTARG} ;;
        *) usage ;;
    esac
done

# pull regions of interest from manta VCF file
bcftools query -f "%CHROM\t%POS\t%END\n" ${v} | \
    bedtools slop -i stdin -g ${g} -b ${m} > ${s}_roi.bed

samtools view -L ${s}_roi.bed -b ${b} > ${o}
rm ${s}_roi.bed
