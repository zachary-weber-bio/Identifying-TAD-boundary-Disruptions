#!/bin/bash
# File Name: runManta.sh
# Created By: ZW
# Created On: 2021-05-09
# Purpose: runs the Illumina MANTA SV caller. on a WGS sample

# Usage message for the runMANTA.sh script
usage () {
    echo "Usage: runManta.sh [-b <BAM file>] [-i <BAI file>]
                    [-s <Sample ID>] [-r <Ref. FASTA>]
                    [-t <Temp Directory>] [-o <Output Directory>]
                    [-j <Jobs>] [-g <Memory(GB)>]" 1>&2;
    exit 1;
}

# getopt setup for named arguments
while getopts ":b:i:s:r:c:t:o:j:g:" x; do
    case "${x}" in
        b) b=${OPTARG} ;;
        i) i=${OPTARG} ;;
	s) s=${OPTARG} ;;
        r) r=${OPTARG} ;;
	c) c=${OPTARG} ;;
        t) t=${OPTARG} ;;
        o) o=${OPTARG} ;;
        j) j=${OPTARG} ;;
        g) g=${OPTARG} ;;
        *) usage ;;
    esac
done

# ----- Run MANTA ------
# switch to temp execution directory and
# run the MANTA configuration step
projectdir=`pwd`/
mkdir ${t}${s}/
cd ${t}${s}/

configManta.py \
    --bam ${b} \
    --referenceFasta ${r} \
    --callRegions ${c}

# run MANTA, then collect and rename results
cd MantaWorkflow/
./runWorkflow.py -j ${j} -g ${g}
mv results/variants/diploidSV.vcf.gz results/variants/${s}_diploidSV.vcf.gz
mv results/variants/diploidSV.vcf.gz.tbi results/variants/${s}_diploidSV.vcf.gz.tbi

cp results/variants/${s}_diploidSV.vcf.gz ${o}
cp results/variants/${s}_diploidSV.vcf.gz.tbi ${o}
