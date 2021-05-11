#!/bin/bash

module load python/anaconda-2020.02
source activate snakemake

snakemake \
	--snakefile Snakefile \
	--use-conda \
	-k -p -j 500 \
	--rerun-incomplete \
	--cluster-config configs/cluster.json \
	-c "sbatch --mem={cluster.memory} \
		--nodes={cluster.n} \
		--time={cluster.time} \
		--tasks-per-node=2 \
		--partition=broadwl \
		--job-name={cluster.name} \
		--output={cluster.output}
		--error={cluster.error}" \
	$*
