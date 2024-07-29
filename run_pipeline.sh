#!/bin/bash

snakemake \
-n \
--snakefile Snakefile_bam_merge.smk \
--configfile settings_rna/config_rna.yaml \
--use-conda \
--conda-frontend conda \
--printshellcmds \
--rerun-incomplete \
-j 50 \
--cluster-config cluster.yaml \
--cluster "sbatch -c {cluster.cpus} \
-t {cluster.time} -N {cluster.nodes} \
--mem {cluster.mem} \
--output {cluster.output} \
--error {cluster.error}" \
--cluster-cancel scancel \
"$@"
