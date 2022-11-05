#!/bin/bash

CODEDIR=/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-nearest_gene.yml" --directory $CODEDIR/workflow --unlock

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-nearest_gene.yml" -j3 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R first

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-nearest_gene.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R second

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-nearest_gene.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R third

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-nearest_gene.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R fourth

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-nearest_gene.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R fifth