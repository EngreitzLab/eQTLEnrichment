#!/bin/bash

CODEDIR=/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-integrated/eQTLEnrichment

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-GM12878_predictors.yml" --directory $CODEDIR/workflow --unlock

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-GM12878_predictors.yml" -j6 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R first

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-GM12878_predictors.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R second

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-GM12878_predictors.yml" -j1 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R third

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-GM12878_predictors.yml" -j1 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R fourth



