#!/bin/bash

CODEDIR=/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment

snakemake -s "$CODEDIR/workflow/eQTLPredictions.snakefile" --configfile "$CODEDIR/config/config.yml" --directory $CODEDIR/ --unlock

snakemake -s "$CODEDIR/workflow/eQTLPredictions.snakefile" --configfile  "$CODEDIR/config/config.yml" -j1 --use-conda --keep-target-files --rerun-incomplete

