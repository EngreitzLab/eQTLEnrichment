#!/bin/bash

CODEDIR=/oak/stanford/groups/engreitz/Users/sheth/eQTLEnrichment-thresholding/eQTLEnrichment

#snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-Ramil.yml" --directory $CODEDIR/workflow --unlock

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-Enformer.yml" -j6 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R first

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-Enformer.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R second

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-Enformer.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R third

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-Enformer.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R fourth

snakemake -s "$CODEDIR/workflow/Snakefile.snakefile" --configfile "$CODEDIR/config/config-hg38-v2-Enformer.yml" -j4 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "$CODEDIR/workflow/envs/" -R fifth


#snakemake -s "workflow/Snakefile.snakefile" --configfile "config/config-hg38-v2-Ramil.yml" -j1 --use-conda --keep-target-files --rerun-incomplete --conda-frontend conda -k --conda-prefix "workflow/envs/" -R second
