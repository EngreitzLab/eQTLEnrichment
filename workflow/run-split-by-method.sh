#!/bin/bash

#quick-sub -m 250G -t 24:00:00 -s dist_to_gene.qsh -o dist_to_gene.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-dist_to_gene.sh"

#quick-sub -m 250G -t 24:00:00 -s dist_to_tss.qsh -o dist_to_tss.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-dist_to_tss.sh"

#quick-sub -m 250G -t 24:00:00 -s nearest_gene.qsh -o nearest_gene.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-nearest_gene.sh"

#quick-sub -m 250G -t 24:00:00 -s nearest_tss.qsh -o nearest_tss.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-nearest_tss.sh"

#quick-sub -m 250G -t 24:00:00 -s within_100kb.qsh -o within_100kb.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-within_100kb_of_tss.sh"

quick-sub -m 250G -t 24:00:00 -s reads_by_dist.qsh -o reads_by_dist.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-reads_by_dist_to_tss.sh"

#quick-sub -m 250G -t 24:00:00 -s norm.qsh -o norm.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-reads_by_dist_to_tss_norm.sh"

quick-sub -m 250G -t 24:00:00 -s H3K_reads_by_dist.qsh -o H3K_reads_by_dist.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-H3K27ac_reads_by_dist_to_tss.sh"

#quick-sub -m 250G -t 24:00:00 -s H3K_norm.qsh -o H3K_norm.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-H3K27ac_reads_by_dist_to_tss_norm.sh"

quick-sub -m 250G -t 24:00:00 -s EpiMap.qsh -o EpiMap.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-EpiMap.sh"

quick-sub -m 250G -t 24:00:00 -s ABC.qsh -o ABC.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-ABC.sh"

quick-sub -m 250G -t 24:00:00 -s GraphRegCl.qsh -o GraphRegCl.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-GraphRegCl.sh"

quick-sub -m 250G -t 24:00:00 -s Ramil.qsh -o Ramil.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-Ramil.sh"

quick-sub -m 250G -t 24:00:00 -s bp.qsh -o bp.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-baseline_predictors.sh"

quick-sub -m 250G -t 24:00:00 -s v2.qsh -o v2.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-LogRegClassifier_v2.sh"

quick-sub -m 250G -t 24:00:00 -s new.qsh -o new.qout "source ~/.bashrc; conda activate eQTLEnv; sh run-new_methods.sh"
