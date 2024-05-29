# Genetic-risk-lifestyles-and-early-onset-CVD
Code for paper entitled "Joint impact of polygenic risk score and lifestyles on early- and late-onset cardiovascular diseases".

1-GWAS_cleaning.R was used to clean the GWAS summary data.
2-1-GWASQC.sh/2-2-GWASQC.pbs/2-3-GWASQC.R were used for QC of the GWAS summary data.
3-1-PRS_CT.sh/3-2-PRS_LS.sh/3-1-PRS_CS.sh/3-4-PRS-catalog.sh were main code to build PRSs using C+T, lassosum, and PRS-CS and derive previous PRSs released in PGS catalog.
4-MetaPRS.R was used to build the MetaPRS using EN.
5-Compare_PRS.R was used to compare asssociations between MetaPRS and trait-specific PRS.
6-1-Age_at_onset.sh/6-2-Age_at_onset.pbs/6-3-Age_at_onset.R were used to estimated the associations between genetic risk and lifestyles with different age-onsets of CVD.
7-Additive_boot.R was used to calculated addtive interaction measures.
8-1-IR_adj_boot.sh/8-2-IR_adj_boot.pbs/8-3-IR_adj_boot.R were used to estimated the incidence rates.
The other codes are branch files for building PRSs.
