#!/bin/zsh
# Regenerate every analysis output, table and figure from the public NHANES files, in dependency order.
#
#   zsh code/run_pipeline.sh              # all stages (several hours; the replicate stages dominate)
#   zsh code/run_pipeline.sh tables figures
#
# Stages: download cohort lca core threestep model4a sensitivity contrasts reconstructed tables figures
# Several scripts write to review/reproduction/output/ or review/reproduction/harmonised/; the copy
# steps below put the files that the tables, figures and manuscript use into supporting/.
# Requirements: R 4.6 (survey, survival, poLCA, rms, nnet); Python 3 with pandas, numpy, matplotlib.
# Figures 1, 2 and S1-S3 are byte-identical to the published PNGs only with matplotlib 3.11.0.
set -e
cd "$(dirname "$0")/.."
[[ -f code/run_pipeline.sh ]] || { echo "run from the repository root" >&2; exit 1; }
# one BLAS thread per process: the replicate stages already use every core through mclapply
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1} OMP_NUM_THREADS=${OMP_NUM_THREADS:-1} VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
R=${RSCRIPT:-Rscript}; PY=${PYTHON:-python3}
L=review/logs; RO=review/reproduction/output; RH=review/reproduction/harmonised
mkdir -p data/derived $L $RO/logs supporting tables supplementary figures
say() { echo "== $* ($(date +%H:%M))"; }
cp_to() { local dst=$1; shift; for f in "$@"; do cp "$f" "$dst"; done; }

download() { say download; $PY code/download_nhanes.py --out data; }

cohort() { say cohort; $PY code/build_cohort.py; }

lca() {
  say "latent class models (writes data/derived/lca_fits_cancer.rds if absent)"
  $R code/reconstructed/fit_indices.R > $RO/logs/fit_indices.log 2>&1
  cp $RO/lca_fit_indices_cancer.csv supporting/
  $R code/analysis_bvr.R > $L/bvr.log 2>&1
  $R code/reconstructed/derive_model4_covariates.R > $RO/logs/derive_model4_covariates.log 2>&1
  $R code/reconstructed/dimensional.R > $RO/logs/dimensional.log 2>&1
  cp $RO/RR_efa_loadings.csv supporting/
  if [[ ! -f $RO/intermediate/lca_fits_all_adults.rds ]]; then     # about 30 minutes; measurement only
    say "latent class model in all adults"; $R code/reconstructed/all_adult_lca.R > $RO/logs/all_adult_lca.log 2>&1
    cp_to supporting/ $RO/lca_fit_all_adults.csv $RO/RR_fullsample_lca_profiles.csv
  fi
}

core() {
  say "full cohort and primary cohort"
  OUT=$RO/full_cohort $R code/analysis_cancer_symptom_mortality.R > $L/full_cohort.log 2>&1
  cp $RO/full_cohort/analysis_frame.rds data/derived/
  cp $RO/full_cohort/lca_profiles_cancer.csv supporting/
  $R code/analysis_primary_exclNMS.R > $L/primary.log 2>&1
}

threestep() {
  for st in fit rep boot sim rep3only; do
    say "three-step $st"; $R code/analysis_three_step.R $st > $L/three_step_$st.log 2>&1
  done
  $R code/export_three_step_checks.R > $L/export_three_step_checks.log 2>&1
}

model4a() {
  for st in fit rep boot; do say "Model 4a $st"; $R code/analysis_model4a.R $st > $L/model4a_$st.log 2>&1; done
}

sensitivity() {
  say "inclusive draws"; $R code/analysis_inclusive_draws_score.R > $L/inclusive_draws_score.log 2>&1
  say "sleep splines"; $R code/analysis_sleep_model3.R > $L/sleep_model3.log 2>&1
  say "measurement model without SLQ050"; $R code/analysis_drop_slq050.R > $L/drop_slq050.log 2>&1
  say "phenotype contrast at a fixed PHQ-9 score"; $R code/analysis_fixed_score_contrast.R > $L/fixed_score_contrast.log 2>&1
  say "excluded versus included survivors"; $R code/analysis_excluded_comparison.R > $L/excluded_comparison.log 2>&1
}

contrasts() {
  say "contrasts, landmark, time since diagnosis, cycle"
  OUT=$RO/reviewer_response $R code/analysis_reviewer_response.R > $L/reviewer_response.log 2>&1
  cp_to supporting/ $RO/reviewer_response/{RR_events_per_phenotype,RR_contrast_vs_insomnia_fatigue,RR_interaction_time_since_dx,RR_landmark_full,RR_cycle_sensitivity}.csv
}

reconstructed() {
  say "pseudo-class draws"; $R code/reconstructed/pseudoclass.R > $RO/logs/pseudoclass.log 2>&1
  cp_to supporting/ $RO/sens_pseudoclass_M100.csv $RO/sens_pseudoclass.csv
  say "prior classification"; $R code/reconstructed/prior_classification.R > $RO/logs/prior_classification.log 2>&1
  cp_to supporting/ $RO/RR_vs_prior_classification.csv $RO/RR_crosstab_vs_prior.csv $RO/RR_prior_joint_tests.csv
  say "measurement sensitivity"; $R code/reconstructed/measurement_sensitivity_RR4.R > $RO/logs/measurement_sensitivity_RR4.log 2>&1
  cp_to supporting/ $RO/RR4_{localdep_cox,threshold_cox,measurement_fit,inclusive_pseudoclass,inclusive_profiles}.csv
  say "Model 4, dimensional scores, all-adult draws, proportional hazards, incremental value"
  $R code/reconstructed/harmonised/run_all_harmonised.R > $L/harmonised.log 2>&1
  cp_to supporting/ $RH/{RR_model4,RR_model4_contrast,RR_model4_pseudoclass,RR_model4_multiplicity,RR_model4_evalues,RR_contrast_pseudoclass,RR_attenuation_comparison,RR_power_confirmatory,RR_dimensional,RR_dimensional_affective,RR_dimensional_multiplicity,RR_fullsample_lca_pseudoclass,ph_tests,ph_period_specific,ph_stratified_baseline,RR_incremental_value,RR_incremental_joint}.csv
  cp $RH/model4_cause_specific.csv supporting/HM_model4_cause_specific.csv
  for f in incremental identifiability cycle; do cp $RH/diagnostics_$f.csv supporting/HM_diagnostics_$f.csv; done
}

tables() {
  say tables
  $PY code/build_table2.py; $PY code/build_table4.py
  $R code/build_tables_S3_S5_S7.R; $R code/build_tables_S11_S12.R
  $PY code/build_tableS14.py; $PY code/build_other_tables.py; $PY code/fix_display_labels.py
}

figures() { say figures; $PY code/make_figures.py; $PY code/make_figure3.py; }

stages=(${@:-download cohort lca core threestep model4a sensitivity contrasts reconstructed tables figures})
for s in $stages; do $s; done
say done
