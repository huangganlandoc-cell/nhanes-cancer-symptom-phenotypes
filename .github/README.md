# Symptom phenotypes and mortality in US cancer survivors

Analysis code for the manuscript *Symptom phenotypes and mortality in US cancer survivors: a
cohort study* (NHANES 2005–2018 with NCHS linked mortality follow-up through 2019).

Archived at Zenodo: https://doi.org/10.5281/zenodo.22970535 (this DOI always resolves to the newest version; each version has its own DOI on the Zenodo page). The code is released under the MIT License.

The repository contains code only. It contains no NHANES data: the survey files and the public-use
linked mortality files are downloaded from the CDC by the first stage of the pipeline.

- NHANES data: https://wwwn.cdc.gov/nchs/nhanes/
- NCHS public-use linked mortality files: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

## Requirements

- R 4.6 with `survey`, `survival`, `poLCA`, `rms`, `nnet` and `parallel`.
- Python 3.11+ with `pandas`, `numpy`, `statsmodels`, `matplotlib` and `requests` (for the download). Figures 1, 2 and S1–S3 are byte-identical
  to the submitted PNGs only with matplotlib 3.11.0; other versions give visually equivalent files.

## Reproducing the analysis

Run from the repository root (all scripts are in `code/`):

```bash
zsh code/run_pipeline.sh
```

This runs every stage in dependency order, from downloading the public files to the tables and
figures. A subset of stages can be named, for example `zsh code/run_pipeline.sh tables figures`.

| Stage | Scripts | What it does |
|---|---|---|
| `download` | `download_nhanes.py` | Retrieves the NHANES files for 2005–2018 and the linked mortality file |
| `cohort` | `build_cohort.py`, `build_derived_inputs.py` | Pools the seven cycles, links mortality, derives the symptom items and covariates, and writes the design frame and the file of all linkage-eligible survivors; then the all-adult measurement-model input, the dimension scores and the raw PHQ-9 item scores |
| `lca` | `reconstructed/fit_indices.R`, `analysis_bvr.R`, `reconstructed/derive_model4_covariates.R` | Latent class models with 1 to 5 classes (30 random starts), bivariate residuals, and the Model 4 covariates from the laboratory, medication, kidney and functioning files |
| `core` | `analysis_cancer_symptom_mortality.R`, `analysis_primary_exclNMS.R` | Full cohort and primary cohort: Models 1 to 3, subgroups, splines, landmark exclusion, cause-specific models |
| `dimensions` | `reconstructed/dimensional.R`, `reconstructed/all_adult_lca.R` | Dimension-score models and factor analysis; the measurement model refitted in all adults (about 30 minutes, skipped if its fits exist) |
| `threestep` | `analysis_three_step.R` (`fit`, `rep`, `boot`, `sim`, `rep3only`), `export_three_step_checks.R` | Maximum-likelihood three-step correction for misclassification, with jackknife (214 replicates) and Rao–Wu bootstrap (500 replicates) standard errors that repeat all three steps, and the simulation with a known hazard ratio |
| `model4a` | `analysis_model4a.R` | Model 3 plus the condition count and eGFR, modal and corrected |
| `sensitivity` | `analysis_inclusive_draws_score.R`, `analysis_sleep_model3.R`, `analysis_drop_slq050.R`, `analysis_fixed_score_contrast.R`, `analysis_excluded_comparison.R` | Inclusive draws, sleep splines, the measurement model without SLQ050, the phenotype contrast at a fixed PHQ-9 score, and excluded versus included survivors |
| `contrasts` | `analysis_reviewer_response.R` | Direct phenotype contrast, full landmark series, time-since-diagnosis interaction, survey-cycle analyses |
| `reconstructed` | `reconstructed/pseudoclass.R`, `reconstructed/prior_classification.R`, `reconstructed/measurement_sensitivity_RR4.R`, `reconstructed/harmonised/run_all_harmonised.R` (runs the scripts in `harmonised/`) | Pseudo-class draws, comparison with a previously published classification, measurement sensitivity analyses, Model 4, dimension scores, proportional hazards and incremental value |
| `tables` | `build_table2.py`, `build_table4.py`, `build_tables_S3_S5_S7.R`, `build_tables_S11_S12.R`, `build_tableS14.py`, `build_other_tables.py`, `fix_display_labels.py` | Every table, formatted from the saved outputs |
| `figures` | `make_figures.py`, `make_figure3.py` | Every figure, drawn from the saved outputs |

The replicate stages are the slow ones. With 16 cores the jackknife takes a few minutes and the
bootstrap about half an hour to an hour; the whole pipeline takes several hours. Set `NCORES` to
control the number of worker processes; the pipeline keeps each process to one BLAS thread.

## Notes

- Reproducibility (checked on 25 September 2026): `zsh code/run_pipeline.sh` in a fresh clone, with a
  fresh download, ran every stage in 1 h 46 min on 18 cores (macOS, R 4.6.1, Python 3.14). Every table
  (Tables 1-4, S1-S15) and every analysis output behind them, including the jackknife and bootstrap
  replicates of the three-step correction, was byte-identical to the files used for the manuscript, and
  every figure was byte-identical when drawn with matplotlib 3.11.0. The exceptions are at the level of
  floating-point rounding: the somatic residual in the dimension scores (identical with statsmodels 0.14.6
  and numpy 2.4.6, up to 4e-12 relative with statsmodels 0.15.0) and the standard errors of the latent
  class item probabilities, which no script uses. The Model 4 design frame differs from the original file
  only in years since diagnosis, corrected in `build_cohort.py`; the analyses take that variable from the
  cohort frame.

- `code/reconstructed/` holds scripts rebuilt from the interactive session in which those analyses
  were first run, each checked against the saved outputs; `compare_outputs.R` and
  `harmonised/compare_harmonised.R` perform that comparison.
- Not run by the pipeline: `reconstructed/model4.R`, `reconstructed/incremental_value.R`,
  `reconstructed/proportional_hazards.R` and `reconstructed/run_all.R` are superseded by the versions in
  `reconstructed/harmonised/`; `analysis_class_assignment.R` and `analysis_inclusive_draws.R` are earlier
  analyses not used for the manuscript (the second could not be completed; its replacement is
  `analysis_inclusive_draws_score.R`); `build_tableS10.py` only points to `build_table4.py`;
  `reconstructed/compare_outputs.R` and `harmonised/compare_harmonised.R` compare with earlier saved
  outputs and need files that the pipeline does not write.
- `code/three_step_functions.R` holds the maximum-likelihood three-step estimator: the
  classification-error matrix, the EM algorithm with a Breslow baseline hazard, the bivariate
  residuals and Rubin's rules.
- Years since diagnosis: the interactive build first read the age at diagnosis for thyroid, uterine
  and other cancers from columns that do not exist in the 2005–2016 files. `build_cohort.py` reads
  the correct columns, and every analysis was rerun on the corrected data (2026-09-25).
- `build_docx.py`, `build_journal_variant.py`, `build_references.py`, `fill_strobe_pages.py` and
  `sync_bilingual.py` assemble the submission files from the manuscript source, which is not part
  of this repository.
