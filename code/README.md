# Symptom phenotypes and mortality in US cancer survivors

Analysis code for the manuscript *Symptom phenotypes and mortality in US cancer survivors: a
cohort study* (NHANES 2005–2018 with NCHS linked mortality follow-up through 2019).

The repository contains code only. It contains no NHANES data: the survey files and the public-use
linked mortality files are downloaded from the CDC by the first stage of the pipeline.

- NHANES data: https://wwwn.cdc.gov/nchs/nhanes/
- NCHS public-use linked mortality files: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm

## Requirements

- R 4.6 with `survey`, `survival`, `poLCA`, `rms`, `nnet` and `parallel`.
- Python 3.11+ with `pandas`, `numpy` and `matplotlib`. Figures 1, 2 and S1–S3 are byte-identical
  to the submitted PNGs only with matplotlib 3.11.0; other versions give visually equivalent files.

## Reproducing the analysis

Run from the repository root:

```bash
zsh code/run_pipeline.sh
```

This runs every stage in dependency order, from downloading the public files to the tables and
figures. A subset of stages can be named, for example `zsh code/run_pipeline.sh tables figures`.

| Stage | Scripts | What it does |
|---|---|---|
| `download` | `download_nhanes.py` | Retrieves the NHANES files for 2005–2018 and the linked mortality file |
| `cohort` | `build_cohort.py` | Pools the seven cycles, links mortality, derives the symptom items and covariates, and writes the design frame and the file of all linkage-eligible survivors |
| `lca` | `reconstructed/fit_indices.R`, `analysis_bvr.R`, `reconstructed/derive_model4_covariates.R`, `reconstructed/dimensional.R`, `reconstructed/all_adult_lca.R` | Latent class models with 1 to 5 classes (30 random starts), bivariate residuals, the Model 4 covariates, dimension scores, and the measurement model refitted in all adults |
| `core` | `analysis_cancer_symptom_mortality.R`, `analysis_primary_exclNMS.R` | Full cohort and primary cohort: Models 1 to 3, subgroups, splines, landmark exclusion, cause-specific models |
| `threestep` | `analysis_three_step.R` (`fit`, `rep`, `boot`, `sim`, `rep3only`), `export_three_step_checks.R` | Maximum-likelihood three-step correction for misclassification, with jackknife (214 replicates) and Rao–Wu bootstrap (500 replicates) standard errors that repeat all three steps, and the simulation with a known hazard ratio |
| `model4a` | `analysis_model4a.R` | Model 3 plus the condition count and eGFR, modal and corrected |
| `sensitivity` | `analysis_inclusive_draws_score.R`, `analysis_sleep_model3.R`, `analysis_drop_slq050.R`, `analysis_fixed_score_contrast.R`, `analysis_excluded_comparison.R` | Inclusive draws, sleep splines, the measurement model without SLQ050, the phenotype contrast at a fixed PHQ-9 score, and excluded versus included survivors |
| `contrasts` | `analysis_reviewer_response.R` | Direct phenotype contrast, full landmark series, time-since-diagnosis interaction, survey-cycle analyses |
| `reconstructed` | `reconstructed/*.R`, `reconstructed/harmonised/*.R` | Pseudo-class draws, comparison with a previously published classification, measurement sensitivity analyses, Model 4, dimension scores, proportional hazards and incremental value |
| `tables` | `build_table2.py`, `build_table4.py`, `build_tables_S3_S5_S7.R`, `build_tables_S11_S12.R`, `build_tableS14.py`, `build_other_tables.py`, `fix_display_labels.py` | Every table, formatted from the saved outputs |
| `figures` | `make_figures.py`, `make_figure3.py` | Every figure, drawn from the saved outputs |

The replicate stages are the slow ones. With 16 cores the jackknife takes a few minutes and the
bootstrap about half an hour to an hour; the whole pipeline takes several hours. Set `NCORES` to
control the number of worker processes; the pipeline keeps each process to one BLAS thread.

## Notes

- `code/reconstructed/` holds scripts rebuilt from the interactive session in which those analyses
  were first run, each checked against the saved outputs; `compare_outputs.R` and
  `harmonised/compare_harmonised.R` perform that comparison.
- `code/three_step_functions.R` holds the maximum-likelihood three-step estimator: the
  classification-error matrix, the EM algorithm with a Breslow baseline hazard, the bivariate
  residuals and Rubin's rules.
- Years since diagnosis: the interactive build first read the age at diagnosis for thyroid, uterine
  and other cancers from columns that do not exist in the 2005–2016 files. `build_cohort.py` reads
  the correct columns, and every analysis was rerun on the corrected data (2026-09-25).
- `build_docx.py`, `build_journal_variant.py`, `build_references.py`, `fill_strobe_pages.py` and
  `sync_bilingual.py` assemble the submission files from the manuscript source, which is not part
  of this repository.
