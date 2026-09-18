# PLATO-Inspired Synthetic MRCT Data — Generation and Provenance

## 1. Purpose

This repository provides a fully synthetic individual-level MRCT dataset inspired by the PLATO trial (NCT00391872) and calibrated to selected published aggregate summaries. It was created for methodological research and software development in MRCTs, including research on regional treatment-effect heterogeneity.

The synthetic dataset represents a trial comparing ticagrelor with clopidogrel, with the United States (US) versus the rest of the world (RoW) as the primary regional comparison. It should not be interpreted as a reconstruction of the original PLATO individual-level dataset.

No PLATO patient-level records were used. The dataset was generated from published aggregate summaries together with the modelling assumptions described below.

## 2. Published information used

The numerical targets were taken from two published PLATO reports:

| Ref | Source | Information used |
|---|---|---|
| **W** | Wallentin L, et al. *Ticagrelor versus Clopidogrel in Patients with Acute Coronary Syndromes.* **N Engl J Med** 2009;361:1045–1057. | Treatment-arm sizes; overall primary-endpoint event counts; 12-month cumulative incidence and overall hazard ratio; Figure 1 number-at-risk table and cumulative-incidence shape. |
| **M** | Mahaffey KW, et al. *Ticagrelor Compared With Clopidogrel by Geographic Region in the PLATO Trial.* **Circulation** 2011;124:544–554. | US-versus-RoW baseline distributions (Table 1); regional primary-endpoint results (Table 2); region × maintenance-aspirin-dose sample sizes, events, and hazard ratios (Figure 3); candidate effect modifiers considered in Figure 2A. |

The fixed reference values used by the generator are documented directly in the R scripts with the corresponding table or figure source.

## 3. Repository files

```text
synthetic_data/
    DATA_GENERATION.md
    generate_plato_synthetic.R
    compare_summaries.R
    used_prompt.txt
    data/
        plato_synthetic.csv
        comparison_summary.csv

figures/
    plato_fidelity.pdf
    plato_fidelity.png
```

`ai_prompt.md` contains the prompt used with Claude Opus 4.8 (Anthropic) to assist development and validation of the synthetic-data generation code. It is included as part of the provenance record for transparency and reproducibility.

From the repository root:

```bash
Rscript synthetic_data/generate_plato_synthetic.R
Rscript synthetic_data/compare_summaries.R
```

The first command generates `synthetic_data/data/plato_synthetic.csv`. The second recomputes the source-versus-synthetic audit table and produces the two-panel fidelity figure used in the Supplement.

The generator requires `mvtnorm`. The comparison script requires `survival` and, for figures, `ggplot2`, `ggsci`, and `patchwork`. Generation uses the fixed seed `set.seed(20090910)`.

### AI-assisted code development

Claude Opus 4.8 (Anthropic) was used to assist with development and validation of the code used to generate the PLATO-inspired synthetic dataset. The prompt used for this assistance is provided in `synthetic_data/ai_prompt.md`.

The final generation procedure, calibration targets, code, and outputs were reviewed by the authors. No individual-level PLATO data were provided to or used by the AI system.

## 4. Variables

| Column | Type | Coding / target source |
|---|---|---|
| `patient_id` | ID | `P#####` |
| `treatment` | factor | `ticagrelor` / `clopidogrel`; W |
| `region` | factor | `US` / `ROW`; M Table 1/2 |
| `asa_dose` | factor | `low` (≤100 mg), `mid` (>100–<300 mg), `high` (≥300 mg), `undetermined`; M Figure 3 |
| `age` | numeric | Region-specific median/IQR; M Table 1 |
| `sex` | factor | Region-specific proportion; M Table 1 |
| `race` | factor | `White` / `Black` / `Asian` / `Other`; M Table 1 |
| `weight_kg` | numeric | Region-specific median/IQR; M Table 1 |
| `bmi` | numeric | Region-specific median/IQR; M Table 1 |
| `smoking` | factor | `Never` / `Ex` / `Current`; M Table 1 |
| `diabetes` | 0/1 | M Table 1 |
| `prior_mi` | 0/1 | M Table 1 |
| `prior_pci` | 0/1 | M Table 1 |
| `prior_cabg` | 0/1 | M Table 1 |
| `beta_blocker` | 0/1 | M Table 1 |
| `stemi` | 0/1 | M Table 1 |
| `troponin_pos` | 0/1 | M Table 1 |
| `time_days` | numeric | Event/censoring time calibrated to W Figure 1 |
| `event` | 0/1 | Primary composite endpoint; W and M |

The baseline covariates are the intersection of variables reported in Mahaffey Table 1 and those included in the Figure 2A effect-modifier assessment: age, sex, race, weight, BMI, prior MI, prior PCI, prior CABG, smoking, diabetes, troponin, beta-blocker use, and STEMI index event.

## 5. Data generation

### Fixed region × aspirin-dose × treatment strata

Patients are first assigned to a fixed cross-tabulation of region, maintenance aspirin dose, and randomized treatment. For patients with a defined median maintenance aspirin dose, cell sample sizes and primary-endpoint event counts are taken from Mahaffey Figure 3. An `undetermined` aspirin-dose category reconciles these cells with the complete regional and treatment-arm totals reported in the published analyses.

Because sample sizes and event counts are fixed within these cells, the overall, treatment-arm, regional, and region × aspirin-dose sample sizes and event counts are reproduced exactly. Cox hazard ratios also depend on event timing and censoring and are therefore checked for close, rather than exact, agreement.

### Baseline covariates

The 13 baseline covariates are generated separately within US and RoW using a Gaussian copula. Continuous margins are calibrated to the published regional medians and interquartile ranges. Binary and categorical margins are assigned by rank thresholding so that the realized regional category counts match the published targets to rounding or exactly, as appropriate.

The latent correlation matrix specifies plausible dependence among the baseline variables, including stronger weight–BMI association and moderate dependence within prior cardiovascular-history variables. These correlations are modelling assumptions because individual-level joint covariate distributions are not reported in the PLATO publications. The matrix is defined explicitly in `generate_plato_synthetic.R`.

Covariates are generated independently of the aspirin-dose and treatment-cell assignment, conditional on region.

### Event allocation and time-to-event data

Within each region × aspirin-dose × treatment cell, the fixed number of endpoint events is allocated using a common prognostic score based on age, diabetes, prior MI, prior CABG, troponin status, STEMI, and weight. The same prognostic model is used for both randomized treatment groups, giving these variables prognostic signal without deliberately introducing additional treatment interactions.

Event times are generated to follow the cumulative-incidence shape in Wallentin Figure 1. Administrative censoring is calibrated to the corresponding published numbers at risk through 12 months. The maximum analysis follow-up is 360 days.

## 6. Fidelity assessment

`synthetic_data/compare_summaries.R` recomputes the published targets and writes `synthetic_data/data/comparison_summary.csv`. The main checks are:

- all prespecified sample sizes and cell event counts match exactly;
- the overall Cox hazard ratio is 0.84 (0.77–0.92) in both the published and synthetic data;
- regional hazard ratios are approximately 1.27 in the US and 0.81 in RoW, closely matching the published estimates;
- the region × aspirin-dose hazard-ratio pattern is closely reproduced;
- 12-month cumulative incidence is approximately 9.8% for ticagrelor and 11.5% for clopidogrel, compared with published values of 9.8% and 11.7%;
- numbers at risk at 2, 4, 6, 8, 10, and 12 months are close to the published Figure 1 values;
- regional baseline distributions reproduce the published Table 1 margins;
- a diagnostic treatment-by-covariate interaction screen includes all candidate variables used in the synthetic workflow analysis.

The script also produces `figures/plato_fidelity.pdf` (and a PNG copy). Panel A compares published and synthetic hazard ratios with 95% confidence intervals. Panel B shows synthetic cumulative-incidence curves, with the published 12-month values shown as horizontal reference lines. Both panels use the NPG colour palette used elsewhere in the manuscript figures.

## 7. Limitations

The synthetic data preserve selected published aggregate features rather than the unobserved patient-level joint distribution from the original PLATO trial. In particular, dependence among baseline covariates is specified rather than estimated from PLATO individual-level data, and individual event histories are generated from aggregate event-timing and number-at-risk information. Only the primary composite endpoint is modelled.

The `undetermined` aspirin-dose category is used to reconcile Mahaffey Figure 3 with the complete trial totals. The resulting dataset is intended for methodological research and software evaluation, not for clinical inference about the original PLATO trial.
