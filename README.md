# PLATO-Inspired Synthetic Data for Multi-Regional Clinical Trial Research

## Overview

This repository provides a fully synthetic individual-level multi-regional clinical trial (MRCT) dataset inspired by the PLATO trial (NCT00391872) and calibrated to selected published aggregate summaries. The dataset is intended for methodological research, software development, and demonstration of statistical methods for MRCTs, particularly methods for investigating regional treatment-effect heterogeneity.

The synthetic dataset represents a trial comparing ticagrelor with clopidogrel, with the United States (US) versus the rest of the world (RoW) as the primary regional comparison. It is entirely synthetic and contains no individual-level data from the original PLATO trial.

The repository also accompanies the methodological work presented in *A Workflow for Evaluating Regional Treatment Effect Heterogeneity in Multi-Regional Clinical Trials* and provides a reproducible implementation of the proposed workflow using the PLATO-inspired synthetic dataset:

> Zhang C, et al. *A Workflow for Evaluating Regional Treatment Effect Heterogeneity in Multi-Regional Clinical Trials*. arXiv:2605.16885.  
> https://doi.org/10.48550/arXiv.2605.16885

## Repository structure

```text
.
├── README.md
├── LICENSE
│
├── synthetic_data/
│   ├── DATA_GENERATION.md
│   ├── generate_plato_synthetic.R
│   ├── compare_summaries.R
│   ├── used_prompt.txt
│   └── data/
│       ├── plato_synthetic.csv
│       └── comparison_summary.csv
│
├── analysis/
│   ├── plato_mrct_workflow.R
│   ├── analysis_util.R
│   └── watch_functions.R
│
├── results/
│   └── ...
│
└── figures/
    ├── plato_fidelity.pdf
    └── figure_plato_workflow.pdf
```

`synthetic_data/ai_prompt.md` contains the prompt used with Claude Opus 4.8 (Anthropic) to assist development and validation of the synthetic-data generation code. It is included for transparency and reproducibility.

## Synthetic dataset

The synthetic dataset contains 18,624 observations and was calibrated to selected published aggregate summaries from the PLATO trial, including treatment allocation, regional sample sizes, primary-endpoint event counts, maintenance-aspirin-dose summaries, selected baseline covariate distributions, and published time-to-event summaries.

The main regional comparison is the US versus RoW. The dataset includes treatment assignment, time-to-event outcome information, maintenance aspirin dose, and demographic and clinical covariates used in the PLATO-inspired illustration.

The data can be loaded directly in R:

```r
dat <- read.csv("synthetic_data/data/plato_synthetic.csv")
```

Detailed information on the calibration targets, generation procedure, fidelity assessment, assumptions, and limitations is provided in [`synthetic_data/DATA_GENERATION.md`](synthetic_data/DATA_GENERATION.md).

## Reproducing the synthetic data

The supplied dataset can be used directly. To regenerate it from the prespecified data-generation procedure, run from the repository root:

```bash
Rscript synthetic_data/generate_plato_synthetic.R
```

A fixed random-number seed is used to make the generation procedure reproducible.

To reproduce the fidelity assessment comparing the synthetic dataset with the published aggregate targets, run:

```bash
Rscript synthetic_data/compare_summaries.R
```

This produces the numerical comparison summary and the corresponding fidelity figure.

## AI-assisted code development

Claude Opus 4.8 (Anthropic) was used to assist with development and validation of the code used to generate the PLATO-inspired synthetic dataset. The prompt used for this assistance is provided in [`synthetic_data/ai_prompt.md`](synthetic_data/ai_prompt.md).

The final generation procedure, calibration targets, code, and outputs were reviewed by the authors. No individual-level PLATO data were used in developing or generating the synthetic dataset.

## Example MRCT regional heterogeneity analysis

The repository contains a worked example applying the Q1--Q4 regional heterogeneity workflow to the synthetic dataset:

- **Q1:** Is there evidence of regional treatment-effect heterogeneity?
- **Q2:** Which measured covariates differ across regions?
- **Q3:** Which measured covariates are associated with treatment-effect heterogeneity?
- **Q4:** How do treatment effects vary along leading region-associated candidate covariates?

From the repository root, run:

```bash
Rscript analysis/plato_mrct_workflow.R
```

The analysis produces numerical results in `results/` and the main workflow figure in `figures/`.

The full conditional-random-forest analysis may take some time to run. A reduced quick-run option is available in `analysis/plato_mrct_workflow.R` for code checking and exploratory use; the full setting should be used to reproduce the reported analysis.

## Intended use and limitations

This dataset is intended as a methodological example and benchmark dataset for MRCT research. It may be useful for developing or comparing methods for regional heterogeneity assessment, covariate prioritization, subgroup analysis, visualization, and related methodological problems.

The data preserve selected published aggregate features of PLATO but do not reconstruct the unobserved patient-level data or joint distribution from the original trial. Individual event histories and the joint covariate dependence structure are synthetic. Results obtained from this dataset should therefore not be interpreted as a reanalysis of individual-level PLATO data or as evidence for clinical or causal conclusions regarding the original trial.

In particular, maintenance aspirin dose in PLATO was selected after randomization and should not be interpreted as a randomized baseline treatment factor.

## References

The synthetic dataset was calibrated using published aggregate information from the PLATO trial, including the following reports:

1. Wallentin L, Becker RC, Budaj A, et al.; PLATO Investigators. **Ticagrelor versus clopidogrel in patients with acute coronary syndromes.** *N Engl J Med.* 2009;361(11):1045-1057. doi: [10.1056/NEJMoa0904327](https://doi.org/10.1056/NEJMoa0904327).

2. Mahaffey KW, Wojdyla DM, Carroll K, et al. **Ticagrelor compared with clopidogrel by geographic region in the Platelet Inhibition and Patient Outcomes (PLATO) trial.** *Circulation.* 2011;124(5):544-554.

Additional details on the specific summaries used for calibration are provided in [`synthetic_data/DATA_GENERATION.md`](synthetic_data/DATA_GENERATION.md).

## Citation

If you use the synthetic dataset or the accompanying MRCT workflow, please cite:

> Zhang C, et al. *A Workflow for Evaluating Regional Treatment Effect Heterogeneity in Multi-Regional Clinical Trials*. arXiv:2605.16885.  
> https://doi.org/10.48550/arXiv.2605.16885

## License

See the repository license files for reuse conditions.
