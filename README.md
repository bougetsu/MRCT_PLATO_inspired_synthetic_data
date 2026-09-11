# PLATO-like Synthetic Data for Multi-Regional Clinical Trial Research

## Overview

This repository provides a fully synthetic individual-level dataset based on selected published aggregate summaries from the PLATO trial. The dataset is intended for methodological research, software development, and demonstration of statistical methods for multi-regional clinical trials (MRCTs), particularly methods for investigating regional treatment-effect heterogeneity.

The repository also includes a worked example implementing the Q1--Q4 workflow described in:

> Zhang C, et al. *A Workflow for Evaluating Regional Treatment Effect Heterogeneity in Multi-Regional Clinical Trials*. arXiv:2605.16885.  
> https://doi.org/10.48550/arXiv.2605.16885

The dataset is entirely synthetic and contains no individual-level data from the original PLATO trial.

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

## Synthetic dataset

The synthetic dataset contains 18,624 observations and was calibrated to selected published aggregate results from the PLATO trial, including treatment allocation, regional sample sizes, primary-endpoint event counts, maintenance-aspirin-dose summaries, selected baseline covariate distributions, and published time-to-event summaries.

The main regional comparison is the United States (US) versus the rest of the world (RoW).

The dataset includes treatment assignment, time-to-event outcome information, maintenance aspirin dose, and demographic and clinical covariates used in the PLATO-based illustration.

The data can be loaded directly in R:

```r
dat <- read.csv("synthetic_data/data/plato_synthetic.csv")
```

Detailed information on the construction of the dataset, calibration targets, generation procedure, fidelity assessment, and limitations is provided in [`synthetic_data/DATA_GENERATION.md`](synthetic_data/DATA_GENERATION.md).

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

## Example MRCT regional heterogeneity analysis

The repository contains a worked example applying the Q1--Q4 regional heterogeneity workflow to the synthetic dataset:

- **Q1:** Is there evidence of regional treatment-effect heterogeneity?
- **Q2:** Which measured covariates differ across regions?
- **Q3:** Which measured covariates modify the treatment effect?
- **Q4:** How do treatment effects vary along leading region-associated effect-modifier candidates?

From the repository root, run:

```bash
Rscript analysis/plato_mrct_workflow.R
```

The analysis produces numerical results in `results/` and the main workflow figure in `figures/`.

The full conditional-random-forest analysis may take some time to run. A reduced quick-run option is available in `analysis/plato_mrct_workflow.R` for code checking and exploratory use; the full setting should be used to reproduce the reported analysis.

## Intended use and limitations

This dataset is intended as a methodological example and benchmark dataset for MRCT research. It may be useful for developing or comparing methods for regional heterogeneity assessment, effect-modifier identification, subgroup analysis, visualization, and related methodological problems.

Because the data were generated from aggregate published information, the individual event histories and joint covariate dependence structure are synthetic rather than observed. Results obtained from this dataset should therefore not be interpreted as a reanalysis of individual-level PLATO data or as evidence for clinical or causal conclusions regarding the original trial.

In particular, maintenance aspirin dose in PLATO was selected after randomization and should not be interpreted as a randomized baseline treatment factor.

## References

The synthetic dataset was calibrated using published aggregate information from the PLATO trial, including the following reports:

1. Wallentin L, Becker RC, Budaj A, et al.; PLATO Investigators. **Ticagrelor versus clopidogrel in patients with acute coronary syndromes.** *N Engl J Med.* 2009;361(11):1045-1057. doi: [10.1056/NEJMoa0904327](https://doi.org/10.1056/NEJMoa0904327).

2. Mahaffey KW, Held C, Wojdyla DM, et al.; PLATO Investigators. **Ticagrelor effects on myocardial infarction and the impact of event adjudication in the PLATO (Platelet Inhibition and Patient Outcomes) trial.** *J Am Coll Cardiol.* 2014;63(15):1493-1499. doi: [10.1016/j.jacc.2014.01.038](https://doi.org/10.1016/j.jacc.2014.01.038).

Additional details on the specific summaries used for calibration are provided in [`synthetic_data/DATA_GENERATION.md`](synthetic_data/DATA_GENERATION.md).

## Citation

If you use the synthetic dataset or the accompanying MRCT workflow, please cite:

> Zhang C, et al. *A Workflow for Evaluating Regional Treatment Effect Heterogeneity in Multi-Regional Clinical Trials*. arXiv:2605.16885.  
> https://doi.org/10.48550/arXiv.2605.16885

## License

See the repository license files for reuse conditions.
