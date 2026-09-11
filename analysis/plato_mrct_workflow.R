###############################################################################
## analysis/plato_mrct_workflow.R
##
## PLATO synthetic data:
## MRCT analysis using the model-based WATCH TEH workflow.
##
## Workflow
## --------
## Step 0:
##   Reproduce the overall and regional PLATO treatment effects using Cox
##   regression. This is a fidelity and consistency check and is not included
##   in the combined WATCH global-test table.
##
## Step 1:
##   Construct the model-based treatment score residual phi using a prognostic
##   ridge Cox risk score and a homogeneous treatment-effect Cox model.
##
## Step 2:
##   Test phi versus Region using a coin maximum-type independence test.
##
## Step 3:
##   Assess Region versus X using:
##     - A coin maximum-type global independence test.
##     - Conditional random forest permutation variable importance.
##
## Step 4:
##   Assess treatment-effect modification using:
##     - A coin maximum-type global independence test of phi versus X.
##     - Conditional random forest permutation variable importance for
##       phi versus X plus Region.
##
## Q4 displays:
##   - Mean phi and 95 percent CI by aspirin-dose category and region.
##   - Aspirin-dose distribution within each region.
##   - Cox hazard ratios by aspirin-dose category and region.
##
## Reproducibility
## ---------------
## Every computationally intensive step is checkpointed with saveRDS() under
## .cache/. These checkpoints are local computational artifacts and need not be
## committed to the public repository. Existing checkpoints are loaded unless:
##   - The corresponding FORCE_RERUN flag is TRUE.
##   - The input data, covariate set, or analysis configuration has changed.
###############################################################################


## ============================================================================
## 0. SETUP
## ============================================================================

rm(list = ls())

SEED <- 42L
ANALYSIS_VERSION <- "plato_workflow_v3"
set.seed(SEED)

## The script is expected to be run from the project root. here::here() also
## permits execution from other working directories when the project has a
## .git directory or .here file.
if (!requireNamespace("here", quietly = TRUE)) {
  stop(
    "Package 'here' is required. Install it before running this script.",
    call. = FALSE
  )
}

required_packages <- c(
  "dplyr",
  "ggplot2",
  "survival",
  "coin",
  "party",
  "permimp",
  "glmnet",
  "sandwich",
  "MASS",
  "patchwork",
  "scales"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "The following packages are required but not installed: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

rm(required_packages, missing_packages)


## ============================================================================
## 1. PATHS AND ANALYSIS CONFIGURATION
## ============================================================================

data_path <- here::here(
  "synthetic_data",
  "data",
  "plato_synthetic.csv"
)

watch_functions_path <- here::here(
  "analysis",
  "watch_functions.R"
)

utils_path <- here::here(
  "analysis",
  "analysis_util.R"
)

if (!file.exists(data_path)) {
  stop(
    "Synthetic data file not found: ",
    data_path,
    call. = FALSE
  )
}

if (!file.exists(watch_functions_path)) {
  stop(
    "WATCH function file not found: ",
    watch_functions_path,
    call. = FALSE
  )
}

if (!file.exists(utils_path)) {
  stop(
    "Utility function file not found: ",
    utils_path,
    call. = FALSE
  )
}

source(watch_functions_path)
source(utils_path)


## Quick versus full CRF analysis.
##
## QUICK:
##   ntree = 100
##   nperm = 5
##
## FULL:
##   ntree = 500
##   nperm = 10
FULL_RUN <- TRUE

RUN_LABEL <- if (FULL_RUN) {
  "full"
} else {
  "quick"
}

cfg <- if (FULL_RUN) {
  list(
    CRF_MTRY_TEH = 5L,
    CRF_MTRY_REGION = 5L,
    CRF_NTREE = 500L,
    CRF_NPERM = 10L
  )
} else {
  list(
    CRF_MTRY_TEH = 5L,
    CRF_MTRY_REGION = 5L,
    CRF_NTREE = 100L,
    CRF_NPERM = 5L
  )
}

## Global tests use the asymptotic reference distribution in the main PLATO
## analysis. An optional Monte Carlo permutation sensitivity check is provided
## below and is disabled by default.
GLOBAL_TEST_DISTRIBUTION <- "asymptotic"
RUN_PERMUTATION_SENSITIVITY <- FALSE
PERM_NRESAMPLE <- 9999L

## NPG palette colours used in the manuscript figures.
REGION_COLOURS <- c(
  "RoW" = "#E64B35FF",
  "US" = "#4DBBD5FF"
)
HIGHLIGHT_COLOUR <- "#00A087FF"


## Public outputs are kept separate from computational checkpoints.
results_dir <- here::here(
  "results",
  RUN_LABEL
)

fig_dir <- here::here(
  "figures",
  RUN_LABEL
)

cache_dir <- here::here(
  ".cache",
  "plato_workflow",
  RUN_LABEL
)

for (path in c(results_dir, fig_dir, cache_dir)) {
  dir.create(
    path,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


## Manual rerun controls.
##
## The script will also rerun all steps automatically when the input data,
## covariate set, or CRF configuration differs from the saved manifest.
FORCE_RERUN <- c(
  step0 = FALSE,
  step1 = FALSE,
  step2 = FALSE,
  step3 = FALSE,
  step4 = FALSE
)


## ============================================================================
## 2. READ AND PREPARE THE DATA
## ============================================================================

d <- utils::read.csv(
  data_path,
  stringsAsFactors = TRUE
)

required_data_columns <- c(
  "treatment",
  "time_days",
  "event",
  "region",
  "asa_dose"
)

missing_data_columns <- setdiff(
  required_data_columns,
  names(d)
)

if (length(missing_data_columns) > 0L) {
  stop(
    "The following required columns are absent from the synthetic data: ",
    paste(missing_data_columns, collapse = ", "),
    call. = FALSE
  )
}


d$event <- as.integer(
  d$event
)

if (!all(d$event %in% c(0L, 1L))) {
  stop(
    "'event' must be coded as 0 or 1.",
    call. = FALSE
  )
}


## Set reference categories explicitly.
d$treatment <- factor(
  d$treatment,
  levels = c(
    "clopidogrel",
    "ticagrelor"
  )
)

## The synthetic-data file uses "ROW"; use "RoW" only as the display label
## in the analysis outputs and figures.
d$region <- as.character(d$region)
d$region[d$region == "ROW"] <- "RoW"
d$region <- factor(
  d$region,
  levels = c(
    "RoW",
    "US"
  )
)


## Use a clinically meaningful order for aspirin-dose categories.
d$asa_dose <- factor(
  d$asa_dose,
  levels = c(
    "undetermined",
    "low",
    "mid",
    "high"
  ),
  ordered = TRUE
)

if (anyNA(d$treatment)) {
  stop(
    "Unexpected treatment labels were found. Expected: clopidogrel and ",
    "ticagrelor.",
    call. = FALSE
  )
}

if (anyNA(d$region)) {
  stop(
    "Unexpected region labels were found. Expected: RoW and US.",
    call. = FALSE
  )
}

if (anyNA(d$asa_dose)) {
  stop(
    "Unexpected aspirin-dose labels were found. Expected: undetermined, low, ",
    "mid, and high.",
    call. = FALSE
  )
}


## Candidate variables used in the PLATO workflow.
##
## X contains the observed candidate variables used for Q2 and the Q3 global
## test. Region is added separately to the Q3 CRF ranking. Maintenance aspirin
## dose is retained as a candidate variable although it was determined after
## randomization in PLATO.
candidate_covariates <- c(
  "age",
  "sex",
  "race",
  "weight_kg",
  "bmi",
  "smoking",
  "diabetes",
  "prior_mi",
  "prior_pci",
  "prior_cabg",
  "beta_blocker",
  "stemi",
  "troponin_pos",
  "asa_dose"
)

required_analysis_covariates <- c(
  candidate_covariates,
  "region"
)

missing_analysis_covariates <- setdiff(
  required_analysis_covariates,
  names(d)
)

if (length(missing_analysis_covariates) > 0L) {
  stop(
    "The following PLATO analysis covariates are missing from the synthetic data: ",
    paste(missing_analysis_covariates, collapse = ", "),
    call. = FALSE
  )
}

covariates <- required_analysis_covariates


## Construct the exact data frame passed to the model-based TEH analysis.
dat_teh <- d[
  ,
  c(
    "treatment",
    "time_days",
    "event",
    covariates
  ),
  drop = FALSE
]


## The extracted WATCH input validator accepts numeric or factor covariates.
## Convert any residual character covariates to factors.
character_covariates <- vapply(
  dat_teh,
  is.character,
  FUN.VALUE = logical(1)
)

dat_teh[
  ,
  character_covariates
] <- lapply(
  dat_teh[
    ,
    character_covariates,
    drop = FALSE
  ],
  factor
)


## Complete-case restriction for the selected analysis variables.
##
## IDA and preprocessing should normally have handled missing covariates before
## this script. This restriction protects the downstream WATCH functions.
n_before <- nrow(dat_teh)

complete_rows <- stats::complete.cases(
  dat_teh
)

dat_teh <- dat_teh[
  complete_rows,
  ,
  drop = FALSE
]

n_after <- nrow(dat_teh)

if (n_after < n_before) {
  warning(
    "Dropped ",
    n_before - n_after,
    " rows with missing analysis variables.",
    call. = FALSE
  )
}


## X excludes Region because Region is the response in Step 3.
X <- dat_teh[
  ,
  candidate_covariates,
  drop = FALSE
]

Region <- factor(
  dat_teh$region,
  levels = c(
    "RoW",
    "US"
  )
)


## Verify that all analysis covariates are numeric or factor.
valid_covariates <- vapply(
  dat_teh[
    ,
    covariates,
    drop = FALSE
  ],
  function(column) {
    is.numeric(column) || is.factor(column)
  },
  FUN.VALUE = logical(1)
)

if (!all(valid_covariates)) {
  stop(
    "The following analysis covariates are neither numeric nor factor: ",
    paste(
      names(valid_covariates)[!valid_covariates],
      collapse = ", "
    ),
    call. = FALSE
  )
}


## Save prepared analysis data.
prepared_path <- file.path(
  cache_dir,
  "00_prepared_data.rds"
)

saveRDS(
  list(
    dat_teh = dat_teh,
    X = X,
    Region = Region,
    covariates = covariates,
    cfg = cfg,
    full_run = FULL_RUN,
    run_label = RUN_LABEL
  ),
  file = prepared_path
)


cat("\n==== Prepared analysis data ====\n")
cat("Rows:", nrow(dat_teh), "\n")
cat("Columns:", ncol(dat_teh), "\n")
cat("Minimum follow-up time:", min(dat_teh$time_days), "\n")

cat("\nTreatment counts:\n")
print(
  table(dat_teh$treatment)
)

cat("\nRegion counts:\n")
print(
  table(dat_teh$region)
)

cat("\nAspirin-dose counts:\n")
print(
  table(dat_teh$asa_dose)
)

cat("\nAspirin dose by region:\n")
print(
  table(
    dat_teh$region,
    dat_teh$asa_dose
  )
)


## ============================================================================
## 3. CHECK WHETHER DATA OR CONFIGURATION HAS CHANGED
## ============================================================================

manifest_path <- file.path(
  cache_dir,
  "analysis_manifest.rds"
)

current_manifest <- list(
  analysis_version = ANALYSIS_VERSION,
  data_path = normalizePath(
    data_path,
    winslash = "/",
    mustWork = TRUE
  ),
  data_md5 = unname(
    tools::md5sum(data_path)
  ),
  covariates = covariates,
  candidate_covariates = candidate_covariates,
  cfg = cfg,
  global_test_distribution = GLOBAL_TEST_DISTRIBUTION,
  full_run = FULL_RUN,
  seed = SEED
)

previous_manifest <- if (file.exists(manifest_path)) {
  readRDS(manifest_path)
} else {
  NULL
}

analysis_changed <- (
  is.null(previous_manifest) ||
    !identical(
      previous_manifest$analysis_version,
      current_manifest$analysis_version
    ) ||
    !identical(
      previous_manifest$data_md5,
      current_manifest$data_md5
    ) ||
    !identical(
      previous_manifest$covariates,
      current_manifest$covariates
    ) ||
    !identical(
      previous_manifest$candidate_covariates,
      current_manifest$candidate_covariates
    ) ||
    !identical(
      previous_manifest$cfg,
      current_manifest$cfg
    ) ||
    !identical(
      previous_manifest$global_test_distribution,
      current_manifest$global_test_distribution
    )
)

if (analysis_changed) {
  message(
    "Input data, covariates, or analysis configuration changed. ",
    "Existing step checkpoints will be recomputed."
  )
}

## Save the current manifest immediately. If the session stops later, completed
## step checkpoints can be loaded on the next run.
saveRDS(
  current_manifest,
  file = manifest_path
)

force_step <- function(step_name) {
  isTRUE(FORCE_RERUN[[step_name]]) ||
    isTRUE(analysis_changed)
}


## ============================================================================
## STEP 0. COX REGIONAL CONSISTENCY
## ============================================================================

step0_path <- file.path(
  cache_dir,
  "step0_cox_regional_consistency.rds"
)

step0 <- run_or_load(
  path = step0_path,
  force = force_step("step0"),
  fun = function() {
    analysis_data <- dat_teh
    
    analysis_data$trt_num <- as.integer(
      analysis_data$treatment == "ticagrelor"
    )
    
    fit_overall <- survival::coxph(
      survival::Surv(time_days, event) ~ trt_num,
      data = analysis_data
    )
    
    fit_region_interaction <- survival::coxph(
      survival::Surv(time_days, event) ~ trt_num * region,
      data = analysis_data
    )
    
    fit_us <- survival::coxph(
      survival::Surv(time_days, event) ~ trt_num,
      data = analysis_data[
        analysis_data$region == "US",
        ,
        drop = FALSE
      ]
    )
    
    fit_RoW <- survival::coxph(
      survival::Surv(time_days, event) ~ trt_num,
      data = analysis_data[
        analysis_data$region == "RoW",
        ,
        drop = FALSE
      ]
    )
    
    extract_hr <- function(
    fit,
    analysis_label,
    term = "trt_num"
    ) {
      fit_summary <- summary(fit)
      
      data.frame(
        analysis = analysis_label,
        HR = unname(
          fit_summary$conf.int[
            term,
            "exp(coef)"
          ]
        ),
        lower_95_CI = unname(
          fit_summary$conf.int[
            term,
            "lower .95"
          ]
        ),
        upper_95_CI = unname(
          fit_summary$conf.int[
            term,
            "upper .95"
          ]
        ),
        p_value = unname(
          fit_summary$coefficients[
            term,
            "Pr(>|z|)"
          ]
        ),
        stringsAsFactors = FALSE
      )
    }
    
    interaction_term <- "trt_num:regionUS"
    
    interaction_summary <- summary(
      fit_region_interaction
    )
    
    interaction_p_value <- unname(
      interaction_summary$coefficients[
        interaction_term,
        "Pr(>|z|)"
      ]
    )
    
    interaction_statistic <- unname(
      interaction_summary$coefficients[
        interaction_term,
        "z"
      ]
    )
    
    hr_table <- dplyr::bind_rows(
      extract_hr(
        fit_overall,
        analysis_label = "Overall"
      ),
      extract_hr(
        fit_us,
        analysis_label = "US"
      ),
      extract_hr(
        fit_RoW,
        analysis_label = "RoW"
      )
    )
    
    hr_table$interaction_p_value <- interaction_p_value
    
    list(
      fit_overall = fit_overall,
      fit_region_interaction = fit_region_interaction,
      fit_us = fit_us,
      fit_RoW = fit_RoW,
      hr_table = hr_table,
      interaction_statistic = interaction_statistic,
      interaction_p_value = interaction_p_value
    )
  }
)

cat("\n==== Step 0: Cox regional consistency ====\n")
print(
  step0$hr_table
)

utils::write.csv(
  step0$hr_table,
  file = file.path(
    results_dir,
    "cox_regional_consistency.csv"
  ),
  row.names = FALSE
)


## ============================================================================
## STEP 1. MODEL-BASED SCORE RESIDUAL PHI AND Q3 CRF RANKING
## ============================================================================

step1_path <- file.path(
  cache_dir,
  "step1_watch_modelbased_result.rds"
)

teh_result <- run_or_load(
  path = step1_path,
  force = force_step("step1"),
  fun = function() {
    set.seed(
      SEED + 101L
    )
    
    exploreTEH_mb(
      data = dat_teh,
      trt_name = "treatment",
      y_name = "time_days",
      event_name = "event",
      drug_name = "ticagrelor",
      control_name = "clopidogrel",
      type = "survival",
      alpha = 0,
      teststat = "maximum",
      mtry = cfg$CRF_MTRY_TEH,
      ntree = cfg$CRF_NTREE,
      nperm = cfg$CRF_NPERM,
      verbose = TRUE,
      include_cforest = FALSE
    )
  }
)

cat("\n==== Step 1: score residual phi and Q3 CRF ranking ====\n")
step1_vi_preview <- get_importance_scores(teh_result)
step1_vi_preview <- step1_vi_preview[
  order(step1_vi_preview$vi, decreasing = TRUE),
  ,
  drop = FALSE
]
print(utils::head(step1_vi_preview, 15), row.names = FALSE)

phi_data <- data.frame(
  phi = teh_result$phi,
  region = Region,
  asa_dose = dat_teh$asa_dose
)

saveRDS(
  phi_data,
  file = file.path(
    cache_dir,
    "step1_phi_data.rds"
  )
)


## ============================================================================
## Q1. PHI VERSUS REGION
## ============================================================================

step2_path <- file.path(
  cache_dir,
  "step2_phi_region_coin.rds"
)

step2 <- run_or_load(
  path = step2_path,
  force = force_step("step2"),
  fun = function() {
    test_data <- data.frame(
      phi = teh_result$phi,
      Region = Region
    )
    
    test_maximum <- coin::independence_test(
      phi ~ Region,
      data = test_data,
      teststat = "maximum",
      distribution = coin::asymptotic()
    )
    
    list(
      test = test_maximum,
      p_value = as.numeric(
        coin::pvalue(test_maximum)
      ),
      statistic = get_coin_scalar_stat(
        test_maximum,
        type = "standardized"
      ),
      statistic_type = "maximum",
      distribution = "asymptotic"
    )
  }
)

cat("\n==== Q1: phi ~ Region ====\n")
cat("Test: coin independence test\n")
cat("Statistic type: maximum\n")
cat("Statistic:", step2$statistic, "\n")
cat("P-value:", step2$p_value, "\n")

step2_summary <- data.frame(
  question = "phi ~ Region",
  test = "coin independence test",
  statistic_type = "maximum",
  statistic = step2$statistic,
  p_value = step2$p_value,
  stringsAsFactors = FALSE
)

utils::write.csv(
  step2_summary,
  file = file.path(
    results_dir,
    "q1_region_global_test.csv"
  ),
  row.names = FALSE
)


## ============================================================================
## Q2. REGION VERSUS X
## ============================================================================

step3_path <- file.path(
  cache_dir,
  "step3_region_x_assess_regional_imbalance.rds"
)

step3 <- run_or_load(
  path = step3_path,
  force = force_step("step3"),
  fun = function() {
    set.seed(
      SEED + 301L
    )
    
    assess_regional_imbalance(
      X = X,
      Region = Region,
      mtry = cfg$CRF_MTRY_REGION,
      ntree = cfg$CRF_NTREE,
      nperm = cfg$CRF_NPERM,
      verbose = TRUE,
      include_cforest = FALSE
    )
  }
)

cat("\n==== Q2: Region ~ X ====\n")
cat("Test: coin independence test\n")
cat("Statistic type:", step3$statistic_type, "\n")
cat("Statistic:", step3$statistic, "\n")
cat("P-value:", step3$region_p_max, "\n")

cat("\nTop 15 region-associated variables:\n")
print(
  utils::head(
    step3$vi_table,
    15
  ),
  row.names = FALSE
)

step3_summary <- data.frame(
  question = "Region ~ X",
  test = "coin independence test",
  statistic_type = "maximum",
  statistic = step3$statistic,
  p_value = step3$region_p_max,
  stringsAsFactors = FALSE
)

utils::write.csv(
  step3_summary,
  file = file.path(
    results_dir,
    "q2_region_global_test.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  step3$vi_table,
  file = file.path(
    results_dir,
    "q2_variable_importance.csv"
  ),
  row.names = FALSE
)


## ============================================================================
## Q3. PHI VERSUS X: GLOBAL TEST; PHI VERSUS X + REGION: CRF RANKING
## ============================================================================

step4_path <- file.path(
  cache_dir,
  "step4_watch_teh_summary.rds"
)

step4 <- run_or_load(
  path = step4_path,
  force = force_step("step4"),
  fun = function() {
    ## Q3 global test follows the manuscript definition: phi versus X.
    ## Region is retained separately in the Q3 CRF ranking.
    q3_test_data <- data.frame(
      phi = teh_result$phi,
      X,
      check.names = FALSE
    )

    q3_global_test <- coin::independence_test(
      phi ~ .,
      data = q3_test_data,
      teststat = "maximum",
      distribution = coin::asymptotic()
    )

    ## Q3 ranking: use the CRF fit from exploreTEH_mb(), which received
    ## baseline covariates plus Region. Region is retained in the ranking as a
    ## diagnostic candidate effect modifier.
    importance_table <- get_importance_scores(
      teh_result
    )

    names(importance_table) <- c(
      "covariate",
      "vi_teh"
    )

    importance_table <- importance_table[
      order(
        importance_table$vi_teh,
        decreasing = TRUE,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    rownames(importance_table) <- NULL

    list(
      test = q3_global_test,
      global_p_value = as.numeric(coin::pvalue(q3_global_test)),
      statistic = get_coin_scalar_stat(
        q3_global_test,
        type = "standardized"
      ),
      statistic_type = "maximum",
      distribution = GLOBAL_TEST_DISTRIBUTION,
      vi_table = importance_table,
      outcome_type = teh_result$outcome_type,
      method = teh_result$method
    )
  }
)

cat("\n==== Q3 global test: phi ~ X ====\n")
cat("Test: coin maximum-type independence test\n")
cat("Reference distribution:", step4$distribution, "\n")
cat("Statistic:", step4$statistic, "\n")
cat("P-value:", step4$global_p_value, "\n")

cat("\nTop 15 Q3 CRF variables (X + Region):\n")
print(
  utils::head(
    step4$vi_table,
    15
  ),
  row.names = FALSE
)

step4_summary <- data.frame(
  question = "phi ~ X",
  test = "coin maximum-type independence test",
  statistic_type = "maximum",
  distribution = step4$distribution,
  statistic = step4$statistic,
  p_value = step4$global_p_value,
  stringsAsFactors = FALSE
)

utils::write.csv(
  step4_summary,
  file = file.path(
    results_dir,
    "q3_teh_global_test.csv"
  ),
  row.names = FALSE
)

utils::write.csv(
  step4$vi_table,
  file = file.path(
    results_dir,
    "q3_variable_importance.csv"
  ),
  row.names = FALSE
)


## ============================================================================
## 4. COMBINED WORKFLOW P-VALUE TABLE
## ============================================================================

workflow_pvalue_summary <- build_workflow_pvalue_summary(
  step2 = step2,
  step3 = step3,
  step4 = step4
)

cat("\n==== Combined WATCH workflow global tests ====\n")
print(
  workflow_pvalue_summary,
  row.names = FALSE
)

utils::write.csv(
  workflow_pvalue_summary,
  file = file.path(
    results_dir,
    "workflow_pvalue_summary.csv"
  ),
  row.names = FALSE
)

saveRDS(
  workflow_pvalue_summary,
  file = file.path(
    cache_dir,
    "workflow_pvalue_summary.rds"
  )
)


## Optional Monte Carlo permutation sensitivity check for the three global tests.
## This is not required for the main PLATO analysis and is disabled by default.
if (isTRUE(RUN_PERMUTATION_SENSITIVITY)) {
  set.seed(SEED + 901L)

  q1_perm_data <- data.frame(
    phi = teh_result$phi,
    Region = Region
  )

  q1_perm <- coin::independence_test(
    phi ~ Region,
    data = q1_perm_data,
    teststat = "maximum",
    distribution = coin::approximate(nresample = PERM_NRESAMPLE)
  )

  q2_perm_data <- data.frame(
    Region = Region,
    X,
    check.names = FALSE
  )

  q2_perm <- coin::independence_test(
    Region ~ .,
    data = q2_perm_data,
    teststat = "maximum",
    distribution = coin::approximate(nresample = PERM_NRESAMPLE)
  )

  q3_perm_data <- data.frame(
    phi = teh_result$phi,
    X,
    check.names = FALSE
  )

  q3_perm <- coin::independence_test(
    phi ~ .,
    data = q3_perm_data,
    teststat = "maximum",
    distribution = coin::approximate(nresample = PERM_NRESAMPLE)
  )

  permutation_sensitivity <- data.frame(
    step = c("Q1", "Q2", "Q3"),
    asymptotic_p = c(
      step2$p_value,
      step3$region_p_max,
      step4$global_p_value
    ),
    permutation_p = c(
      as.numeric(coin::pvalue(q1_perm)),
      as.numeric(coin::pvalue(q2_perm)),
      as.numeric(coin::pvalue(q3_perm))
    ),
    nresample = PERM_NRESAMPLE,
    stringsAsFactors = FALSE
  )

  print(permutation_sensitivity, row.names = FALSE)

  utils::write.csv(
    permutation_sensitivity,
    file = file.path(
      results_dir,
      "global_test_permutation_sensitivity.csv"
    ),
    row.names = FALSE
  )
}


## ============================================================================
## 5. Q2 AND Q3 RANK COMPARISON
## ============================================================================

## Retain the complete Q2 and Q3 rankings. No fixed top-K threshold is used to
## select Q4 candidates; prominence in both rankings is considered together
## with clinical and external evidence.
region_ranking <- step3$vi_table
region_ranking$region_rank <- seq_len(nrow(region_ranking))

teh_ranking <- step4$vi_table
teh_ranking$teh_rank <- seq_len(nrow(teh_ranking))

rank_comparison <- merge(
  region_ranking[
    ,
    c(
      "covariate",
      "vi_region",
      "region_rank"
    ),
    drop = FALSE
  ],
  teh_ranking[
    ,
    c(
      "covariate",
      "vi_teh",
      "teh_rank"
    ),
    drop = FALSE
  ],
  by = "covariate",
  all = TRUE
)

rank_comparison$combined_rank <-
  rank_comparison$region_rank + rank_comparison$teh_rank

rank_comparison <- rank_comparison[
  order(
    rank_comparison$combined_rank,
    rank_comparison$region_rank,
    rank_comparison$teh_rank,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

rownames(rank_comparison) <- NULL

cat("\n==== Q2/Q3 rank comparison ====\n")
print(utils::head(rank_comparison, 15), row.names = FALSE)

utils::write.csv(
  rank_comparison,
  file = file.path(
    results_dir,
    "q2_q3_rank_comparison.csv"
  ),
  row.names = FALSE
)


## ============================================================================
## 6. VARIABLE-IMPORTANCE FIGURES
## ============================================================================

TOP_N_PLOT <- 15L

region_vi_plot <- plot_variable_importance(
  vi_table = step3$vi_table,
  value_col = "vi_region",
  top_n = TOP_N_PLOT,
  highlight = "asa_dose",
  title = "Regional covariate imbalance",
  subtitle = "Conditional random forest variable importance for Region versus X",
  y_label = "Permutation variable importance"
)

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "q2_variable_importance_top15.pdf"
  ),
  plot = region_vi_plot,
  width = 8,
  height = 6,
  units = "in"
)


teh_vi_plot <- plot_variable_importance(
  vi_table = step4$vi_table,
  value_col = "vi_teh",
  top_n = TOP_N_PLOT,
  highlight = c(
    "asa_dose",
    "region"
  ),
  title = "Treatment effect heterogeneity",
  subtitle = paste(
    "Conditional random forest variable importance for",
    "treatment score residual versus Region plus X"
  ),
  y_label = "Permutation variable importance"
)

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "q3_variable_importance_top15.pdf"
  ),
  plot = teh_vi_plot,
  width = 8,
  height = 6,
  units = "in"
)



## ============================================================================
## 7. Q4 ASPIRIN-DOSE DISPLAYS
## ============================================================================

## Q4A: Score residual phi version (workflow-native metric)
q4_phi_plot <- plot_q4_display_categorical(
  data = dat_teh,
  phi = teh_result$phi,
  x_name = "asa_dose",
  region_name = "region",
  x_levels = c(
    "undetermined",
    "low",
    "mid",
    "high"
  ),
  x_labels = c(
    "Low",
    "Intermediate",
    "High"
  ),
  exclude_levels = "undetermined",
  metric = "phi",
  region_colours = REGION_COLOURS,
  title = "Treatment effect heterogeneity and aspirin dose by region"
)

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "q4_asa_phi.pdf"
  ),
  plot = q4_phi_plot,
  width = 8.5,
  height = 8,
  units = "in"
)

## ============================================================================
## 8. TREATMENT HAZARD RATIO BY REGION AND ASPIRIN DOSE
## ============================================================================

hr_by_aspirin <- estimate_hr_by_subgroup(
  data = dat_teh,
  time_name = "time_days",
  event_name = "event",
  treatment_name = "treatment",
  drug_name = "ticagrelor",
  control_name = "clopidogrel",
  region_name = "region",
  subgroup_name = "asa_dose",
  subgroup_levels = c(
    "undetermined",
    "low",
    "mid",
    "high"
  )
)

utils::write.csv(
  hr_by_aspirin,
  file = file.path(
    results_dir,
    "hr_by_region_asa.csv"
  ),
  row.names = FALSE
)

saveRDS(
  hr_by_aspirin,
  file = file.path(
    cache_dir,
    "step4_hr_by_region_and_asa_dose.rds"
  )
)

hr_aspirin_plot <- plot_hr_by_subgroup(
  hr_table = hr_by_aspirin,
  treatment_label = "Hazard ratio for ticagrelor versus clopidogrel",
  title = "Treatment effect by region and aspirin dose",
  subtitle = "Separate Cox models within each region and aspirin-dose category",
  region_colours = REGION_COLOURS,
  log_scale = TRUE
)

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "hr_by_region_asa.pdf"
  ),
  plot = hr_aspirin_plot,
  width = 8.5,
  height = 6,
  units = "in"
)

## Region overall HR from Step 0 Cox models
extract_overall_hr_from_cox <- function(fit, region_label) {
  fit_summary <- summary(fit)
  
  data.frame(
    region = region_label,
    HR = unname(
      fit_summary$conf.int[
        "trt_num",
        "exp(coef)"
      ]
    ),
    lower = unname(
      fit_summary$conf.int[
        "trt_num",
        "lower .95"
      ]
    ),
    upper = unname(
      fit_summary$conf.int[
        "trt_num",
        "upper .95"
      ]
    ),
    stringsAsFactors = FALSE
  )
}

region_overall_hr <- dplyr::bind_rows(
  extract_overall_hr_from_cox(step0$fit_RoW, "RoW"),
  extract_overall_hr_from_cox(step0$fit_us, "US")
)

## Q4B: Hazard ratio version (clinical metric)
q4_hr_plot <- plot_q4_display_categorical(
  data = dat_teh,
  phi = teh_result$phi,
  x_name = "asa_dose",
  region_name = "region",
  x_levels = c(
    "undetermined",
    "low",
    "mid",
    "high"
  ),
  x_labels = c(
    "Low",
    "Intermediate",
    "High"
  ),
  exclude_levels = c("undetermined"),
  metric = "hr",
  hr_table = hr_by_aspirin,
  region_overall_hr = region_overall_hr,
  min_events = 5L,
  y_breaks = c(0.3, 0.5, 0.8, 1, 1.5, 2, 3),
  region_colours = REGION_COLOURS,
  title = "Treatment effect (HR) and aspirin dose by region"
)

ggplot2::ggsave(
  filename = file.path(
    fig_dir,
    "q4_asa_hr.pdf"
  ),
  plot = q4_hr_plot,
  width = 8.5,
  height = 8,
  units = "in"
)


## ============================================================================
## 9. FINAL SUMMARY OBJECT
## ============================================================================

final_summary <- list(
  manifest = current_manifest,
  prepared_data_path = prepared_path,
  
  ## Published-result consistency analysis
  step0_cox = step0$hr_table,
  step0_interaction_statistic = step0$interaction_statistic,
  step0_interaction_p_value = step0$interaction_p_value,
  
  ## WATCH global tests
  workflow_pvalue_summary = workflow_pvalue_summary,
  
  ## Regional imbalance
  step3_top_region_variables = utils::head(
    step3$vi_table,
    20
  ),
  
  ## Treatment effect heterogeneity
  step4_top_teh_variables = utils::head(
    step4$vi_table,
    20
  ),
  
  ## Q2/Q3 complete rank comparison
  q2_q3_rank_comparison = rank_comparison,
  
  ## Aspirin-dose clinical summary
  hr_by_region_and_asa_dose = hr_by_aspirin,
  
  ## Configuration
  cfg = cfg,
  full_run = FULL_RUN,
  run_label = RUN_LABEL,
  seed = SEED
)

saveRDS(
  final_summary,
  file = file.path(
    cache_dir,
    "final_main_workflow_summary.rds"
  )
)



TOP_N_PLOT <- 10L

region_vi_plot_ms <- plot_variable_importance(
  vi_table = step3$vi_table,
  value_col = "vi_region",
  top_n = TOP_N_PLOT,
  highlight = "asa_dose",
  title = "Q2: Regional covariate ranking",
  subtitle = NULL,
  y_label = "Permutation variable importance",
  highlight_colour = HIGHLIGHT_COLOUR,
  other_colour = "grey75"
) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    axis.title = ggplot2::element_text(size = 9.5),
    axis.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(5.5, 8, 5.5, 5.5)
  ) +ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 10.5),
    axis.text.x = ggplot2::element_text(size = 9),
    axis.title.x = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(
      size = 12,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = ggplot2::element_text(size = 9)
  )

teh_vi_plot_ms <- plot_variable_importance(
  vi_table = step4$vi_table,
  value_col = "vi_teh",
  top_n = TOP_N_PLOT,
  highlight = c("asa_dose", "region"),
  title = "Q3: Effect modifier ranking",
  subtitle = NULL,
  y_label = "Permutation variable importance",
  highlight_colour = HIGHLIGHT_COLOUR,
  other_colour = "grey75"
) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      size = 11
    ),
    axis.title = ggplot2::element_text(size = 9.5),
    axis.text = ggplot2::element_text(size = 9),
    plot.margin = ggplot2::margin(5.5, 8, 5.5, 5.5)
  )+ggplot2::theme(
    panel.grid.major.y = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = 10.5),
    axis.text.x = ggplot2::element_text(size = 9),
    axis.title.x = ggplot2::element_text(size = 10),
    plot.title = ggplot2::element_text(
      size = 12,
      face = "bold",
      hjust = 0
    ),
    plot.subtitle = ggplot2::element_text(size = 9)
  )


q4_hr_plot <- plot_q4_display_categorical(
  data = dat_teh,
  phi = teh_result$phi,
  x_name = "asa_dose",
  region_name = "region",
  x_levels = c(
    "undetermined",
    "low",
    "mid",
    "high"
  ),
  x_labels = c(
    "Low",
    "Intermediate",
    "High"
  ),
  exclude_levels = "undetermined",
  metric = "hr",
  hr_table = hr_by_aspirin,
  region_overall_hr = region_overall_hr,
  y_breaks = c(0.3, 0.5, 0.8, 1, 1.5, 2, 3),
  min_events = 5L,
  region_colours = REGION_COLOURS,
  title = "Q4: Treatment effects and aspirin dose by region"
)



title_theme <- ggplot2::theme(
  plot.title = ggplot2::element_text(
    size = 12,
    face = "bold",
    hjust = 0,
    margin = ggplot2::margin(b = 5)
  )
)


region_vi_plot_ms <- region_vi_plot_ms + title_theme
teh_vi_plot_ms <- teh_vi_plot_ms + title_theme


region_vi_plot_ms <- region_vi_plot_ms +
  ggplot2::labs(tag = "A")

teh_vi_plot_ms <- teh_vi_plot_ms +
  ggplot2::labs(tag = "B")

panel_tag_theme <- ggplot2::theme(
  plot.tag = ggplot2::element_text(
    size = 14,
    face = "bold"
  ),
  plot.tag.position = c(0, 1)
)

region_vi_plot_ms <- region_vi_plot_ms + panel_tag_theme
teh_vi_plot_ms <- teh_vi_plot_ms + panel_tag_theme

q4_hr_plot[[1]] <- q4_hr_plot[[1]] +
  ggplot2::labs(tag = "C") +
  panel_tag_theme

q4_hr_plot_ms <- q4_hr_plot & title_theme


library(patchwork)

left_vi <- patchwork::wrap_plots(
  region_vi_plot_ms,
  teh_vi_plot_ms,
  ncol = 1,
  heights = c(1, 1)
)

plato_figure_layout_a <- patchwork::wrap_plots(
  left_vi,
  q4_hr_plot_ms,
  ncol = 2,
  widths = c(1.0, 1.45)
) +
  patchwork::plot_annotation(
    theme = ggplot2::theme(
      plot.margin = ggplot2::margin(5, 5, 5, 5)
    )
  )

pdf(
  file.path(fig_dir, "figure_plato_workflow.pdf"),
  width = 13,
  height = 8.5
)
print(plato_figure_layout_a)
dev.off()


## ============================================================================
## 10. COMPLETION MESSAGE
## ============================================================================

cat("\n============================================================\n")
cat("PLATO MRCT WATCH analysis completed.\n")
cat("============================================================\n")
cat("Run type:       ", RUN_LABEL, "\n", sep = "")
cat("Results folder: ", results_dir, "\n", sep = "")
cat("Figures folder: ", fig_dir, "\n", sep = "")
cat("Cache folder:   ", cache_dir, "\n", sep = "")
cat("\nKey outputs:\n")
cat("  - cox_regional_consistency.csv\n")
cat("  - workflow_pvalue_summary.csv\n")
cat("  - q2_variable_importance.csv\n")
cat("  - q3_variable_importance.csv\n")
cat("  - q2_q3_rank_comparison.csv\n")
cat("  - hr_by_region_asa.csv\n")
cat("\nKey figures:\n")
cat("  - q2_variable_importance_top15.pdf\n")
cat("  - q3_variable_importance_top15.pdf\n")
cat("  - q4_asa_phi.pdf\n")
cat("  - hr_by_region_asa.pdf\n")
cat("  - q4_asa_hr.pdf\n")
cat("  - figure_plato_workflow.pdf\n")
cat("============================================================\n")
