###############################################################################
## data_analysis/util/watch_functions.R
##
## Model-based treatment effect heterogeneity functions extracted from the
## pre-release WATCH TEH package.
##
## Purpose
## -------
## This file provides the functions required to run the model-based WATCH
## workflow without installing or loading the development version of the TEH
## package.
##
## The implementation is kept as close as possible to the original package
## implementation. Explicit namespace prefixes have been added so that the
## functions can be sourced directly.
##
## Main function
## -------------
## exploreTEH_mb()
##   Runs model-based treatment effect heterogeneity analysis:
##     1. Validates the input data.
##     2. Constructs a prognostic risk score using elastic-net regression.
##     3. Fits a homogeneous treatment-effect model.
##     4. Extracts treatment score residuals, phi.
##     5. Performs a global coin independence test.
##     6. Fits a conditional random forest and computes permutation variable
##        importance.
##
## Supporting functions
## --------------------
## validate_inputs()
##   Validates treatment, outcome, event, and covariate inputs.
##
## check_outcome_type()
##   Checks the declared outcome type and returns the corresponding model family.
##
## riskPred()
##   Constructs a prognostic risk score using glmnet. For survival outcomes,
##   a ridge or elastic-net Cox model is fitted.
##
## teh_varimp()
##   Fits a conditional random forest for phi as a function of the covariates
##   and computes permutation variable importance.
##
## new_TEH()
##   Constructs the returned TEH S3 object.
##
## validate_TEH()
##   Checks the structure and required fields of a TEH object.
##
## is_TEH()
##   Tests whether an object inherits from class TEH.
##
## get_importance_scores()
##   Extracts variable importance scores from a TEH object as a data frame.
##
## s_value()
##   Converts a p-value to a surprise value in bits.
##
## evidence_label()
##   Maps a surprise value to a verbal evidence category.
##
## print.TEH()
##   Prints the global heterogeneity result and leading variable importance
##   scores.
##
## Required packages
## -----------------
## glmnet
## survival
## sandwich
## coin
## party
## permimp
## MASS
##
## Usage
## -----
## source(here::here("data_analysis", "util", "watch_functions.R"))
##
## result <- exploreTEH_mb(
##   data = dat_teh,
##   trt_name = "treatment",
##   y_name = "time_days",
##   event_name = "event",
##   drug_name = "ticagrelor",
##   control_name = "clopidogrel",
##   type = "survival"
## )
###############################################################################

## Check required packages without attaching them to the search path.
.watch_required_packages <- c(
  "glmnet",
  "survival",
  "sandwich",
  "coin",
  "party",
  "permimp",
  "MASS"
)

.watch_missing_packages <- .watch_required_packages[
  !vapply(
    .watch_required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(.watch_missing_packages) > 0L) {
  stop(
    "The following packages are required by watch_functions.R but are not ",
    "installed: ",
    paste(.watch_missing_packages, collapse = ", "),
    call. = FALSE
  )
}

rm(.watch_required_packages, .watch_missing_packages)


## ============================================================================
## 1. TEH OBJECT CONSTRUCTOR AND VALIDATION
## ============================================================================

#' Internal constructor for TEH objects
#'
#' Creates a treatment effect heterogeneity result object containing the
#' score residuals, global test result, and variable importance object.
#'
#' @param phi Numeric vector of treatment score residuals.
#' @param p_value Numeric scalar containing the global heterogeneity p-value.
#' @param method Character string describing the analysis method.
#' @param outcome_type Character outcome type.
#' @param cforest_fit Optional fitted conditional random forest.
#' @param varimp A permimp VarImp object.
#' @param details List of analysis parameters.
#' @param call Matched function call.
#' @param treatment Treatment variable name.
#' @param outcome Outcome variable name.
#' @param event Event-indicator variable name.
#' @param treatment_labels Named treatment labels.
#'
#' @return An object of class TEH.
new_TEH <- function(
    phi,
    p_value,
    method,
    outcome_type = c("continuous", "binary", "count", "survival"),
    cforest_fit = NULL,
    varimp,
    details = list(),
    call = match.call(),
    treatment = NULL,
    outcome = NULL,
    event = NULL,
    treatment_labels = NULL
) {
  outcome_type <- match.arg(outcome_type)
  
  obj <- list(
    phi = phi,
    p_value = p_value,
    method = method,
    outcome_type = outcome_type,
    cforest_fit = cforest_fit,
    varimp = varimp,
    details = details,
    call = call,
    treatment = treatment,
    outcome = outcome,
    event = event,
    treatment_labels = treatment_labels
  )
  
  class(obj) <- c("TEH", "list")
  
  validate_TEH(obj)
  obj
}


#' Validate a TEH object
#'
#' Checks required fields, classes, and dimensions of a TEH result object.
#'
#' @param x Object to validate.
#'
#' @return Invisibly returns TRUE if the object is valid.
validate_TEH <- function(x) {
  if (!inherits(x, "TEH")) {
    stop("Object is not of class 'TEH'.", call. = FALSE)
  }
  
  required <- c(
    "phi",
    "p_value",
    "varimp",
    "method",
    "outcome_type",
    "cforest_fit",
    "details",
    "call"
  )
  
  missing_fields <- setdiff(required, names(x))
  
  if (length(missing_fields) > 0L) {
    stop(
      "TEH object missing required fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }
  
  if (!is.numeric(x$phi) || length(x$phi) == 0L) {
    stop(
      "'phi' must be a non-empty numeric vector.",
      call. = FALSE
    )
  }
  
  if (
    !is.numeric(x$p_value) ||
    length(x$p_value) != 1L ||
    !is.finite(x$p_value) ||
    x$p_value < 0 ||
    x$p_value > 1
  ) {
    stop(
      "'p_value' must be a single finite numeric value between 0 and 1.",
      call. = FALSE
    )
  }
  
  if (!inherits(x$varimp, "VarImp")) {
    stop(
      "'varimp' must inherit class 'VarImp' from permimp.",
      call. = FALSE
    )
  }
  
  if (!is.character(x$method) || length(x$method) == 0L) {
    stop(
      "'method' must be a non-empty character string.",
      call. = FALSE
    )
  }
  
  if (
    !x$outcome_type %in%
    c("continuous", "binary", "count", "survival")
  ) {
    stop(
      "'outcome_type' must be one of continuous, binary, count, or survival.",
      call. = FALSE
    )
  }
  
  if (!is.list(x$details)) {
    stop(
      "'details' must be a list.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}


#' Test whether an object is a TEH object
#'
#' @param x Object to test.
#'
#' @return Logical scalar.
is_TEH <- function(x) {
  inherits(x, "TEH")
}


## ============================================================================
## 2. INPUT VALIDATION
## ============================================================================

#' Validate inputs for model-based TEH exploration
#'
#' Validates the input data, treatment labels, outcome, event indicator, and
#' covariate types. It also recodes the investigational treatment to 1 and the
#' control treatment to 0.
#'
#' @param data Input data frame.
#' @param trt_name Treatment column name.
#' @param y_name Outcome column name.
#' @param event_name Event-indicator column name for survival outcomes.
#' @param type Outcome type.
#' @param drug_name Investigational treatment label.
#' @param control_name Control treatment label.
#' @param allow_missing_y Whether missing outcomes are permitted.
#'
#' @return A list containing the filtered data, X, trt, y, and event.
validate_inputs <- function(
    data,
    trt_name,
    y_name,
    event_name = NULL,
    type,
    drug_name = "1",
    control_name = "0",
    allow_missing_y = FALSE
) {
  if (!is.data.frame(data)) {
    stop("data must be a data.frame")
  }
  
  if (!(trt_name %in% colnames(data))) {
    stop(
      trt_name,
      " column not identified in the data provided"
    )
  }
  
  if (!(y_name %in% colnames(data))) {
    stop(
      y_name,
      " column not identified in the data provided"
    )
  }
  
  if (anyNA(data[[trt_name]])) {
    stop("trt must not contain missing values")
  }
  
  trt_category <- table(data[[trt_name]])
  
  if (length(trt_category) == 1L) {
    stop(
      "Two treatment levels are needed; only one was identified."
    )
  }
  
  trt_levels <- names(trt_category)
  
  if (!all(c(drug_name, control_name) %in% trt_levels)) {
    stop(
      "Please specify drug and control names correctly.\n",
      "Treatment levels identified are: ",
      paste(trt_levels, collapse = ", "),
      ".\nDrug name provided: ",
      drug_name,
      ".\nControl name provided: ",
      control_name,
      "."
    )
  }
  
  if (length(trt_category) > 2L) {
    warning(
      paste(
        "More than two treatment levels identified; filtering to include",
        drug_name,
        "and",
        control_name
      ),
      call. = FALSE
    )
  }
  
  data <- data[
    data[[trt_name]] %in% c(drug_name, control_name),
    ,
    drop = FALSE
  ]
  
  excluded_names <- c(trt_name, y_name, event_name)
  X <- data[
    ,
    !(names(data) %in% excluded_names),
    drop = FALSE
  ]
  
  y <- data[[y_name]]
  
  y <- as.numeric(
    if (is.data.frame(y)) {
      y[[1]]
    } else if (is.matrix(y)) {
      y[, 1]
    } else {
      y
    }
  )
  
  trt_raw <- data[[trt_name]]
  
  trt <- as.numeric(
    ifelse(
      trt_raw == drug_name,
      1,
      ifelse(trt_raw == control_name, 0, NA)
    )
  )
  
  if (!is.null(event_name) && type == "continuous") {
    stop(
      "event_name was supplied but type = 'continuous'; ",
      "use type = 'survival'."
    )
  }
  
  if (type == "survival") {
    if (is.null(event_name)) {
      stop(
        "event_name must be provided for survival outcomes"
      )
    }
    
    if (!(event_name %in% colnames(data))) {
      stop(
        event_name,
        " column not identified in the data provided"
      )
    }
    
    event <- data[[event_name]]
    
    event <- as.numeric(
      if (is.data.frame(event)) {
        event[[1]]
      } else if (is.matrix(event)) {
        event[, 1]
      } else {
        event
      }
    )
    
    if (anyNA(event)) {
      stop("event must not contain missing values")
    }
  } else {
    event <- NULL
  }
  
  if (anyNA(X)) {
    stop(
      "Covariates input must not contain missing values"
    )
  }
  
  is_valid <- vapply(
    X,
    function(column) {
      is.numeric(column) || is.factor(column)
    },
    FUN.VALUE = logical(1)
  )
  
  if (!all(is_valid)) {
    bad_names <- names(X)[!is_valid]
    max_show <- 5L
    
    if (length(bad_names) > max_show) {
      bad_list <- paste0(
        paste(head(bad_names, max_show), collapse = ", "),
        ", ... and ",
        length(bad_names) - max_show,
        " more"
      )
    } else {
      bad_list <- paste(bad_names, collapse = ", ")
    }
    
    stop(
      "Input validation failed: covariate columns must be numeric or factor.\n",
      "Unsupported columns: ",
      bad_list,
      call. = FALSE
    )
  }
  
  if (!allow_missing_y && anyNA(y)) {
    stop(
      "y_name column must not contain missing values"
    )
  }
  
  if (!is.numeric(y)) {
    stop(
      "y_name column must be numeric"
    )
  }
  
  if (length(trt) != nrow(X)) {
    stop(
      "Treatment and covariate data must have the same number of rows."
    )
  }
  
  if (length(y) != nrow(X)) {
    stop(
      "Outcome and covariate data must have the same number of rows."
    )
  }
  
  if (type == "survival" && length(event) != nrow(X)) {
    stop(
      "Event indicator and covariate data must have the same number of rows."
    )
  }
  
  list(
    data = data,
    X = X,
    trt = trt,
    y = y,
    event = event
  )
}


## ============================================================================
## 3. OUTCOME-TYPE VALIDATION
## ============================================================================

#' Check and standardize outcome type
#'
#' Validates the outcome according to the declared outcome type and returns the
#' corresponding GLM family where appropriate.
#'
#' @param type One of continuous, binary, count, or survival.
#' @param y Numeric outcome.
#' @param event Optional survival event indicator.
#' @param allow_missing Whether missing outcome values are permitted.
#'
#' @return A list containing type and family.
check_outcome_type <- function(
    type,
    y,
    event = NULL,
    allow_missing = FALSE
) {
  type <- match.arg(
    type,
    choices = c(
      "continuous",
      "binary",
      "count",
      "survival"
    )
  )
  
  if (!allow_missing && anyNA(y)) {
    stop(
      "y must not contain missing values for this method"
    )
  }
  
  if (type == "binary") {
    if (!all(y[!is.na(y)] %in% c(0, 1))) {
      stop(
        "For type = 'binary', y must be coded strictly as 0 or 1."
      )
    }
    
    family_obj <- stats::binomial()
  } else if (type == "count") {
    if (any(y < 0)) {
      stop(
        "For type = 'count', y must be non-negative."
      )
    }
    
    if (!all(y == floor(y))) {
      stop(
        "For type = 'count', y must be integer-valued."
      )
    }
    
    family_obj <- NULL
  } else if (type == "survival") {
    if (any(y < 0)) {
      stop(
        "For type = 'survival', time y must be non-negative."
      )
    }
    
    if (!is.null(event) && !all(event %in% c(0, 1))) {
      stop(
        "For type = 'survival', event must be coded as 0 or 1."
      )
    }
    
    family_obj <- NULL
  } else {
    if (all(y %in% c(0, 1))) {
      stop(
        "y looks binary but type = 'continuous'; use type = 'binary'."
      )
    }
    
    if (all(y >= 0) && all(y == floor(y))) {
      stop(
        "y is non-negative and integer-valued but type = 'continuous'; ",
        "use type = 'count'."
      )
    }
    
    family_obj <- stats::gaussian()
  }
  
  list(
    type = type,
    family = family_obj
  )
}


## ============================================================================
## 4. PROGNOSTIC RISK SCORE
## ============================================================================

#' Build a prognostic risk score using elastic net
#'
#' Encodes categorical covariates using a model matrix and fits a cross-validated
#' elastic-net model. For survival outcomes, a Cox model is fitted.
#'
#' @param X Data frame of covariates.
#' @param Y Numeric outcome or survival time.
#' @param event Survival event indicator.
#' @param outcome_info Result from check_outcome_type().
#' @param alpha Elastic-net mixing parameter. Zero corresponds to ridge.
#' @param relax Whether to use a relaxed elastic-net fit.
#' @param over_dispersion Count-outcome overdispersion setting.
#'
#' @return Numeric prognostic risk score on the link scale.
riskPred <- function(
    X,
    Y,
    event = NULL,
    outcome_info,
    alpha = 0,
    relax = FALSE,
    over_dispersion = NULL
) {
  factor_cols <- names(X)[
    vapply(X, is.factor, FUN.VALUE = logical(1))
  ]
  
  if (length(factor_cols) > 0L) {
    contrast_list <- lapply(
      factor_cols,
      function(column_name) {
        stats::contrasts(
          X[[column_name]],
          contrasts = FALSE
        )
      }
    )
    
    names(contrast_list) <- factor_cols
  } else {
    contrast_list <- NULL
  }
  
  X_mm <- stats::model.matrix(
    ~ 0 + .,
    data = X,
    contrasts.arg = contrast_list
  )
  
  y_type <- outcome_info$type
  
  if (y_type == "survival") {
    if (is.null(event)) {
      stop(
        "event must be provided for survival outcomes"
      )
    }
    
    cv_fit <- glmnet::cv.glmnet(
      x = X_mm,
      y = survival::Surv(Y, event),
      family = "cox",
      alpha = alpha,
      relax = relax
    )
  } else if (y_type == "count") {
    if (isFALSE(over_dispersion)) {
      message(
        "Using Poisson regression for risk modelling."
      )
      
      cv_fit <- glmnet::cv.glmnet(
        x = X_mm,
        y = Y,
        family = "poisson",
        alpha = alpha,
        relax = relax
      )
    } else {
      if (is.null(over_dispersion)) {
        message(
          "over_dispersion unspecified; using Negative Binomial ",
          "regression for risk modelling."
        )
      } else {
        message(
          "Using Negative Binomial regression for risk modelling."
        )
      }
      
      nb_data <- data.frame(
        Y = Y,
        X_mm,
        check.names = FALSE
      )
      
      theta_hat <- MASS::glm.nb(
        Y ~ .,
        data = nb_data
      )$theta
      
      family_nb <- MASS::negative.binomial(
        theta = theta_hat
      )
      
      cv_fit <- glmnet::cv.glmnet(
        x = X_mm,
        y = Y,
        family = family_nb,
        alpha = alpha,
        relax = relax
      )
    }
  } else {
    cv_fit <- glmnet::cv.glmnet(
      x = X_mm,
      y = Y,
      family = outcome_info$family,
      alpha = alpha,
      relax = relax
    )
  }
  
  if (isTRUE(relax)) {
    risk <- stats::predict(
      cv_fit,
      newx = X_mm,
      s = "lambda.min",
      type = "link",
      gamma = "gamma.min"
    )
  } else {
    risk <- stats::predict(
      cv_fit,
      newx = X_mm,
      s = "lambda.min",
      type = "link"
    )
  }
  
  as.numeric(risk)
}


## ============================================================================
## 5. CONDITIONAL RANDOM FOREST VARIABLE IMPORTANCE
## ============================================================================

#' Estimate permutation variable importance from TEH score residuals
#'
#' Fits a conditional random forest with phi as the response and computes
#' unconditional permutation variable importance.
#'
#' @param X Covariate data frame.
#' @param phi Numeric score residuals.
#' @param mtry Number of variables considered at each split.
#' @param ntree Number of trees.
#' @param nperm Number of variable-importance permutations.
#' @param verbose Whether to display a progress bar.
#' @param include_cforest Whether to retain the fitted forest.
#' @param cl Optional cluster supplied to permimp.
#'
#' @return A list containing varimp and optionally cforest_fit.
teh_varimp <- function(
    X,
    phi,
    mtry = 5,
    ntree = 500,
    nperm = 10,
    verbose = FALSE,
    include_cforest = FALSE,
    cl = NULL
) {
  dat <- data.frame(
    X,
    y = phi,
    check.names = FALSE
  )
  
  control_cforest <- party::cforest_unbiased(
    mtry = mtry,
    ntree = ntree
  )
  
  cforest_fit <- party::cforest(
    y ~ .,
    data = dat,
    control = control_cforest
  )
  
  if (is.null(cl)) {
    cf_vi <- permimp::permimp(
      cforest_fit,
      conditional = FALSE,
      nperm = nperm,
      progressBar = verbose
    )
  } else {
    cf_vi <- permimp::permimp(
      cforest_fit,
      conditional = FALSE,
      nperm = nperm,
      progressBar = verbose,
      cl = cl
    )
  }
  
  list(
    varimp = cf_vi,
    cforest_fit = if (include_cforest) {
      cforest_fit
    } else {
      NULL
    }
  )
}


## ============================================================================
## 6. VARIABLE-IMPORTANCE EXTRACTION
## ============================================================================

#' Extract variable importance scores from a TEH object
#'
#' @param x TEH object.
#'
#' @return Data frame with variable and vi columns.
get_importance_scores <- function(x) {
  if (!inherits(x, "TEH")) {
    stop(
      "'x' must be a TEH object",
      call. = FALSE
    )
  }
  
  vi_values <- x$varimp$values
  
  data.frame(
    variable = names(vi_values),
    vi = as.numeric(vi_values),
    stringsAsFactors = FALSE
  )
}


## ============================================================================
## 7. SURPRISE VALUE AND EVIDENCE LABEL
## ============================================================================

#' Compute a surprise value in bits
#'
#' Computes S = -log2(p), with clamping to ensure a finite value.
#'
#' @param p Numeric scalar p-value.
#'
#' @return Numeric surprise value.
s_value <- function(p) {
  if (!is.numeric(p) || length(p) != 1L) {
    stop(
      "'p' must be a single numeric value.",
      call. = FALSE
    )
  }
  
  p <- max(
    min(p, 1),
    .Machine$double.eps
  )
  
  -log2(p)
}


#' Map a surprise value to an evidence category
#'
#' @param s Numeric surprise value.
#'
#' @return Character evidence label.
evidence_label <- function(s) {
  if (
    !is.numeric(s) ||
    length(s) != 1L ||
    !is.finite(s)
  ) {
    return(NA_character_)
  }
  
  if (s <= 2) {
    "low"
  } else if (s <= 4) {
    "moderate"
  } else if (s <= 7) {
    "noteworthy"
  } else if (s <= 10) {
    "strong"
  } else {
    "very strong"
  }
}


## ============================================================================
## 8. PRINT METHOD
## ============================================================================

#' Print method for a TEH object
#'
#' Prints the global heterogeneity p-value, surprise value, evidence category,
#' and leading variable importance scores.
#'
#' @param x TEH object.
#' @param top_n Number of top variables to display.
#' @param digits_s Decimal places for the surprise value.
#' @param digits_vi Significant digits for variable importance.
#' @param ... Additional print arguments.
#'
#' @return The TEH object invisibly.
print.TEH <- function(
    x,
    top_n = 5,
    digits_s = 1,
    digits_vi = 3,
    ...
) {
  stopifnot(inherits(x, "TEH"))
  
  cat("TEH Exploration\n")
  
  cat("\nGlobal heterogeneity test:\n")
  
  p <- x$p_value
  s <- s_value(p)
  evidence <- evidence_label(s)
  
  s_formatted <- formatC(
    s,
    format = "f",
    digits = digits_s
  )
  
  p_formatted <- format.pval(
    p,
    digits = 4,
    eps = 1e-4
  )
  
  cat(
    sprintf(
      "  P-value: %s (Surprise value: %s bits)\n",
      p_formatted,
      s_formatted
    )
  )
  
  cat(
    sprintf(
      "  Evidence against homogeneity: %s\n",
      evidence
    )
  )
  
  cat(
    sprintf(
      "\nVariable Importance (top %s variables):\n",
      top_n
    )
  )
  
  importance <- get_importance_scores(x)
  
  importance <- importance[
    order(importance$vi, decreasing = TRUE),
    ,
    drop = FALSE
  ]
  
  n_show <- min(
    top_n,
    nrow(importance)
  )
  
  print_data <- importance[
    seq_len(n_show),
    ,
    drop = FALSE
  ]
  
  print_data$vi <- formatC(
    print_data$vi,
    format = "g",
    digits = digits_vi
  )
  
  print(
    print_data,
    row.names = FALSE,
    right = FALSE
  )
  
  invisible(x)
}


## ============================================================================
## 9. MODEL-BASED TEH EXPLORATION
## ============================================================================

#' Explore treatment effect heterogeneity using a model-based approach
#'
#' The function first creates a prognostic risk score. A homogeneous treatment
#' effect model is then fitted using the risk score and a centered treatment
#' indicator. Treatment score residuals are extracted and assessed for global
#' heterogeneity and variable importance.
#'
#' @param data Input data frame. All columns other than treatment, outcome, and
#'   event are treated as covariates.
#' @param trt_name Treatment column name.
#' @param y_name Outcome or survival-time column name.
#' @param event_name Event-indicator column for survival outcomes.
#' @param drug_name Investigational treatment label.
#' @param control_name Control treatment label.
#' @param type One of continuous, binary, count, or survival.
#' @param over_dispersion Count-outcome overdispersion setting.
#' @param alpha Elastic-net mixing parameter.
#' @param relax Whether to use relaxed elastic net.
#' @param teststat coin test statistic, usually maximum.
#' @param mtry Conditional random forest mtry.
#' @param ntree Number of trees.
#' @param nperm Number of variable-importance permutations.
#' @param verbose Whether to display progress bars.
#' @param include_cforest Whether to retain the fitted forest.
#' @param cl Optional cluster for permimp.
#'
#' @return A TEH object.
exploreTEH_mb <- function(
    data,
    trt_name,
    y_name,
    event_name = NULL,
    drug_name = "1",
    control_name = "0",
    type,
    over_dispersion = NULL,
    alpha = 0,
    relax = FALSE,
    teststat = "maximum",
    mtry = 5,
    ntree = 500,
    nperm = 10,
    verbose = FALSE,
    include_cforest = FALSE,
    cl = NULL
) {
  if (missing(type) || is.null(type)) {
    stop(
      "Argument 'type' is required and must be one of: ",
      "'continuous', 'binary', 'count', or 'survival'."
    )
  }
  
  type <- match.arg(
    type,
    choices = c(
      "continuous",
      "binary",
      "count",
      "survival"
    )
  )
  
  call <- match.call()
  
  input_results <- validate_inputs(
    data = data,
    trt_name = trt_name,
    y_name = y_name,
    event_name = event_name,
    drug_name = drug_name,
    type = type,
    control_name = control_name,
    allow_missing_y = FALSE
  )
  
  data <- input_results$data
  X <- input_results$X
  y <- input_results$y
  trt <- input_results$trt
  event <- input_results$event
  
  if (type == "survival") {
    if (
      !is.numeric(event) ||
      !all(event %in% c(0, 1))
    ) {
      stop(
        "For type = 'survival', event must be coded as numeric 0 or 1."
      )
    }
    
    if (any(y <= 0)) {
      stop(
        "For type = 'survival', all event and censoring times must be positive."
      )
    }
  }
  
  if (type != "count" && !is.null(over_dispersion)) {
    message(
      "Argument 'over_dispersion' is only relevant for count outcomes ",
      "and will be ignored for type = '",
      type,
      "'."
    )
    
    over_dispersion <- NULL
  }
  
  if (type == "count" && !is.null(over_dispersion)) {
    if (
      !is.logical(over_dispersion) ||
      length(over_dispersion) != 1L
    ) {
      stop(
        "For count outcomes, 'over_dispersion' must be TRUE, FALSE, or NULL."
      )
    }
  }
  
  outcome_info <- check_outcome_type(
    type = type,
    y = y,
    event = event,
    allow_missing = FALSE
  )
  
  risk <- riskPred(
    X = X,
    Y = y,
    event = event,
    outcome_info = outcome_info,
    alpha = alpha,
    relax = relax,
    over_dispersion = if (type == "count") {
      over_dispersion
    } else {
      NULL
    }
  )
  
  trt_centered <- trt - mean(trt)
  
  if (type == "survival") {
    fit_homogeneous <- survival::coxph(
      survival::Surv(y, event) ~ trt_centered + risk
    )
    
    score_matrix <- sandwich::estfun(
      fit_homogeneous
    )
    
    score_residual <- as.numeric(
      score_matrix[, "trt_centered"]
    )
  } else if (type == "count") {
    if (isFALSE(over_dispersion)) {
      message(
        "Using Poisson regression for the homogeneous treatment-effect model."
      )
      
      fit_homogeneous <- stats::glm(
        y ~ trt_centered + risk,
        family = stats::poisson()
      )
    } else {
      if (is.null(over_dispersion)) {
        message(
          "over_dispersion unspecified; using Negative Binomial regression ",
          "for the homogeneous treatment-effect model."
        )
      } else {
        message(
          "Using Negative Binomial regression for the homogeneous ",
          "treatment-effect model."
        )
      }
      
      fit_homogeneous <- MASS::glm.nb(
        y ~ trt_centered + risk
      )
    }
    
    score_matrix <- sandwich::estfun(
      fit_homogeneous
    )
    
    score_residual <- as.numeric(
      score_matrix[, "trt_centered"]
    )
  } else {
    fit_homogeneous <- stats::glm(
      y ~ trt_centered + risk,
      family = outcome_info$family
    )
    
    score_matrix <- sandwich::estfun(
      fit_homogeneous
    )
    
    score_residual <- as.numeric(
      score_matrix[, "trt_centered"]
    )
  }
  
  global_test_data <- data.frame(
    X,
    y = score_residual,
    check.names = FALSE
  )
  
  global_test <- coin::independence_test(
    y ~ .,
    data = global_test_data,
    teststat = teststat,
    distribution = coin::asymptotic()
  )
  
  p_value <- as.numeric(
    coin::pvalue(global_test)
  )
  
  variable_importance <- teh_varimp(
    X = X,
    phi = score_residual,
    mtry = mtry,
    ntree = ntree,
    nperm = nperm,
    verbose = verbose,
    include_cforest = include_cforest,
    cl = cl
  )
  
  details <- list(
    alpha = alpha,
    relax = relax,
    teststat = teststat,
    test_distribution = "asymptotic",
    test_statistic = as.numeric(
      coin::statistic(global_test)
    ),
    mtry = mtry,
    ntree = ntree,
    nperm = nperm,
    over_dispersion = over_dispersion
  )
  
  new_TEH(
    phi = score_residual,
    p_value = p_value,
    varimp = variable_importance$varimp,
    method = "model-based",
    outcome_type = type,
    cforest_fit = variable_importance$cforest_fit,
    details = details,
    call = call,
    treatment = trt_name,
    outcome = y_name,
    event = event_name,
    treatment_labels = c(
      "0" = control_name,
      "1" = drug_name
    )
  )
}