###############################################################################
## data_analysis/util/utils.R
##
## Utility functions for the PLATO MRCT WATCH analysis.
##
## Purpose
## -------
## This file contains analysis helpers that are not part of the extracted
## WATCH model-based TEH implementation.
##
## Included functions
## ------------------
## run_or_load()
##   Runs an analysis step and saves the result as an RDS checkpoint, or loads
##   an existing checkpoint when rerunning the analysis.
##
## get_coin_scalar_stat()
##   Extracts one scalar statistic from a coin independence-test object.
##   For a maximum test, the largest absolute standardized statistic is used.
##
## assess_regional_imbalance()
##   Implements the Region versus X analysis:
##     1. Global coin independence test using the maximum statistic and
##        asymptotic distribution.
##     2. Conditional random forest.
##     3. Permutation variable importance.
##
## build_workflow_pvalue_summary()
##   Combines the p-values and maximum statistics from:
##     - phi versus Region
##     - Region versus X
##     - phi versus Region plus X
##
## plot_variable_importance()
##   Produces a horizontal variable-importance bar plot.
##
## plot_q4_display_categorical()
##   Produces the main Q4 display for a categorical candidate effect modifier:
##     - Upper panel: mean phi and 95 percent CI by covariate level and region.
##     - Lower panel: covariate-level proportions within each region.
##
## estimate_hr_by_subgroup()
##   Estimates treatment hazard ratios separately within each combination of
##   region and categorical candidate effect-modifier level.
##
## plot_hr_by_subgroup()
##   Produces a forest plot of subgroup-specific treatment hazard ratios.
##
## Required packages
## -----------------
## dplyr
## ggplot2
## coin
## party
## permimp
## survival
## patchwork
##
## Usage
## -----
## source(
##   here::here(
##     "data_analysis",
##     "util",
##     "utils.R"
##   )
## )
###############################################################################

.utils_required_packages <- c(
  "dplyr",
  "ggplot2",
  "coin",
  "party",
  "permimp",
  "survival",
  "patchwork"
)

.utils_missing_packages <- .utils_required_packages[
  !vapply(
    .utils_required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(.utils_missing_packages) > 0L) {
  stop(
    "The following packages are required by utils.R but are not installed: ",
    paste(.utils_missing_packages, collapse = ", "),
    call. = FALSE
  )
}

rm(.utils_required_packages, .utils_missing_packages)


## ============================================================================
## 1. CHECKPOINT MANAGEMENT
## ============================================================================

#' Run an analysis step or load its existing checkpoint
#'
#' @param path Path to the RDS checkpoint.
#' @param fun Function with no arguments that runs the analysis step.
#' @param force Logical; rerun the step even if its checkpoint exists.
#'
#' @return The saved or newly computed result.
run_or_load <- function(path, fun, force = FALSE) {
  if (!is.character(path) || length(path) != 1L) {
    stop(
      "'path' must be a single character string.",
      call. = FALSE
    )
  }
  
  if (!is.function(fun)) {
    stop(
      "'fun' must be a function.",
      call. = FALSE
    )
  }
  
  if (file.exists(path) && !isTRUE(force)) {
    message("Loading existing result: ", path)
    return(readRDS(path))
  }
  
  output_directory <- dirname(path)
  
  if (!dir.exists(output_directory)) {
    dir.create(
      output_directory,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  message("Running step: ", basename(path))
  
  result <- fun()
  
  ## Save to a temporary file first to reduce the risk of leaving a corrupted
  ## checkpoint if the R session stops during saveRDS().
  temporary_path <- paste0(
    path,
    ".tmp"
  )
  
  saveRDS(
    result,
    file = temporary_path
  )
  
  if (file.exists(path)) {
    file.remove(path)
  }
  
  renamed <- file.rename(
    from = temporary_path,
    to = path
  )
  
  if (!renamed) {
    file.copy(
      from = temporary_path,
      to = path,
      overwrite = TRUE
    )
    
    file.remove(temporary_path)
  }
  
  message("Saved result: ", path)
  
  result
}


## ============================================================================
## 2. COIN TEST HELPERS
## ============================================================================

#' Extract a scalar statistic from a coin test object
#'
#' For multivariate statistics, returns the component with the largest absolute
#' value. This is appropriate for summarizing a maximum-type test.
#'
#' @param test A coin independence-test object.
#' @param type Statistic type passed to coin::statistic(). Use "standardized"
#'   for the maximum test statistic.
#'
#' @return Numeric scalar.
get_coin_scalar_stat <- function(
    test,
    type = "standardized"
) {
  statistic_value <- if (is.null(type)) {
    coin::statistic(test)
  } else {
    coin::statistic(
      test,
      type = type
    )
  }
  
  statistic_value <- as.numeric(
    statistic_value
  )
  
  statistic_value <- statistic_value[
    is.finite(statistic_value)
  ]
  
  if (length(statistic_value) == 0L) {
    return(NA_real_)
  }
  
  statistic_value[
    which.max(abs(statistic_value))
  ]
}


#' Create a compact summary of a coin maximum test
#'
#' @param test A coin independence-test object.
#' @param question Character description of the analysis question.
#'
#' @return One-row data frame.
summarise_coin_test <- function(
    test,
    question
) {
  data.frame(
    question = question,
    test = "coin independence test",
    statistic_type = "maximum",
    statistic = get_coin_scalar_stat(
      test,
      type = "standardized"
    ),
    p_value = as.numeric(
      coin::pvalue(test)
    ),
    stringsAsFactors = FALSE
  )
}


## ============================================================================
## 3. REGIONAL IMBALANCE ASSESSMENT
## ============================================================================

#' Assess associations between region and covariates
#'
#' Performs a global coin independence test of Region versus X and estimates
#' permutation variable importance from a conditional random forest.
#'
#' @param X Data frame containing candidate covariates.
#' @param Region Factor defining the regions.
#' @param mtry Number of variables considered at each forest split.
#' @param ntree Number of trees.
#' @param nperm Number of variable-importance permutations.
#' @param verbose Logical; show the permimp progress bar.
#' @param cl Optional cluster passed to permimp.
#' @param include_cforest Logical; retain the fitted forest in the result.
#'
#' @return A list containing the coin test, p-value, maximum statistic,
#'   variable-importance object and table, forest settings, and optionally the
#'   fitted forest.
assess_regional_imbalance <- function(
    X,
    Region,
    mtry = 5,
    ntree = 100,
    nperm = 5,
    verbose = FALSE,
    cl = NULL,
    include_cforest = FALSE
) {
  if (!is.data.frame(X)) {
    X <- as.data.frame(X)
  }
  
  if (nrow(X) != length(Region)) {
    stop(
      "X and Region must contain the same number of observations.",
      call. = FALSE
    )
  }
  
  if (anyNA(X)) {
    stop(
      "X must not contain missing values.",
      call. = FALSE
    )
  }
  
  if (anyNA(Region)) {
    stop(
      "Region must not contain missing values.",
      call. = FALSE
    )
  }
  
  Region <- factor(Region)
  
  if (nlevels(Region) < 2L) {
    stop(
      "Region must contain at least two levels.",
      call. = FALSE
    )
  }
  
  test_data <- data.frame(
    Region = Region,
    X,
    check.names = FALSE
  )
  
  ## Global test of association between Region and the joint covariate set.
  test_maximum <- coin::independence_test(
    Region ~ .,
    data = test_data,
    teststat = "maximum",
    distribution = coin::asymptotic()
  )
  
  ## Conditional random forest for regional imbalance ranking.
  forest_control <- party::cforest_unbiased(
    mtry = mtry,
    ntree = ntree
  )
  
  forest_fit <- party::cforest(
    Region ~ .,
    data = test_data,
    control = forest_control
  )
  
  if (is.null(cl)) {
    variable_importance <- permimp::permimp(
      forest_fit,
      conditional = FALSE,
      nperm = nperm,
      progressBar = verbose
    )
  } else {
    variable_importance <- permimp::permimp(
      forest_fit,
      conditional = FALSE,
      nperm = nperm,
      progressBar = verbose,
      cl = cl
    )
  }
  
  importance_values <- variable_importance$values
  
  importance_table <- data.frame(
    covariate = names(importance_values),
    vi_region = as.numeric(importance_values),
    stringsAsFactors = FALSE
  )
  
  importance_table <- importance_table[
    order(
      importance_table$vi_region,
      decreasing = TRUE,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  rownames(importance_table) <- NULL
  
  list(
    test = test_maximum,
    region_p_max = as.numeric(
      coin::pvalue(test_maximum)
    ),
    statistic = get_coin_scalar_stat(
      test_maximum,
      type = "standardized"
    ),
    statistic_type = "maximum",
    test_distribution = "asymptotic",
    varimp = variable_importance,
    vi_region = importance_values,
    vi_table = importance_table,
    mtry = mtry,
    ntree = ntree,
    nperm = nperm,
    cforest_fit = if (include_cforest) {
      forest_fit
    } else {
      NULL
    }
  )
}


## ============================================================================
## 4. WORKFLOW P-VALUE SUMMARY
## ============================================================================

#' Build a combined p-value summary for the WATCH workflow
#'
#' Combines the three global tests from:
#'   1. phi versus Region.
#'   2. Region versus X.
#'   3. phi versus Region plus X.
#'
#' The Cox consistency analysis is intentionally excluded because it reproduces
#' published PLATO results and is not one of these three WATCH global tests.
#'
#' @param step2 Result object for phi versus Region. It should contain either a
#'   coin test in $test or scalar values in $p_value and $statistic.
#' @param step3 Result from assess_regional_imbalance().
#' @param step4 WATCH TEH object or a list containing global_p_value and
#'   statistic.
#'
#' @return Three-row data frame.
build_workflow_pvalue_summary <- function(
    step2,
    step3,
    step4
) {
  ## Step 2: phi versus Region
  if (!is.null(step2$test)) {
    step2_row <- summarise_coin_test(
      step2$test,
      question = "phi ~ Region"
    )
  } else {
    step2_row <- data.frame(
      question = "phi ~ Region",
      test = "coin independence test",
      statistic_type = "maximum",
      statistic = as.numeric(step2$statistic),
      p_value = as.numeric(step2$p_value),
      stringsAsFactors = FALSE
    )
  }
  
  ## Step 3: Region versus X
  if (!is.null(step3$test)) {
    step3_row <- summarise_coin_test(
      step3$test,
      question = "Region ~ X"
    )
  } else {
    step3_row <- data.frame(
      question = "Region ~ X",
      test = "coin independence test",
      statistic_type = "maximum",
      statistic = as.numeric(step3$statistic),
      p_value = as.numeric(step3$region_p_max),
      stringsAsFactors = FALSE
    )
  }
  
  ## Step 4: phi versus Region plus X
  if (inherits(step4, "TEH")) {
    step4_p_value <- step4$p_value
    
    step4_statistic <- step4$details$test_statistic
    
    if (is.null(step4_statistic)) {
      step4_statistic <- NA_real_
    }
    
    step4_statistic <- as.numeric(
      step4_statistic
    )
    
    if (length(step4_statistic) > 1L) {
      step4_statistic <- step4_statistic[
        which.max(abs(step4_statistic))
      ]
    }
  } else {
    step4_p_value <- if (!is.null(step4$global_p_value)) {
      step4$global_p_value
    } else {
      step4$p_value
    }
    
    step4_statistic <- if (!is.null(step4$statistic)) {
      step4$statistic
    } else if (!is.null(step4$details$test_statistic)) {
      step4$details$test_statistic
    } else {
      NA_real_
    }
    
    step4_statistic <- as.numeric(
      step4_statistic
    )
    
    if (length(step4_statistic) > 1L) {
      step4_statistic <- step4_statistic[
        which.max(abs(step4_statistic))
      ]
    }
  }
  
  step4_row <- data.frame(
    question = "phi ~ Region + X",
    test = "coin independence test",
    statistic_type = "maximum",
    statistic = step4_statistic,
    p_value = as.numeric(step4_p_value),
    stringsAsFactors = FALSE
  )
  
  result <- rbind(
    step2_row,
    step3_row,
    step4_row
  )
  
  result$step <- c(
    "Q1",
    "Q2",
    "Q3"
  )
  
  result <- result[
    ,
    c(
      "step",
      "question",
      "test",
      "statistic_type",
      "statistic",
      "p_value"
    ),
    drop = FALSE
  ]
  
  rownames(result) <- NULL
  
  result
}


## ============================================================================
## 5. VARIABLE-IMPORTANCE PLOT
## ============================================================================

#' Plot leading variable-importance values
#'
#' @param vi_table Data frame containing a covariate column and a variable-
#'   importance column.
#' @param value_col Name of the variable-importance column.
#' @param top_n Number of variables to display.
#' @param highlight Character vector of variables to highlight.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param y_label Axis label for variable importance.
#' @param highlight_colour Colour used for highlighted variables.
#' @param other_colour Colour used for other variables.
#'
#' @return ggplot object.
plot_variable_importance <- function(
    vi_table,
    value_col,
    top_n = 15,
    highlight = character(),
    title = NULL,
    subtitle = NULL,
    y_label = "Variable importance",
    highlight_colour = "#1B7837",
    other_colour = "grey70"
) {
  required_columns <- c(
    "covariate",
    value_col
  )
  
  if (!all(required_columns %in% names(vi_table))) {
    stop(
      "vi_table must contain: ",
      paste(required_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  plot_data <- vi_table[
    order(
      vi_table[[value_col]],
      decreasing = TRUE,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  plot_data <- utils::head(
    plot_data,
    top_n
  )
  
  plot_data$covariate_plot <- factor(
    plot_data$covariate,
    levels = rev(plot_data$covariate)
  )
  
  plot_data$highlighted <- plot_data$covariate %in% highlight
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = covariate_plot,
      y = .data[[value_col]],
      fill = highlighted
    )
  ) +
    ggplot2::geom_col(
      width = 0.75
    ) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(
      values = c(
        "FALSE" = other_colour,
        "TRUE" = highlight_colour
      ),
      guide = "none"
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_label,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()
    )
}


## ============================================================================
## 6. Q4 CATEGORICAL DISPLAY
## ============================================================================

#' Create the Q4 display for a categorical candidate effect modifier
#'
#' The upper panel displays either the mean treatment score residual (phi) or
#' the Cox hazard ratio (HR) with 95 percent confidence intervals, by
#' candidate covariate level and region. Region-specific overall values are
#' shown as horizontal dotted reference lines with labelled values. The lower
#' panel displays the proportion of each covariate level within each region,
#' with sample counts annotated above each bar.
#'
#' @param data Analysis data frame.
#' @param phi Numeric score residual vector (required for lower panel and
#'   for metric = "phi").
#' @param x_name Name of the categorical candidate effect modifier.
#' @param region_name Name of the region variable.
#' @param x_levels Optional order for the candidate covariate levels.
#' @param x_labels Optional labels corresponding to x_levels after excluding.
#' @param exclude_levels Character vector of x levels to exclude from the plot.
#' @param metric Either "phi" (mean score residual) or "hr" (Cox hazard ratio).
#' @param hr_table Required when metric = "hr". A data frame from
#'   estimate_hr_by_subgroup() with columns region, subgroup, HR, lower, upper,
#'   events.
#' @param region_overall_hr Required when metric = "hr". A data frame with
#'   columns region, HR, lower, upper representing each region's overall HR.
#' @param min_events Minimum number of events required in a cell to display its
#'   HR estimate. Cells with fewer events are hidden. Default 5.
#' @param region_colours Named vector of region colours.
#' @param title Overall title.
#'
#' @return A patchwork combined plot.
plot_q4_display_categorical <- function(
    data,
    phi,
    x_name,
    region_name = "region",
    x_levels = NULL,
    x_labels = NULL,
    exclude_levels = NULL,
    metric = c("phi", "hr"),
    hr_table = NULL,
    region_overall_hr = NULL,
    y_breaks = NULL,
    min_events = 5L,
    region_colours = c(
      "ROW" = "#D95F5F",
      "US" = "#3C91E6"
    ),
    title = NULL
) {
  metric <- match.arg(metric)
  
  if (!(x_name %in% names(data))) {
    stop(x_name, " was not found in data.", call. = FALSE)
  }
  
  if (!(region_name %in% names(data))) {
    stop(region_name, " was not found in data.", call. = FALSE)
  }
  
  if (length(phi) != nrow(data)) {
    stop("phi must have one value per row of data.", call. = FALSE)
  }
  
  if (metric == "hr") {
    if (is.null(hr_table)) {
      stop("hr_table must be provided when metric = 'hr'.", call. = FALSE)
    }
    
    if (is.null(region_overall_hr)) {
      stop(
        "region_overall_hr must be provided when metric = 'hr'.",
        call. = FALSE
      )
    }
    
    required_hr_columns <- c(
      "region",
      "subgroup",
      "HR",
      "lower",
      "upper",
      "events"
    )
    
    missing_hr_columns <- setdiff(required_hr_columns, names(hr_table))
    
    if (length(missing_hr_columns) > 0L) {
      stop(
        "hr_table is missing required columns: ",
        paste(missing_hr_columns, collapse = ", "),
        call. = FALSE
      )
    }
    
    required_overall_columns <- c("region", "HR", "lower", "upper")
    
    missing_overall_columns <- setdiff(
      required_overall_columns,
      names(region_overall_hr)
    )
    
    if (length(missing_overall_columns) > 0L) {
      stop(
        "region_overall_hr is missing required columns: ",
        paste(missing_overall_columns, collapse = ", "),
        call. = FALSE
      )
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Assemble base data
  ## --------------------------------------------------------------------------
  plot_data <- data.frame(
    x = data[[x_name]],
    region = data[[region_name]],
    phi = as.numeric(phi)
  )
  
  plot_data <- plot_data[stats::complete.cases(plot_data), , drop = FALSE]
  
  plot_data$region <- factor(plot_data$region)
  
  if (!is.null(x_levels)) {
    plot_data$x <- factor(
      plot_data$x,
      levels = x_levels,
      ordered = TRUE
    )
  } else {
    plot_data$x <- factor(plot_data$x)
  }
  
  ## Exclude specified levels (e.g., undetermined).
  if (!is.null(exclude_levels)) {
    plot_data <- plot_data[
      !(as.character(plot_data$x) %in% exclude_levels),
      ,
      drop = FALSE
    ]
    
    plot_data$x <- droplevels(plot_data$x)
  }
  
  if (!is.null(x_labels)) {
    if (length(x_labels) != nlevels(plot_data$x)) {
      stop(
        "x_labels must have one value for each retained x level. ",
        "Currently ",
        nlevels(plot_data$x),
        " levels remain.",
        call. = FALSE
      )
    }
    
    levels(plot_data$x) <- x_labels
  }
  
  ## Region colour completion
  observed_regions <- levels(plot_data$region)
  
  missing_colours <- setdiff(observed_regions, names(region_colours))
  
  if (length(missing_colours) > 0L) {
    additional_colours <- grDevices::hcl.colors(
      length(missing_colours),
      palette = "Dark 3"
    )
    
    names(additional_colours) <- missing_colours
    
    region_colours <- c(region_colours, additional_colours)
  }
  
  n_x_levels <- nlevels(plot_data$x)
  label_x_position <- n_x_levels + 0.35
  
  ## --------------------------------------------------------------------------
  ## Upper panel: metric-specific
  ## --------------------------------------------------------------------------
  if (metric == "phi") {
    upper_summary <- plot_data |>
      dplyr::group_by(x, region) |>
      dplyr::summarise(
        n = dplyr::n(),
        estimate = mean(phi),
        sd_val = stats::sd(phi),
        se_val = sd_val / sqrt(n),
        lower = estimate - stats::qnorm(0.975) * se_val,
        upper = estimate + stats::qnorm(0.975) * se_val,
        .groups = "drop"
      )
    
    region_overall <- plot_data |>
      dplyr::group_by(region) |>
      dplyr::summarise(
        estimate_overall = mean(phi),
        .groups = "drop"
      )
    
    region_overall$label <- sprintf(
      "%s overall mean \u03c6 = %.4f",
      region_overall$region,
      region_overall$estimate_overall
    )
    
    y_axis_label <- expression(
      paste("Mean treatment score residual ", phi)
    )
    
    reference_intercept <- 0
    use_log_scale <- FALSE
    
  } else {
    ## metric == "hr"
    upper_summary <- hr_table
    
    ## Exclude specified levels
    if (!is.null(exclude_levels)) {
      upper_summary <- upper_summary[
        !(as.character(upper_summary$subgroup) %in% exclude_levels),
        ,
        drop = FALSE
      ]
    }
    
    ## Set x factor with the same levels retained in plot_data
    retained_original_levels <- if (!is.null(x_levels)) {
      setdiff(x_levels, exclude_levels)
    } else {
      levels(plot_data$x)
    }
    
    upper_summary$x <- factor(
      as.character(upper_summary$subgroup),
      levels = retained_original_levels,
      ordered = TRUE
    )
    
    ## Relabel to display labels
    if (!is.null(x_labels)) {
      levels(upper_summary$x) <- x_labels
    }
    
    upper_summary$region <- factor(
      as.character(upper_summary$region),
      levels = levels(plot_data$region)
    )
    
    upper_summary$estimate <- upper_summary$HR
    
    ## Apply minimum-events threshold: hide unreliable estimates
    unreliable <- is.na(upper_summary$events) |
      upper_summary$events < min_events
    
    if (any(unreliable)) {
      message(
        "Hiding ",
        sum(unreliable),
        " subgroup(s) with fewer than ",
        min_events,
        " events."
      )
    }
    
    upper_summary$estimate[unreliable] <- NA_real_
    upper_summary$lower[unreliable] <- NA_real_
    upper_summary$upper[unreliable] <- NA_real_
    
    region_overall <- region_overall_hr
    
    region_overall$region <- factor(
      as.character(region_overall$region),
      levels = levels(plot_data$region)
    )
    
    region_overall$estimate_overall <- region_overall$HR
    
    region_overall$label <- sprintf(
      "%s overall HR = %.2f",
      region_overall$region,
      region_overall$estimate_overall
    )
    
    y_axis_label <- "Hazard ratio (ticagrelor vs clopidogrel)"
    
    reference_intercept <- 1
    use_log_scale <- TRUE
  }
  
  ## --------------------------------------------------------------------------
  ## Auto-adjust label vjust to avoid overlap
  ## --------------------------------------------------------------------------
  region_overall <- region_overall[
    order(region_overall$estimate_overall),
    ,
    drop = FALSE
  ]
  
  region_overall$vjust_label <- 1.4
  
  if (nrow(region_overall) >= 2L) {
    region_overall$vjust_label[1] <- 1.4
    region_overall$vjust_label[2] <- -0.6
    
    finite_bounds <- is.finite(upper_summary$lower) &
      is.finite(upper_summary$upper)
    
    if (any(finite_bounds)) {
      y_range <- range(
        c(
          upper_summary$lower[finite_bounds],
          upper_summary$upper[finite_bounds]
        ),
        na.rm = TRUE
      )
      
      y_gap <- diff(y_range)
      
      label_gap <- diff(region_overall$estimate_overall)
      
      close_labels <- length(label_gap) > 0L &&
        !is.na(label_gap[1]) &&
        y_gap > .Machine$double.eps &&
        (abs(label_gap[1]) / y_gap) < 0.08
      
      if (close_labels) {
        region_overall$vjust_label[1] <- 1.7
        region_overall$vjust_label[2] <- -0.9
      }
    }
  }
  
  ## --------------------------------------------------------------------------
  ## Upper panel plot
  ## --------------------------------------------------------------------------
  upper_panel <- ggplot2::ggplot(
    upper_summary,
    ggplot2::aes(
      x = x,
      y = estimate,
      colour = region,
      group = region
    )
  ) +
    ggplot2::geom_hline(
      yintercept = reference_intercept,
      linetype = "dashed",
      colour = "grey60"
    ) +
    ggplot2::geom_hline(
      data = region_overall,
      ggplot2::aes(
        yintercept = estimate_overall,
        colour = region
      ),
      linetype = "dotted",
      linewidth = 0.7,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = lower,
        ymax = upper
      ),
      position = ggplot2::position_dodge(width = 0.45),
      width = 0.12,
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_dodge(width = 0.45),
      size = 2.8,
      na.rm = TRUE
    ) +
    ggplot2::geom_text(
      data = region_overall,
      ggplot2::aes(
        x = label_x_position,
        y = estimate_overall,
        colour = region,
        label = label,
        vjust = vjust_label
      ),
      hjust = 1,
      size = 3.2,
      show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    ggplot2::scale_colour_manual(values = region_colours) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(add = c(0.5, 0.6))
    ) +
    ggplot2::labs(
      x = NULL,
      y = y_axis_label,
      title = title,
      colour = "Region"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (isTRUE(use_log_scale)) {
    if (!is.null(y_breaks)) {
      ## Ensure limits cover both the data range and the user-specified breaks.
      data_range <- range(
        c(upper_summary$lower, upper_summary$upper),
        na.rm = TRUE,
        finite = TRUE
      )
      
      y_limits <- c(
        min(data_range[1], min(y_breaks)),
        max(data_range[2], max(y_breaks))
      )
      
      upper_panel <- upper_panel +
        ggplot2::scale_y_continuous(
          trans = "log10",
          breaks = y_breaks,
          labels = as.character(y_breaks),
          limits = y_limits,
          minor_breaks = NULL,
          oob = scales::squish_infinite
        )
    } else {
      upper_panel <- upper_panel +
        ggplot2::scale_y_continuous(
          trans = "log10"
        )
    }
  } else if (!is.null(y_breaks)) {
    data_range <- range(
      c(upper_summary$lower, upper_summary$upper),
      na.rm = TRUE,
      finite = TRUE
    )
    
    y_limits <- c(
      min(data_range[1], min(y_breaks)),
      max(data_range[2], max(y_breaks))
    )
    
    upper_panel <- upper_panel +
      ggplot2::scale_y_continuous(
        breaks = y_breaks,
        limits = y_limits,
        minor_breaks = NULL
      )
  }
  
  ## --------------------------------------------------------------------------
  ## Lower panel: aspirin-dose proportions within region, with counts labelled
  ## --------------------------------------------------------------------------
  distribution_summary <- plot_data |>
    dplyr::count(region, x, name = "n") |>
    dplyr::group_by(region) |>
    dplyr::mutate(proportion = n / sum(n)) |>
    dplyr::ungroup()
  
  lower_panel <- ggplot2::ggplot(
    distribution_summary,
    ggplot2::aes(
      x = x,
      y = proportion,
      fill = region
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(width = 0.8),
      width = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = format(
          n,
          big.mark = ",",
          scientific = FALSE
        )
      ),
      position = ggplot2::position_dodge(width = 0.8),
      vjust = -0.4,
      size = 3.0,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = region_colours) +
    ggplot2::scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      expand = ggplot2::expansion(mult = c(0, 0.14))
    ) +
    ggplot2::scale_x_discrete(
      expand = ggplot2::expansion(add = c(0.5, 0.6))
    ) +
    ggplot2::labs(
      x = x_name,
      y = "Proportion within region",
      fill = "Region"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 20, hjust = 1)
    )
  
  upper_panel /
    lower_panel +
    patchwork::plot_layout(heights = c(2.1, 1))
}

## ============================================================================
## 7. SUBGROUP-SPECIFIC COX HAZARD RATIOS
## ============================================================================

#' Estimate treatment hazard ratios by region and categorical subgroup
#'
#' Fits a separate Cox model within each combination of region and categorical
#' candidate effect-modifier level.
#'
#' @param data Analysis data frame.
#' @param time_name Survival-time column name.
#' @param event_name Event-indicator column name.
#' @param treatment_name Treatment column name.
#' @param drug_name Investigational treatment label.
#' @param control_name Control treatment label.
#' @param region_name Region column name.
#' @param subgroup_name Categorical subgroup column name.
#' @param subgroup_levels Optional subgroup ordering.
#'
#' @return Data frame of subgroup-specific HRs, confidence intervals, p-values,
#'   sample sizes, and event counts.
estimate_hr_by_subgroup <- function(
    data,
    time_name,
    event_name,
    treatment_name,
    drug_name,
    control_name,
    region_name = "region",
    subgroup_name,
    subgroup_levels = NULL
) {
  required_columns <- c(
    time_name,
    event_name,
    treatment_name,
    region_name,
    subgroup_name
  )
  
  missing_columns <- setdiff(
    required_columns,
    names(data)
  )
  
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  analysis_data <- data.frame(
    time = as.numeric(data[[time_name]]),
    event = as.integer(data[[event_name]]),
    treatment = data[[treatment_name]],
    region = data[[region_name]],
    subgroup = data[[subgroup_name]]
  )
  
  analysis_data <- analysis_data[
    stats::complete.cases(analysis_data),
    ,
    drop = FALSE
  ]
  
  if (any(analysis_data$time <= 0)) {
    stop(
      "All survival times must be positive.",
      call. = FALSE
    )
  }
  
  analysis_data <- analysis_data[
    analysis_data$treatment %in% c(
      control_name,
      drug_name
    ),
    ,
    drop = FALSE
  ]
  
  analysis_data$trt_num <- as.integer(
    analysis_data$treatment == drug_name
  )
  
  analysis_data$region <- factor(
    analysis_data$region
  )
  
  if (!is.null(subgroup_levels)) {
    analysis_data$subgroup <- factor(
      analysis_data$subgroup,
      levels = subgroup_levels,
      ordered = TRUE
    )
  } else {
    analysis_data$subgroup <- factor(
      analysis_data$subgroup
    )
  }
  
  combinations <- unique(
    analysis_data[
      ,
      c(
        "region",
        "subgroup"
      ),
      drop = FALSE
    ]
  )
  
  combinations <- combinations[
    order(
      combinations$region,
      combinations$subgroup
    ),
    ,
    drop = FALSE
  ]
  
  estimates <- lapply(
    seq_len(nrow(combinations)),
    function(index) {
      selected_region <- combinations$region[index]
      selected_subgroup <- combinations$subgroup[index]
      
      subset_data <- analysis_data[
        analysis_data$region == selected_region &
          analysis_data$subgroup == selected_subgroup,
        ,
        drop = FALSE
      ]
      
      n_total <- nrow(subset_data)
      n_control <- sum(
        subset_data$trt_num == 0
      )
      n_drug <- sum(
        subset_data$trt_num == 1
      )
      events_total <- sum(
        subset_data$event
      )
      events_control <- sum(
        subset_data$event[
          subset_data$trt_num == 0
        ]
      )
      events_drug <- sum(
        subset_data$event[
          subset_data$trt_num == 1
        ]
      )
      
      estimable <- (
        n_control > 0 &&
          n_drug > 0 &&
          events_control > 0 &&
          events_drug > 0
      )
      
      if (!estimable) {
        return(
          data.frame(
            region = as.character(selected_region),
            subgroup = as.character(selected_subgroup),
            n = n_total,
            n_control = n_control,
            n_drug = n_drug,
            events = events_total,
            events_control = events_control,
            events_drug = events_drug,
            HR = NA_real_,
            lower = NA_real_,
            upper = NA_real_,
            p_value = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
      
      fit <- try(
        survival::coxph(
          survival::Surv(time, event) ~ trt_num,
          data = subset_data
        ),
        silent = TRUE
      )
      
      if (inherits(fit, "try-error")) {
        return(
          data.frame(
            region = as.character(selected_region),
            subgroup = as.character(selected_subgroup),
            n = n_total,
            n_control = n_control,
            n_drug = n_drug,
            events = events_total,
            events_control = events_control,
            events_drug = events_drug,
            HR = NA_real_,
            lower = NA_real_,
            upper = NA_real_,
            p_value = NA_real_,
            stringsAsFactors = FALSE
          )
        )
      }
      
      fit_summary <- summary(fit)
      
      data.frame(
        region = as.character(selected_region),
        subgroup = as.character(selected_subgroup),
        n = n_total,
        n_control = n_control,
        n_drug = n_drug,
        events = events_total,
        events_control = events_control,
        events_drug = events_drug,
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
        p_value = unname(
          fit_summary$coefficients[
            "trt_num",
            "Pr(>|z|)"
          ]
        ),
        stringsAsFactors = FALSE
      )
    }
  )
  
  result <- do.call(
    rbind,
    estimates
  )
  
  if (!is.null(subgroup_levels)) {
    result$subgroup <- factor(
      result$subgroup,
      levels = subgroup_levels,
      ordered = TRUE
    )
  }
  
  result$region <- factor(
    result$region,
    levels = levels(analysis_data$region)
  )
  
  result
}


#' Plot subgroup-specific treatment hazard ratios
#'
#' @param hr_table Result from estimate_hr_by_subgroup().
#' @param treatment_label Label describing the hazard ratio.
#' @param title Plot title.
#' @param subtitle Plot subtitle.
#' @param region_colours Named vector of region colours.
#' @param x_limits Optional HR-axis limits.
#' @param log_scale Logical; display HR on a logarithmic axis.
#'
#' @return ggplot object.
plot_hr_by_subgroup <- function(
    hr_table,
    treatment_label = "Hazard ratio for investigational treatment versus control",
    title = "Treatment effect by region and aspirin dose",
    subtitle = NULL,
    region_colours = c(
      "ROW" = "#D95F5F",
      "US" = "#3C91E6"
    ),
    x_limits = NULL,
    log_scale = TRUE
) {
  required_columns <- c(
    "region",
    "subgroup",
    "HR",
    "lower",
    "upper"
  )
  
  if (!all(required_columns %in% names(hr_table))) {
    stop(
      "hr_table must contain: ",
      paste(required_columns, collapse = ", "),
      call. = FALSE
    )
  }
  
  plot_data <- hr_table
  
  plot_data$label <- paste(
    plot_data$region,
    plot_data$subgroup,
    sep = "  "
  )
  
  plot_data <- plot_data[
    order(
      plot_data$region,
      plot_data$subgroup
    ),
    ,
    drop = FALSE
  ]
  
  plot_data$label <- factor(
    plot_data$label,
    levels = rev(plot_data$label)
  )
  
  observed_regions <- unique(
    as.character(plot_data$region)
  )
  
  missing_colours <- setdiff(
    observed_regions,
    names(region_colours)
  )
  
  if (length(missing_colours) > 0L) {
    additional_colours <- grDevices::hcl.colors(
      length(missing_colours),
      palette = "Dark 3"
    )
    
    names(additional_colours) <- missing_colours
    
    region_colours <- c(
      region_colours,
      additional_colours
    )
  }
  
  plot_object <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = HR,
      y = label,
      colour = region
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey50"
    ) +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = lower,
        xmax = upper
      ),
      height = 0.15,
      linewidth = 0.7,
      na.rm = TRUE
    ) +
    ggplot2::geom_point(
      size = 2.6,
      na.rm = TRUE
    ) +
    ggplot2::scale_colour_manual(
      values = region_colours
    ) +
    ggplot2::labs(
      x = treatment_label,
      y = NULL,
      title = title,
      subtitle = subtitle,
      colour = "Region"
    ) +
    ggplot2::theme_bw(
      base_size = 12
    ) +
    ggplot2::theme(
      legend.position = "top",
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank()
    )
  
  if (isTRUE(log_scale)) {
    plot_object <- plot_object +
      ggplot2::scale_x_log10(
        limits = x_limits
      )
  } else if (!is.null(x_limits)) {
    plot_object <- plot_object +
      ggplot2::coord_cartesian(
        xlim = x_limits
      )
  }
  
  plot_object
}