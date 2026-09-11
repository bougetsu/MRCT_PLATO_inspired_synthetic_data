# =============================================================================
# compare_summaries.R
#
# Recomputes key summaries from the PLATO-like synthetic dataset and compares
# them with published aggregate results from Wallentin et al. (2009) and
# Mahaffey et al. (2011).
#
# Outputs:
#   synthetic_data/data/comparison_summary.csv
#   figures/plato_fidelity.pdf
#   figures/plato_fidelity.png
#
# Run from the repository root.
# =============================================================================

suppressPackageStartupMessages(library(survival))

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
data_path <- file.path("synthetic_data", "data", "plato_synthetic.csv")
comparison_path <- file.path("synthetic_data", "data", "comparison_summary.csv")
figure_dir <- "figures"

if (!file.exists(data_path)) {
  stop("Synthetic data file not found: ", data_path, call. = FALSE)
}
if (!dir.exists(dirname(comparison_path))) {
  dir.create(dirname(comparison_path), recursive = TRUE)
}
if (!dir.exists(figure_dir)) {
  dir.create(figure_dir, recursive = TRUE)
}

# -----------------------------------------------------------------------------
# Read data and helpers
# -----------------------------------------------------------------------------
d <- read.csv(data_path, stringsAsFactors = FALSE)
d$treatment <- factor(d$treatment, levels = c("clopidogrel", "ticagrelor"))
d$region <- factor(d$region, levels = c("ROW", "US"))

cmp <- list()
add <- function(metric, source, synthetic) {
  cmp[[length(cmp) + 1L]] <<- data.frame(
    metric = metric,
    source = source,
    synthetic = synthetic,
    stringsAsFactors = FALSE
  )
}

hr_ci <- function(fit) {
  s <- summary(fit)
  c(
    hr = unname(s$conf.int[1, "exp(coef)"]),
    lo = unname(s$conf.int[1, "lower .95"]),
    hi = unname(s$conf.int[1, "upper .95"])
  )
}

fmt_hr <- function(x) {
  sprintf("%.2f (%.2f-%.2f)", x["hr"], x["lo"], x["hi"])
}

km_inc <- function(dat, t = 360) {
  sf <- survfit(Surv(time_days, event) ~ 1, data = dat)
  1 - summary(sf, times = t, extend = TRUE)$surv
}

# -----------------------------------------------------------------------------
# Sample sizes and event counts
# -----------------------------------------------------------------------------
add("N total", "18624", as.character(nrow(d)))
add("N ticagrelor", "9333", as.character(sum(d$treatment == "ticagrelor")))
add("N clopidogrel", "9291", as.character(sum(d$treatment == "clopidogrel")))
add("N US", "1413", as.character(sum(d$region == "US")))
add("N ROW", "17211", as.character(sum(d$region == "ROW")))

fig3 <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
region asa_dose treatment    N     E
US     high     ticagrelor   324   40
US     high     clopidogrel  352   27
US     mid      ticagrelor   22    2
US     mid      clopidogrel  16    2
US     low      ticagrelor   284   19
US     low      clopidogrel  263   24
US     undetermined ticagrelor 77  23
US     undetermined clopidogrel 75 14
ROW    high     ticagrelor   140   28
ROW    high     clopidogrel  140   23
ROW    mid      ticagrelor   503   62
ROW    mid      clopidogrel  511   63
ROW    low      ticagrelor   7449  546
ROW    low      clopidogrel  7443  699
ROW    undetermined ticagrelor 534 144
ROW    undetermined clopidogrel 491 162
")

for (i in seq_len(nrow(fig3))) {
  sub <- d[
    d$region == fig3$region[i] &
      d$asa_dose == fig3$asa_dose[i] &
      d$treatment == fig3$treatment[i],
    ,
    drop = FALSE
  ]
  tag <- paste(fig3$region[i], fig3$asa_dose[i], substr(fig3$treatment[i], 1, 4))
  add(paste0("N [", tag, "]"), as.character(fig3$N[i]), as.character(nrow(sub)))
  add(paste0("E [", tag, "]"), as.character(fig3$E[i]), as.character(sum(sub$event)))
}

# -----------------------------------------------------------------------------
# 12-month cumulative incidence
# -----------------------------------------------------------------------------
add(
  "Cum. incidence @360d, ticagrelor",
  "9.8%",
  sprintf("%.1f%%", 100 * km_inc(d[d$treatment == "ticagrelor", ]))
)
add(
  "Cum. incidence @360d, clopidogrel",
  "11.7%",
  sprintf("%.1f%%", 100 * km_inc(d[d$treatment == "clopidogrel", ]))
)

# -----------------------------------------------------------------------------
# Hazard ratios
# -----------------------------------------------------------------------------
add(
  "HR overall (Wallentin Fig 1)",
  "0.84 (0.77-0.92)",
  fmt_hr(hr_ci(coxph(Surv(time_days, event) ~ treatment, data = d)))
)
add(
  "HR US (Mahaffey Table 2)",
  "1.27 (0.92-1.75)",
  fmt_hr(hr_ci(coxph(Surv(time_days, event) ~ treatment, data = d[d$region == "US", ])))
)
add(
  "HR ROW (Mahaffey Table 2)",
  "0.81 (0.74-0.90)",
  fmt_hr(hr_ci(coxph(Surv(time_days, event) ~ treatment, data = d[d$region == "ROW", ])))
)

fig3_hr <- list(
  c("US", "high", "1.62 (0.99-2.64)"),
  c("US", "low",  "0.73 (0.40-1.33)"),
  c("ROW", "high", "1.23 (0.71-2.14)"),
  c("ROW", "mid",  "1.00 (0.71-1.42)"),
  c("ROW", "low",  "0.78 (0.69-0.87)")
)
for (x in fig3_hr) {
  sub <- d[d$region == x[1] & d$asa_dose == x[2], ]
  add(
    sprintf("HR %s ASA=%s (Mahaffey Fig 3)", x[1], x[2]),
    x[3],
    fmt_hr(hr_ci(coxph(Surv(time_days, event) ~ treatment, data = sub)))
  )
}

add(
  "HR ASA high pooled",
  "1.45 (1.01-2.09)",
  fmt_hr(hr_ci(coxph(Surv(time_days, event) ~ treatment, data = d[d$asa_dose == "high", ])))
)
add(
  "HR ASA low pooled",
  "0.77 (0.69-0.86)",
  fmt_hr(hr_ci(coxph(Surv(time_days, event) ~ treatment, data = d[d$asa_dose == "low", ])))
)

# Published interaction summaries, included in the audit table but not used as
# direct construction targets for the fidelity figure.
di <- d[d$asa_dose %in% c("low", "high"), ]
di$asa_dose <- factor(di$asa_dose, levels = c("low", "high"))
m_int <- coxph(Surv(time_days, event) ~ treatment * asa_dose, data = di)
p_int <- summary(m_int)$coefficients[
  "treatmentticagrelor:asa_dosehigh", "Pr(>|z|)"
]
add("Trt x ASA interaction p (Mahaffey)", "0.00006", sprintf("%.5f", p_int))

m_reg <- coxph(Surv(time_days, event) ~ treatment * region, data = d)
p_reg <- summary(m_reg)$coefficients[
  "treatmentticagrelor:regionUS", "Pr(>|z|)"
]
add("Trt x region interaction p (Mahaffey ~0.01)", "~0.01", sprintf("%.4f", p_reg))

# -----------------------------------------------------------------------------
# Number at risk
# -----------------------------------------------------------------------------
risk_days <- c(60, 120, 180, 240, 300, 360)
src_atrisk <- c(
  8628 + 8521,
  8460 + 8362,
  8219 + 8124,
  6743 + 6650,
  5161 + 5096,
  4147 + 4047
)
syn_atrisk <- sapply(risk_days, function(t) sum(d$time_days >= t))
for (i in seq_along(risk_days)) {
  add(
    sprintf("At risk @%dd", risk_days[i]),
    as.character(src_atrisk[i]),
    as.character(syn_atrisk[i])
  )
}

# -----------------------------------------------------------------------------
# Baseline covariate distributions by region
# -----------------------------------------------------------------------------
pct <- function(x) sprintf("%.1f%%", 100 * mean(x))
med_iqr <- function(x) {
  sprintf(
    "%.0f (%.0f-%.0f)",
    median(x), quantile(x, .25), quantile(x, .75)
  )
}

t1 <- list(
  list("Age US",         "61 (53-70)",  function(z) med_iqr(z$age)),
  list("Age ROW",        "62 (54-71)",  function(z) med_iqr(z$age)),
  list("Weight US",      "87 (75-100)", function(z) med_iqr(z$weight_kg)),
  list("Weight ROW",     "80 (70-89)",  function(z) med_iqr(z$weight_kg)),
  list("BMI US",         "29.1",        function(z) sprintf("%.1f", median(z$bmi))),
  list("BMI ROW",        "27.3",        function(z) sprintf("%.1f", median(z$bmi))),
  list("Female US",      "28.7%",       function(z) pct(z$sex == "female")),
  list("Female ROW",     "28.4%",       function(z) pct(z$sex == "female")),
  list("Diabetes US",    "33.4%",       function(z) pct(z$diabetes)),
  list("Diabetes ROW",   "24.4%",       function(z) pct(z$diabetes)),
  list("Prior PCI US",   "29.4%",       function(z) pct(z$prior_pci)),
  list("Prior PCI ROW",  "12.1%",       function(z) pct(z$prior_pci)),
  list("Prior CABG US",  "16.7%",       function(z) pct(z$prior_cabg)),
  list("Prior CABG ROW", "5.1%",        function(z) pct(z$prior_cabg)),
  list("STEMI US",       "15.7%",       function(z) pct(z$stemi)),
  list("STEMI ROW",      "39.6%",       function(z) pct(z$stemi)),
  list("Beta-blocker US","80.9%",       function(z) pct(z$beta_blocker)),
  list("Beta-blocker ROW","70.8%",      function(z) pct(z$beta_blocker)),
  list("White race US",  "89.3%",       function(z) pct(z$race == "White")),
  list("White race ROW", "91.9%",       function(z) pct(z$race == "White"))
)

for (row in t1) {
  rg <- if (grepl("US$", row[[1]])) "US" else "ROW"
  add(paste0("Table1: ", row[[1]]), row[[2]], row[[3]](d[d$region == rg, ]))
}

# -----------------------------------------------------------------------------
# Figure 2A-style interaction screen (diagnostic only)
# -----------------------------------------------------------------------------
cat("\n--- Treatment-by-covariate interaction screen ---\n")
covs <- c(
  "asa_dose", "region", "age", "sex", "race", "weight_kg", "bmi",
  "smoking", "diabetes", "prior_mi", "prior_pci", "prior_cabg",
  "beta_blocker", "stemi", "troponin_pos"
)
int_tab <- data.frame(covariate = covs, interaction_p = NA_real_)
for (i in seq_along(covs)) {
  cv <- covs[i]
  dd <- d
  if (cv == "asa_dose") {
    dd <- d[d$asa_dose %in% c("low", "high"), ]
    dd$asa_dose <- droplevels(factor(dd$asa_dose))
  }
  m0 <- coxph(
    as.formula(sprintf("Surv(time_days, event) ~ treatment + %s", cv)),
    data = dd
  )
  m1 <- coxph(
    as.formula(sprintf("Surv(time_days, event) ~ treatment * %s", cv)),
    data = dd
  )
  int_tab$interaction_p[i] <- anova(m0, m1, test = "LRT")$`Pr(>|Chi|)`[2]
}
int_tab <- int_tab[order(int_tab$interaction_p), ]
print(int_tab, row.names = FALSE)

# -----------------------------------------------------------------------------
# Write audit table
# -----------------------------------------------------------------------------
res <- do.call(rbind, cmp)
write.csv(res, comparison_path, row.names = FALSE)
cat("\nWrote ", comparison_path, "\n", sep = "")

# -----------------------------------------------------------------------------
# Fidelity figure: Panel A (HRs) + Panel B (cumulative incidence)
# -----------------------------------------------------------------------------
plot_packages <- c("ggplot2", "ggsci", "patchwork")
missing_plot_packages <- plot_packages[
  !vapply(plot_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_plot_packages) > 0L) {
  warning(
    "Skipping fidelity figure because these packages are not installed: ",
    paste(missing_plot_packages, collapse = ", "),
    call. = FALSE
  )
} else {
  library(ggplot2)
  library(patchwork)

  npg <- ggsci::pal_npg("nrc")(10)
  comparison_colours <- c("Published" = npg[1], "Synthetic" = npg[2])
  treatment_colours <- c("clopidogrel" = npg[2], "ticagrelor" = npg[1])

  parse_hr_ci <- function(x) {
    vals <- regmatches(x, gregexpr("[0-9]+(?:\\.[0-9]+)?", x, perl = TRUE))[[1]]
    vals <- as.numeric(vals)
    if (length(vals) < 3L) {
      return(c(hr = NA_real_, lo = NA_real_, hi = NA_real_))
    }
    c(hr = vals[1], lo = vals[2], hi = vals[3])
  }

  hr_metrics <- c(
    "HR overall (Wallentin Fig 1)",
    "HR US (Mahaffey Table 2)",
    "HR ROW (Mahaffey Table 2)",
    "HR US ASA=low (Mahaffey Fig 3)",
    "HR US ASA=high (Mahaffey Fig 3)",
    "HR ROW ASA=low (Mahaffey Fig 3)",
    "HR ROW ASA=mid (Mahaffey Fig 3)",
    "HR ROW ASA=high (Mahaffey Fig 3)"
  )
  hr_labels <- c(
    "Overall",
    "US overall",
    "RoW overall",
    "US: low ASA",
    "US: high ASA",
    "RoW: low ASA",
    "RoW: intermediate ASA",
    "RoW: high ASA"
  )

  hr_res <- res[match(hr_metrics, res$metric), ]
  hr_long <- do.call(
    rbind,
    lapply(seq_len(nrow(hr_res)), function(i) {
      src <- parse_hr_ci(hr_res$source[i])
      syn <- parse_hr_ci(hr_res$synthetic[i])
      rbind(
        data.frame(
          subgroup = hr_labels[i], estimate = src["hr"], lower = src["lo"],
          upper = src["hi"], series = "Published", stringsAsFactors = FALSE
        ),
        data.frame(
          subgroup = hr_labels[i], estimate = syn["hr"], lower = syn["lo"],
          upper = syn["hi"], series = "Synthetic", stringsAsFactors = FALSE
        )
      )
    })
  )
  hr_long$subgroup <- factor(hr_long$subgroup, levels = rev(hr_labels))
  hr_long$series <- factor(hr_long$series, levels = c("Published", "Synthetic"))

  dodge <- position_dodge(width = 0.55)
  p_hr <- ggplot(
    hr_long,
    aes(x = subgroup, y = estimate, colour = series, shape = series)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_errorbar(
      aes(ymin = lower, ymax = upper),
      position = dodge,
      width = 0.16,
      linewidth = 0.55
    ) +
    geom_point(position = dodge, size = 2.5) +
    coord_flip() +
    scale_y_log10() +
    scale_colour_manual(values = comparison_colours) +
    scale_shape_manual(values = c("Published" = 16, "Synthetic" = 17)) +
    labs(x = NULL, y = "Hazard ratio", colour = NULL, shape = NULL) +
    theme_bw(base_size = 10.5) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank()
    )

  sf <- survfit(Surv(time_days, event) ~ treatment, data = d)
  km_df <- data.frame(
    time = sf$time,
    incidence = 100 * (1 - sf$surv),
    treatment = sub(
      "^treatment=", "",
      rep(names(sf$strata), as.integer(sf$strata))
    ),
    stringsAsFactors = FALSE
  )
  km_df$treatment <- factor(
    km_df$treatment,
    levels = c("clopidogrel", "ticagrelor")
  )

  published_12m <- data.frame(
    treatment = factor(
      c("clopidogrel", "ticagrelor"),
      levels = c("clopidogrel", "ticagrelor")
    ),
    incidence = c(11.7, 9.8)
  )

  p_km <- ggplot(km_df, aes(x = time, y = incidence, colour = treatment)) +
    geom_step(linewidth = 0.75) +
    geom_hline(
      data = published_12m,
      aes(yintercept = incidence, colour = treatment),
      linetype = "dotted",
      linewidth = 0.65,
      show.legend = FALSE
    ) +
    scale_colour_manual(
      values = treatment_colours,
      labels = c("Clopidogrel", "Ticagrelor")
    ) +
    scale_x_continuous(breaks = seq(0, 360, by = 60), limits = c(0, 360)) +
    labs(
      x = "Days from randomization",
      y = "Cumulative incidence (%)",
      colour = NULL
    ) +
    theme_bw(base_size = 10.5) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )

  fidelity_figure <- p_hr + p_km +
    plot_layout(widths = c(1.15, 1)) +
    plot_annotation(tag_levels = "A") &
    theme(
      plot.tag = element_text(face = "bold", size = 12),
      plot.tag.position = c(0, 1)
    )

  pdf_path <- file.path(figure_dir, "plato_fidelity.pdf")
  png_path <- file.path(figure_dir, "plato_fidelity.png")

  ggsave(pdf_path, fidelity_figure, width = 11.5, height = 5.3, units = "in")
  ggsave(png_path, fidelity_figure, width = 11.5, height = 5.3, units = "in", dpi = 300)

  cat("Wrote ", pdf_path, "\n", sep = "")
  cat("Wrote ", png_path, "\n", sep = "")
}
