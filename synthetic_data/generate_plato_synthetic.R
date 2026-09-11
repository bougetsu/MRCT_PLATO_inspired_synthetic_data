# =============================================================================
# generate_plato_synthetic.R
#
# Generates a synthetic individual-patient dataset for the PLATO trial primary
# composite endpoint (cardiovascular death / MI / stroke), reproducing the
# published summary statistics of:
#   * Wallentin et al. NEJM 2009;361:1045-57  (primary results, Table 1, Fig 1)
#   * Mahaffey et al. Circulation 2011;124:544-54 (regional analysis: Table 1,
#     Table 2, Figure 2A, Figure 3)
#
# The data are SYNTHETIC. No real patient records are used. See DATA_GENERATION.md
# for the full provenance of every number and modelling choice.
#
# Design summary
# --------------
# 1. Patients are laid out on a FIXED cross-tabulation of
#        region (US / ROW) x maintenance-aspirin-dose x treatment arm,
#    with cell sample sizes (N) and event counts (E) taken verbatim from
#    Mahaffey Figure 3 + Table 2 (reconciled with an "undetermined ASA dose"
#    cell so all margins equal the published totals). Because both N and E are
#    fixed per cell, sample sizes and event counts are reproduced exactly; Cox
#    hazard ratios are reproduced closely but are not constrained to match exactly.
# 2. Baseline covariates are drawn WITHIN region from a Gaussian copula that
#    preserves plausible clinical correlations, calibrated to the exact
#    US-vs-ROW marginal distributions of Mahaffey Table 1.
# 3. Which patients within a cell experience the (fixed number of) events is
#    chosen by a prognostic score, giving covariates a realistic prognostic
#    role without altering fixed cell event counts or deliberately introducing
#    additional treatment interactions.
# 4. Event times follow the Wallentin Figure 1 cumulative-incidence shape;
#    administrative censoring reproduces the Figure 1 number-at-risk table.
# =============================================================================

set.seed(20090910)  # NEJM publication date of PLATO primary results

# Run from the repository root.
data_dir <- file.path("synthetic_data", "data")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
out_csv <- file.path(data_dir, "plato_synthetic.csv")

# -----------------------------------------------------------------------------
# 1. Fixed strata: region x ASA maintenance dose x arm, with exact N and E
#    Source: Mahaffey Figure 3 (N, E per cell) + Table 2 (region totals) +
#    Wallentin (arm totals). The "undetermined" ASA cell reconciles Figure 3
#    (which only covers patients with a defined median maintenance dose) with
#    the full arm/region/event totals. Derivation is documented in
#    DATA_GENERATION.md.
# -----------------------------------------------------------------------------
strata <- read.table(header = TRUE, stringsAsFactors = FALSE, text = "
region asa_dose        treatment    N     E
US     high            ticagrelor   324   40
US     high            clopidogrel  352   27
US     mid             ticagrelor   22    2
US     mid             clopidogrel  16    2
US     low             ticagrelor   284   19
US     low             clopidogrel  263   24
US     undetermined    ticagrelor   77    23
US     undetermined    clopidogrel  75    14
ROW    high            ticagrelor   140   28
ROW    high            clopidogrel  140   23
ROW    mid             ticagrelor   503   62
ROW    mid             clopidogrel  511   63
ROW    low             ticagrelor   7449  546
ROW    low             clopidogrel  7443  699
ROW    undetermined    ticagrelor   534   144
ROW    undetermined    clopidogrel  491   162
")

stopifnot(sum(strata$N) == 18624)
stopifnot(sum(strata$N[strata$treatment == "ticagrelor"]) == 9333)
stopifnot(sum(strata$N[strata$treatment == "clopidogrel"]) == 9291)
stopifnot(sum(strata$E[strata$treatment == "ticagrelor"]) == 864)
stopifnot(sum(strata$E[strata$treatment == "clopidogrel"]) == 1014)

# Expand cells into one row per patient
pat <- strata[rep(seq_len(nrow(strata)), strata$N), c("region", "asa_dose", "treatment")]
pat$cell_E <- rep(strata$E, strata$N)
pat$cell_id <- rep(seq_len(nrow(strata)), strata$N)
rownames(pat) <- NULL
n_total <- nrow(pat)

# -----------------------------------------------------------------------------
# 2. Baseline covariates via a within-region Gaussian copula.
#    Covariate set = intersection of (Mahaffey Table 1) and (Mahaffey Figure 2A
#    effect-modifier list): age, sex, race, weight, BMI, prior MI, prior PCI,
#    prior CABG, smoking status, diabetes, troponin, beta-blocker, STEMI.
#    Marginals (US vs ROW) are taken from Mahaffey Table 1.
# -----------------------------------------------------------------------------
library(mvtnorm)

# Region-specific marginal targets (from Mahaffey Table 1).
# Continuous: median and IQR (Q1, Q3); SD approximated as IQR / 1.349.
marg <- list(
  US = list(
    age    = c(m = 61, q1 = 53, q3 = 70),
    weight = c(m = 87, q1 = 75, q3 = 100),
    bmi    = c(m = 29.1, q1 = 25.7, q3 = 33.1),
    p_female   = 406 / 1413,
    p_diabetes = 472 / 1413,
    p_prior_mi = 387 / 1413,
    p_prior_pci = 415 / 1413,
    p_prior_cabg = 236 / 1413,
    p_beta_blocker = 1142 / 1413,
    p_stemi    = 222 / 1413,
    p_troponin_pos = 1176 / 1413,
    race_counts = c(White = 1262, Black = 137, Asian = 9, Other = 5),
    smoke_counts = c(Never = 416, Ex = 481, Current = 515)
  ),
  ROW = list(
    age    = c(m = 62, q1 = 54, q3 = 71),
    weight = c(m = 80, q1 = 70, q3 = 89),
    bmi    = c(m = 27.3, q1 = 24.7, q3 = 30.2),
    p_female   = 4882 / 17211,
    p_diabetes = 4190 / 17211,
    p_prior_mi = 3437 / 17211,
    p_prior_pci = 2077 / 17211,
    p_prior_cabg = 870 / 17211,
    p_beta_blocker = 12169 / 17211,
    p_stemi    = 6804 / 17211,
    p_troponin_pos = 13913 / 17211,
    race_counts = c(White = 15815, Black = 92, Asian = 1087, Other = 216),
    smoke_counts = c(Never = 6840, Ex = 4195, Current = 6163)
  )
)

# Latent variable order for the copula
lat_vars <- c("age", "weight", "bmi", "female", "diabetes", "prior_mi",
              "prior_pci", "prior_cabg", "beta_blocker", "stemi",
              "troponin_pos", "race", "smoke")
p <- length(lat_vars)

# Plausible clinical correlation matrix (documented, not published).
R <- diag(p)
dimnames(R) <- list(lat_vars, lat_vars)
set_cor <- function(a, b, r) { R[a, b] <<- r; R[b, a] <<- r }
set_cor("weight", "bmi", 0.70)
set_cor("age", "weight", -0.10)
set_cor("age", "diabetes", 0.15)
set_cor("bmi", "diabetes", 0.25)
set_cor("prior_mi", "prior_pci", 0.40)
set_cor("prior_mi", "prior_cabg", 0.30)
set_cor("prior_pci", "prior_cabg", 0.35)
set_cor("female", "weight", -0.30)
set_cor("stemi", "troponin_pos", 0.30)
set_cor("stemi", "prior_mi", -0.10)
set_cor("age", "stemi", -0.05)

# Ensure positive-definiteness (shrink toward identity if needed)
ev <- min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
if (ev <= 1e-6) {
  a <- 0.98
  R <- a * R + (1 - a) * diag(p)
}

# Map a standard-normal latent vector to a binary variable with exactly k ones
# (patients with the largest latent values become 1), preserving copula rank
# correlation.
assign_binary_exact <- function(z, k) {
  out <- integer(length(z))
  if (k > 0) out[order(z, decreasing = TRUE)[seq_len(k)]] <- 1L
  out
}

# Map a latent vector to categories with exact counts (ordered by latent value)
assign_categorical_exact <- function(z, counts) {
  ord <- order(z)
  lab <- rep(names(counts), counts)
  out <- character(length(z))
  out[ord] <- lab
  out
}

gen_region <- function(rg, idx) {
  n <- length(idx)
  m <- marg[[rg]]
  Z <- rmvnorm(n, sigma = R)
  colnames(Z) <- lat_vars
  df <- data.frame(row = idx)

  # Continuous margins: median -> mean, IQR -> SD, mapped through the Gaussian
  sd_from_iqr <- function(v) unname((v["q3"] - v["q1"]) / 1.349)
  df$age    <- round(pmin(95, pmax(30, m$age["m"]    + sd_from_iqr(m$age)    * Z[, "age"])))
  df$weight <- round(pmin(180, pmax(40, m$weight["m"] + sd_from_iqr(m$weight) * Z[, "weight"])), 1)
  df$bmi    <- round(pmin(60, pmax(15, m$bmi["m"]     + sd_from_iqr(m$bmi)    * Z[, "bmi"])), 1)

  # Binary margins with exact counts
  df$sex          <- ifelse(assign_binary_exact(Z[, "female"],       round(m$p_female       * n)) == 1, "female", "male")
  df$diabetes     <- assign_binary_exact(Z[, "diabetes"],     round(m$p_diabetes     * n))
  df$prior_mi     <- assign_binary_exact(Z[, "prior_mi"],     round(m$p_prior_mi     * n))
  df$prior_pci    <- assign_binary_exact(Z[, "prior_pci"],    round(m$p_prior_pci    * n))
  df$prior_cabg   <- assign_binary_exact(Z[, "prior_cabg"],   round(m$p_prior_cabg   * n))
  df$beta_blocker <- assign_binary_exact(Z[, "beta_blocker"], round(m$p_beta_blocker * n))
  df$stemi        <- assign_binary_exact(Z[, "stemi"],        round(m$p_stemi        * n))
  df$troponin_pos <- assign_binary_exact(Z[, "troponin_pos"], round(m$p_troponin_pos * n))

  # Categorical margins with exact counts (rescaled to region n)
  rc <- round(m$race_counts / sum(m$race_counts) * n)
  rc[1] <- rc[1] + (n - sum(rc))                      # absorb rounding in the largest cell
  df$race <- assign_categorical_exact(Z[, "race"], rc)
  sc <- round(m$smoke_counts / sum(m$smoke_counts) * n)
  sc[3] <- sc[3] + (n - sum(sc))
  df$smoking <- assign_categorical_exact(Z[, "smoke"], sc)

  df
}

pat$order_id <- seq_len(n_total)
cov_list <- lapply(c("US", "ROW"), function(rg) {
  idx <- which(pat$region == rg)
  # shuffle so covariates are independent of ASA-dose/arm cell assignment
  gen_region(rg, sample(idx))
})
cov <- do.call(rbind, cov_list)
cov <- cov[order(cov$row), ]
pat <- cbind(pat, cov[, setdiff(names(cov), "row")])

# -----------------------------------------------------------------------------
# 3. Event assignment: within each cell, the fixed number of events (cell_E) is
#    allocated to patients with probability proportional to a prognostic score.
#    Applied identically in both arms, so it adds prognostic signal without
#    creating any treatment interaction or changing counts/hazard ratios.
# -----------------------------------------------------------------------------
lp_prog <- with(pat,
  0.030 * (age - 62) +
  0.35  * diabetes +
  0.30  * prior_mi +
  0.20  * prior_cabg +
  0.25  * (troponin_pos) +
  0.20  * stemi -
  0.010 * (weight - 80))
pat$event <- 0L
for (cid in unique(pat$cell_id)) {
  rows <- which(pat$cell_id == cid)
  e <- pat$cell_E[rows[1]]
  if (e > 0) {
    w <- exp(lp_prog[rows])
    winners <- sample(rows, size = e, prob = w)
    pat$event[winners] <- 1L
  }
}
stopifnot(sum(pat$event) == 1878)  # 864 + 1014

# -----------------------------------------------------------------------------
# 4. Times.
#    Administrative censoring reproduces the Wallentin Figure 1 number-at-risk
#    table; event times follow the Figure 1 cumulative-incidence shape.
# -----------------------------------------------------------------------------
# Number-at-risk fractions (average of the two arms), months 0..12 -> days
risk_days <- c(0, 60, 120, 180, 240, 300, 360)
risk_frac <- c(1, 0.9208, 0.9033, 0.8775, 0.7191, 0.5508, 0.4400)

# Target cumulative fraction of the 12-month events by day, read from the
# Figure 1 curve shape (steep rise in the first weeks, then flattening).
inc_days <- c(0, 7, 15, 30, 60, 90, 120, 180, 240, 300, 360)
inc_frac <- c(0, 0.12, 0.26, 0.46, 0.62, 0.71, 0.78, 0.88, 0.94, 0.975, 1.00)

# Recover a pure-censoring survival curve (add back the events removed from the
# risk set) and sample potential follow-up A_i by inverse transform.
overall_event_rate <- sum(pat$event) / n_total
# Administrative-censoring survival S_C(t) = P(admin follow-up >= t), recovered
# from the Figure 1 risk set. The risk set satisfies
#   risk_frac(t) = rate*(1 - F_event(t)) + (1 - rate)*S_C(t),
# because an at-risk patient is either a not-yet-event case or a not-yet-censored
# non-event; inverting gives S_C exactly (rate = overall event rate).
overall_event_rate <- sum(pat$event) / n_total
F_event_at_risk <- approx(inc_days, inc_frac, risk_days, rule = 2)$y
cens_surv <- (risk_frac - overall_event_rate * (1 - F_event_at_risk)) /
             (1 - overall_event_rate)
cens_surv <- pmin(1, pmax(0, cens_surv))
cens_surv[1] <- 1
sample_admin <- function(n) {
  u <- runif(n)
  # invert S_C: patients with u below the 12-month floor reach the 360-day cap
  out <- approx(x = rev(cens_surv), y = rev(risk_days), xout = u, rule = 2)$y
  pmin(360, pmax(1, out))
}
pat$admin_fu <- sample_admin(n_total)

# Event times: drawn independently from the Figure 1 cumulative-incidence shape.
# A patient counted as a primary event necessarily evented before their own
# administrative censoring, so event times are NOT truncated to admin_fu
# (truncation would make censoring informative and bias the KM curve upward).
sample_event_time <- function(m) {
  u <- runif(m)
  pmin(360, pmax(0.5, approx(inc_frac, inc_days, u, rule = 2)$y))
}
ev <- pat$event == 1L
pat$time_days <- pat$admin_fu               # non-events: administrative censoring
pat$time_days[ev] <- sample_event_time(sum(ev))
pat$time_days <- round(pat$time_days, 1)

# -----------------------------------------------------------------------------
# 5. Assemble and write the CSV
# -----------------------------------------------------------------------------
final <- data.frame(
  patient_id   = sprintf("P%05d", seq_len(n_total)),
  treatment    = pat$treatment,
  region       = pat$region,
  asa_dose     = factor(pat$asa_dose, levels = c("low", "mid", "high", "undetermined")),
  age          = pat$age,
  sex          = pat$sex,
  race         = pat$race,
  weight_kg    = pat$weight,
  bmi          = pat$bmi,
  smoking      = pat$smoking,
  diabetes     = pat$diabetes,
  prior_mi     = pat$prior_mi,
  prior_pci    = pat$prior_pci,
  prior_cabg   = pat$prior_cabg,
  beta_blocker = pat$beta_blocker,
  stemi        = pat$stemi,
  troponin_pos = pat$troponin_pos,
  time_days    = pat$time_days,
  event        = pat$event,
  stringsAsFactors = FALSE
)
# shuffle row order so structure is not encoded in row position
final <- final[sample(nrow(final)), ]
final$patient_id <- sprintf("P%05d", seq_len(n_total))

write.csv(final, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "with", nrow(final), "rows and", ncol(final), "columns\n")
