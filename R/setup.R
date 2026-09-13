## Shared setup: dependencies, data loading, labels and model helpers.

required_packages <- c(
  "survival", "rms", "gtsummary", "survminer", "ggplot2", "dplyr", "scales"
)

load_dependencies <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace,
                                       logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing packages: ", paste(missing, collapse = ", "),
         "\nInstall with install.packages(c(",
         paste(sprintf('\"%s\"', missing), collapse = ", "), "))")
  }
  invisible(lapply(required_packages, function(p)
    suppressPackageStartupMessages(library(p, character.only = TRUE))))
}

load_septralu <- function(path = NULL) {
  if (is.null(path)) {
    candidates <- c("septralu_retreatment.rds",
                    file.path("..", "septralu_retreatment.rds"))
    path <- candidates[file.exists(candidates)][1]
    if (is.na(path)) stop("septralu_retreatment.rds not found; pass 'path'.")
  }
  readRDS(path)
}

variable_labels <- list(
  age                               = "Age at re-treatment, years",
  sex                               = "Sex",
  primary_site                      = "Primary tumour site",
  ki67                              = "Ki-67 index, %",
  grade                             = "WHO grade",
  ecog                              = "ECOG PS",
  functioning_tumor                 = "Functioning tumour",
  peritoneal_mets                   = "Peritoneal metastases",
  liver_mets                        = "Liver metastases",
  n_metastatic_sites                = "Number of metastatic sites",
  pet_ga_heterogeneity              = "Intratumoural 68Ga-PET heterogeneity",
  retreatment_interval_months       = "Interval from last I-PRRT to R-PRRT, months",
  metastasis_to_retreatment_months  = "Interval from metastasis to R-PRRT, months",
  n_cycles                          = "Number of R-PRRT cycles",
  dose_reduced_c1                   = "Cycle 1 dose",
  discontinuation_reason            = "Reason for end of treatment",
  recist                            = "RECIST 1.1 response",
  clinical_response                 = "Clinical response",
  biochemical_response_cga          = "Biochemical response (chromogranin A)",
  biochemical_response_secreted     = "Biochemical response (secreted substance)"
)

## Patients with a recorded RECIST category. Response rates are computed among
## these; in the remaining patients the registry holds no response category.
recist_assessed <- function(data) {
  data[!is.na(data$recist) & data$recist != "Not recorded", ]
}

## Parsimonious multivariable model (appropriate for the number of events):
## proliferation (Ki-67), performance status, disease burden (metastatic sites)
## and primary tumour site (pancreas vs others).
## ECOG is modelled as an ordinal score (per 1-point increase). A 0-1 vs 2+
## dichotomy is uninformative in this cohort, where only 6 patients have a
## performance status above 1 and none above 2.
cox_covariates <- c("ki67", "ecog", "n_metastatic_sites", "pancreas")

## Analysis dataset for one endpoint: patients with follow-up, ECOG and number
## of metastatic sites, with the Nelson-Aalen cumulative hazard used for imputation.
model_data <- function(data, endpoint = c("os", "pfs")) {
  endpoint <- match.arg(endpoint)
  d <- data.frame(
    time               = data[[paste0(endpoint, "_time")]],
    event              = data[[paste0(endpoint, "_event")]],
    ki67               = data$ki67,
    ecog               = as.numeric(as.character(data$ecog)),
    n_metastatic_sites = as.numeric(data$n_metastatic_sites),
    pancreas           = as.numeric(data$primary_site_pancreas == "Pancreas"),
    grade              = data$grade
  )
  d <- d[stats::complete.cases(d[, c("time", "event", "ecog", "n_metastatic_sites")]), ]
  d$ecog_1plus <- as.numeric(d$ecog >= 1)
  h <- survival::basehaz(survival::coxph(survival::Surv(time, event) ~ 1, data = d))
  d$cumhaz <- h$hazard[match(d$time, h$time)]
  d
}

## Multiple imputation of Ki-67 (50 imputations, predictive mean matching).
## The imputation model includes the other model covariates, WHO grade and the
## outcome (cumulative hazard and event indicator).
impute_ki67 <- function(d, n_impute = 50, seed = 2026) {
  set.seed(seed)
  Hmisc::aregImpute(~ ki67 + I(ecog) + I(n_metastatic_sites) + pancreas + grade +
                      cumhaz + event, data = d, n.impute = n_impute, pr = FALSE)
}

cox_formula <- function(covariates) {
  stats::as.formula(paste("survival::Surv(time, event) ~",
                          paste(covariates, collapse = " + ")))
}

## Hazard ratios over the cohort interquartile range for Ki-67 (15 vs 2%) and
## the number of metastatic sites (4 vs 2), per 1 point for ECOG, and for
## pancreas vs others.
hazard_ratios <- function(fit, covariates) {
  limits <- list(ki67 = c(2, 15), ecog = c(0, 1), ecog_1plus = c(0, 1),
                 n_metastatic_sites = c(2, 4), pancreas = c(0, 1))
  s <- do.call(summary, c(list(fit), limits[covariates]))
  s <- s[s[, "Type"] == 2, , drop = FALSE]
  data.frame(term = covariates, HR = s[, "Effect"], low = s[, "Lower 0.95"],
             high = s[, "Upper 0.95"], p = stats::anova(fit)[covariates, "P"],
             n = fit$stats[["Obs"]], events = fit$stats[["Events"]], row.names = NULL)
}

## Cox model after multiple imputation, combined with Rubin's rules. Harrell's C
## is averaged across imputations and the global proportional-hazards test is
## summarised by its median p-value.
imputed_cox <- function(d, imp, covariates) {
  options(datadist = rms::datadist(d))
  fit <- Hmisc::fit.mult.impute(cox_formula(covariates), rms::cph, imp, data = d,
                                fitargs = list(x = TRUE, y = TRUE), pr = FALSE,
                                fun = function(f) c(C = (f$stats[["Dxy"]] + 1) / 2,
                                                    ph = survival::cox.zph(f)$table["GLOBAL", "p"]))
  stopifnot(fit$stats[["Obs"]] == nrow(d))
  by_imputation <- do.call(rbind, fit$funresults)
  cbind(hazard_ratios(fit, covariates), C_index = mean(by_imputation[, "C"]),
        PH_global_p = stats::median(by_imputation[, "ph"]))
}

## Cox model restricted to patients with observed Ki-67.
complete_case_cox <- function(d, covariates) {
  options(datadist = rms::datadist(d))
  fit <- rms::cph(cox_formula(covariates), data = d[!is.na(d$ki67), ], x = TRUE, y = TRUE)
  cbind(hazard_ratios(fit, covariates), C_index = (fit$stats[["Dxy"]] + 1) / 2,
        PH_global_p = survival::cox.zph(fit)$table["GLOBAL", "p"])
}
