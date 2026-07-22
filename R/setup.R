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
  discontinuation_reason            = "Reason for discontinuation",
  recist                            = "RECIST 1.1 response",
  clinical_response                 = "Clinical response",
  biochemical_response_cga          = "Biochemical response (chromogranin A)",
  biochemical_response_secreted     = "Biochemical response (secreted substance)"
)

## Derived model covariates used in the multivariable Cox models.
## ECOG is modelled as an ordinal score (per 1-point increase). A 0-1 vs 2+
## dichotomy is uninformative in this cohort, where only 6 patients have a
## performance status above 1 and none above 2.
add_model_covariates <- function(data) {
  data$ki67_imputed <- ifelse(is.na(data$ki67),
                              mean(data$ki67, na.rm = TRUE), data$ki67)
  data$ecog_linear <- as.numeric(as.character(data$ecog))
  data
}

## Parsimonious multivariable model (appropriate for the number of events):
## proliferation (Ki-67), performance status, disease burden (metastatic sites)
## and primary tumour site.
cox_covariates <- c(
  "ki67_imputed", "ecog_linear", "n_metastatic_sites", "primary_site_pancreas"
)

## Reference distribution for the multivariable models. Continuous covariates are
## contrasted over their interquartile range; ECOG is contrasted over a single
## point, which is the clinically interpretable increment.
model_datadist <- function(data) {
  dd <- rms::datadist(data[, cox_covariates])
  dd$limits["Low:effect",  "ecog_linear"] <- 0
  dd$limits["High:effect", "ecog_linear"] <- 1
  dd
}
