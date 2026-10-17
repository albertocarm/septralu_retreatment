## Sensitivity analyses for the multivariable Cox models (Online Resource 7):
## the primary model (multiple imputation of Ki-67), a complete-case analysis,
## and the primary model with ECOG dichotomised (>= 1 vs 0). A further analysis
## excludes patients with an unknown primary site.
sensitivity_table <- function(data = load_septralu(), n_impute = 50, seed = 2026) {
  load_dependencies()
  ecog_binary <- replace(cox_covariates, cox_covariates == "ecog", "ecog_1plus")
  res <- list()
  for (endpoint in c("os", "pfs")) {
    d <- model_data(data, endpoint)
    imp <- impute_ki67(d, n_impute, seed)
    known <- model_data(data[!is.na(data$primary_site), ], endpoint)
    res[[endpoint]] <- list(
      primary       = imputed_cox(d, imp, cox_covariates),
      complete_case = complete_case_cox(d, cox_covariates),
      ecog_binary   = imputed_cox(d, imp, ecog_binary),
      known_primary = imputed_cox(known, impute_ki67(known, n_impute, seed),
                                  cox_covariates)
    )
  }
  res
}

print_sensitivity_table <- function(res) {
  for (endpoint in names(res)) {
    cat("\n== ", toupper(endpoint), " ==\n", sep = "")
    for (model in names(res[[endpoint]])) {
      r <- res[[endpoint]][[model]]
      cat(sprintf("\n%s (n = %d, events = %d, C = %.3f, PH p = %.2f)\n", model,
                  r$n[1], r$events[1], r$C_index[1], r$PH_global_p[1]))
      p <- ifelse(r$p < 0.001, "< 0.001", sprintf("= %.3f", r$p))
      cat(sprintf("  %-20s %.2f (%.2f-%.2f); p %s\n",
                  r$term, r$HR, r$low, r$high, p), sep = "")
    }
  }
  invisible(res)
}
