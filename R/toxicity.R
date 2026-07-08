## Maximum toxicity grade per patient across all R-PRRT cycles (CTCAE).
toxicity_table <- function(data = load_septralu()) {
  d <- data[data$tox_evaluable %in% TRUE, ]
  types <- c(tox_nausea = "Nausea", tox_vomiting = "Vomiting",
             tox_urticaria = "Urticaria", tox_hematologic = "Haematologic",
             tox_nephrotoxicity = "Nephrotoxicity", tox_other = "Other",
             tox_overall_max = "Overall maximum")
  grid <- 0:4
  res <- do.call(rbind, lapply(names(types), function(col) {
    g <- factor(d[[col]], levels = grid)
    n <- table(g)
    row <- data.frame(Toxicity = types[[col]])
    for (k in grid) {
      cnt <- as.integer(n[as.character(k)])
      row[[paste0("Grade", k)]] <- sprintf("%d (%.1f%%)", cnt, 100 * cnt / sum(n))
    }
    row
  }))
  rownames(res) <- NULL
  attr(res, "n_evaluable") <- nrow(d)
  res
}

## Characteristics of patients who developed a therapy-related myeloid neoplasm.
myeloid_neoplasm_cases <- function(data = load_septralu()) {
  cases <- data[data$therapy_related_myeloid_neoplasm == "Yes", ]
  cases <- cases[order(cases$cumulative_activity_gbq), ]
  list(
    n = nrow(cases),
    percent = round(100 * nrow(cases) / nrow(data), 1),
    without_alkylating = sum(cases$alkylating_exposure == "No", na.rm = TRUE),
    cumulative_activity_range = range(cases$cumulative_activity_gbq, na.rm = TRUE),
    cases = data.frame(
      age = cases$age, sex = cases$sex, primary_site = cases$primary_site,
      ki67 = cases$ki67, iprrt_cycles = cases$iprrt_cycles,
      rprrt_cycles = cases$n_cycles,
      cumulative_activity_gbq = cases$cumulative_activity_gbq,
      alkylating_exposure = cases$alkylating_exposure)
  )
}
