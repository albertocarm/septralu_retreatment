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
