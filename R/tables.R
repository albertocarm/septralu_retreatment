## Tables 1-3.

## Table 1. Baseline characteristics, stratified by primary tumour site.
table1 <- function(data = load_septralu()) {
  load_dependencies()
  vars <- c("age", "sex", "ki67", "grade", "ecog", "functioning_tumor",
            "peritoneal_mets", "liver_mets", "n_metastatic_sites",
            "pet_ga_heterogeneity",
            "metastasis_to_retreatment_months", "retreatment_interval_months")
  d <- data[, c("primary_site", vars)]
  d$primary_site <- factor(d$primary_site, levels = c("Pancreas", "Ileum", "Lung", "Other"))
  gtsummary::tbl_summary(
    d,
    by = primary_site,
    label = variable_labels[vars],
    type = list(ki67 ~ "continuous",
                n_metastatic_sites ~ "continuous",
                retreatment_interval_months ~ "continuous",
                metastasis_to_retreatment_months ~ "continuous"),
    statistic = list(gtsummary::all_continuous() ~ "{median} [{p25}, {p75}]",
                     gtsummary::all_categorical() ~ "{n} ({p}%)"),
    digits = gtsummary::all_continuous() ~ 1,
    missing = "ifany"
  ) |>
    gtsummary::add_overall() |>
    gtsummary::modify_header(label ~ "**Characteristic**") |>
    gtsummary::bold_labels()
}

## Table 2. Treatment response.
table2 <- function(data = load_septralu()) {
  load_dependencies()
  vars <- c("recist", "clinical_response",
            "biochemical_response_cga", "biochemical_response_secreted")
  tbl <- gtsummary::tbl_summary(
    data[, vars],
    label = variable_labels[vars],
    statistic = gtsummary::all_categorical() ~ "{n} ({p}%)",
    missing = "no"
  ) |>
    gtsummary::modify_header(label ~ "**Response**") |>
    gtsummary::bold_labels()

  ev  <- data$recist[!is.na(data$recist)]
  dcr <- mean(ev %in% c("CR", "PR", "SD")) * 100
  orr <- mean(ev %in% c("CR", "PR")) * 100
  message(sprintf("Evaluable n = %d | DCR = %.1f%% | ORR = %.1f%%",
                  length(ev), dcr, orr))
  tbl
}

## Table 3. Kaplan-Meier estimates for OS and PFS, overall and by primary site.
table3 <- function(data = load_septralu()) {
  load_dependencies()

  km_row <- function(fit, group) {
    s <- summary(fit)$table
    if (is.null(dim(s))) s <- t(as.matrix(s))
    data.frame(
      Group   = group,
      N       = as.integer(s[, "records"]),
      Events  = as.integer(s[, "events"]),
      Median  = round(s[, "median"], 1),
      CI_low  = round(s[, "0.95LCL"], 1),
      CI_high = round(s[, "0.95UCL"], 1),
      row.names = NULL
    )
  }

  endpoint <- function(time, event, label) {
    ok <- !is.na(data[[time]]) & !is.na(data[[event]])
    d  <- data[ok, ]
    d$site <- addNA(d$primary_site)
    levels(d$site)[is.na(levels(d$site))] <- "Unknown primary"
    overall <- km_row(survival::survfit(
      survival::Surv(d[[time]], d[[event]]) ~ 1), "Overall")
    by_site <- km_row(survival::survfit(
      survival::Surv(d[[time]], d[[event]]) ~ site, data = d),
      levels(droplevels(d$site)))
    p <- survival::survdiff(
      survival::Surv(d[[time]], d[[event]]) ~ site, data = d)
    p <- 1 - stats::pchisq(p$chisq, length(p$n) - 1)
    res <- rbind(overall, by_site)
    res$Endpoint <- label
    res$logrank_p <- c(NA, rep(round(p, 3), nrow(by_site)))
    res
  }

  rbind(
    endpoint("os_time",  "os_event",  "Overall survival"),
    endpoint("pfs_time", "pfs_event", "Progression-free survival")
  )[, c("Endpoint", "Group", "N", "Events", "Median", "CI_low", "CI_high", "logrank_p")]
}
