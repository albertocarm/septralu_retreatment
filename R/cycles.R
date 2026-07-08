## Exploratory analyses of number of R-PRRT cycles and of the I-PRRT to
## R-PRRT interval on survival and response.

km_summary <- function(time, event, group, data) {
  fit <- survival::survfit(survival::Surv(data[[time]], data[[event]]) ~ group,
                           data = data)
  s <- summary(fit)$table
  if (is.null(dim(s))) s <- t(as.matrix(s))
  sd <- survival::survdiff(survival::Surv(data[[time]], data[[event]]) ~ group,
                           data = data)
  p <- 1 - stats::pchisq(sd$chisq, length(sd$n) - 1)
  out <- data.frame(
    Group = sub("group=", "", rownames(s)),
    N = as.integer(s[, "records"]), Events = as.integer(s[, "events"]),
    Median = round(s[, "median"], 1),
    CI_low = round(s[, "0.95LCL"], 1), CI_high = round(s[, "0.95UCL"], 1),
    row.names = NULL)
  out$logrank_p <- round(p, 3)
  out
}

## Number of cycles: primary comparison 2 vs 4 cycles, plus all groups.
cycles_analysis <- function(data = load_septralu()) {
  load_dependencies()
  d <- data[!is.na(data$n_cycles), ]
  d$cycles <- factor(d$n_cycles, levels = 1:4,
                     labels = c("1", "2", "3", "4"))

  two_four <- d[d$n_cycles %in% c(2, 4), ]
  two_four$cycles <- droplevels(two_four$cycles)

  response_by_group <- function(df) {
    ev <- df[!is.na(df$recist), ]
    tab <- tapply(ev$recist, droplevels(ev$cycles), function(r)
      c(n = length(r),
        DCR = round(mean(r %in% c("CR", "PR", "SD")) * 100, 1),
        ORR = round(mean(r %in% c("CR", "PR")) * 100, 1)))
    do.call(rbind, tab)
  }
  chisq_dcr <- function(df) {
    ev <- df[!is.na(df$recist), ]
    dcr <- ev$recist %in% c("CR", "PR", "SD")
    suppressWarnings(stats::chisq.test(table(droplevels(ev$cycles), dcr))$p.value)
  }

  list(
    os_2v4        = km_summary("os_time", "os_event", two_four$cycles, two_four),
    pfs_2v4       = km_summary("pfs_time", "pfs_event", two_four$cycles, two_four),
    os_all        = km_summary("os_time", "os_event", d$cycles, d),
    pfs_all       = km_summary("pfs_time", "pfs_event", d$cycles, d),
    response_2v4  = response_by_group(two_four),
    response_all  = response_by_group(d),
    dcr_p_2v4     = round(chisq_dcr(two_four), 3),
    dcr_p_all     = round(chisq_dcr(d), 3)
  )
}

## Number of cycles modelled as a time-dependent covariate.
## Cumulative cycles received change value at each cycle date, so a patient
## contributes person-time under 1, then 2, ... cycles as they are delivered.
## This removes the guarantee-time bias of the naive baseline comparison.
cycles_timedep <- function(data = load_septralu()) {
  load_dependencies()

  fit_td <- function(time, event) {
    base <- data.frame(id = seq_len(nrow(data)),
                       T = data[[time]], E = data[[event]],
                       t2 = data$cycle2_month, t3 = data$cycle3_month,
                       t4 = data$cycle4_month)
    base <- base[!is.na(base$T) & !is.na(base$E) & base$T > 0, ]
    dt <- survival::tmerge(base, base, id = id, endpt = event(T, E))
    dt <- survival::tmerge(dt, base, id = id, c2 = tdc(t2))
    dt <- survival::tmerge(dt, base, id = id, c3 = tdc(t3))
    dt <- survival::tmerge(dt, base, id = id, c4 = tdc(t4))
    dt$n_cycles_received <- 1 + dt$c2 + dt$c3 + dt$c4
    fit <- survival::coxph(
      survival::Surv(tstart, tstop, endpt) ~ n_cycles_received, data = dt)
    s <- summary(fit)
    c(HR = round(s$conf.int[1, 1], 2), low = round(s$conf.int[1, 3], 2),
      high = round(s$conf.int[1, 4], 2), p = round(s$coefficients[1, 5], 3),
      patients = length(unique(dt$id)), events = sum(dt$endpt))
  }

  fit_naive <- function(time, event) {
    s <- data[data$n_cycles %in% c(2, 4) & !is.na(data[[time]]), ]
    s$g <- factor(ifelse(s$n_cycles == 4, "4", "2"))
    fit <- survival::coxph(survival::Surv(s[[time]], s[[event]]) ~ g, data = s)
    cf <- summary(fit)
    c(HR = round(cf$conf.int[1, 1], 2), low = round(cf$conf.int[1, 3], 2),
      high = round(cf$conf.int[1, 4], 2), p = round(cf$coefficients[1, 5], 3))
  }

  list(
    timedep_os  = fit_td("os_time", "os_event"),
    timedep_pfs = fit_td("pfs_time", "pfs_event"),
    naive_os    = fit_naive("os_time", "os_event"),
    naive_pfs   = fit_naive("pfs_time", "pfs_event")
  )
}

## Systemic antitumour therapy given between the two PRRT courses.
sequence_analysis <- function(data = load_septralu()) {
  load_dependencies()
  km <- function(time, event) {
    d <- data[!is.na(data[[time]]), ]
    fit <- survival::survfit(
      survival::Surv(d[[time]], d[[event]]) ~ intercurrent_systemic_therapy, data = d)
    s <- summary(fit)$table
    sd <- survival::survdiff(
      survival::Surv(d[[time]], d[[event]]) ~ intercurrent_systemic_therapy, data = d)
    data.frame(Group = c("No", "Yes"), N = as.integer(s[, "records"]),
               Events = as.integer(s[, "events"]), Median = round(s[, "median"], 1),
               logrank_p = round(1 - stats::pchisq(sd$chisq, length(sd$n) - 1), 3),
               row.names = NULL)
  }
  list(prevalence = table(data$intercurrent_systemic_therapy),
       os = km("os_time", "os_event"), pfs = km("pfs_time", "pfs_event"))
}

## Interval from last I-PRRT cycle to first R-PRRT cycle.
interval_analysis <- function(data = load_septralu()) {
  load_dependencies()
  d <- data[!is.na(data$retreatment_interval_months), ]

  cutpoint_row <- function(cut, time, event, label) {
    dd <- d[!is.na(d[[time]]), ]
    dd$g <- factor(ifelse(dd$retreatment_interval_months < cut, "short", "long"),
                   levels = c("long", "short"))
    sd <- survival::survdiff(survival::Surv(dd[[time]], dd[[event]]) ~ g, data = dd)
    ci <- summary(survival::coxph(
      survival::Surv(dd[[time]], dd[[event]]) ~ g, data = dd))$conf.int[1, c(1, 3, 4)]
    data.frame(cutpoint = cut, endpoint = label, n_short = sum(dd$g == "short"),
               HR_short_vs_long = round(ci[1], 2), low = round(ci[2], 2),
               high = round(ci[3], 2),
               logrank_p = round(1 - stats::pchisq(sd$chisq, 1), 3), row.names = NULL)
  }
  cutpoints <- do.call(rbind, unlist(lapply(c(12, 18, 24), function(cc)
    list(cutpoint_row(cc, "os_time", "os_event", "OS"),
         cutpoint_row(cc, "pfs_time", "pfs_event", "PFS"))), recursive = FALSE))

  cox_uni <- function(time, event) {
    dd <- rms::datadist(d[, "retreatment_interval_months", drop = FALSE])
    options(datadist = dd); on.exit(options(datadist = NULL))
    fit <- rms::cph(survival::Surv(d[[time]], d[[event]]) ~
                    retreatment_interval_months, data = d)
    s <- summary(fit)
    hr <- s[s[, "Type"] == 2, c("Effect", "Lower 0.95", "Upper 0.95")]
    co <- s[s[, "Type"] == 1, ]
    c(HR = round(hr["Effect"], 2), low = round(hr["Lower 0.95"], 2),
      high = round(hr["Upper 0.95"], 2),
      p = round(2 * stats::pnorm(-abs(co["Effect"] / co["S.E."])), 3))
  }

  cuts <- stats::quantile(d$retreatment_interval_months, c(0, 1/3, 2/3, 1),
                          na.rm = TRUE)
  d$interval_tertile <- cut(d$retreatment_interval_months, breaks = cuts,
                            include.lowest = TRUE,
                            labels = c("Short", "Intermediate", "Long"))

  list(
    median_interval = round(stats::median(d$retreatment_interval_months), 1),
    tertile_cutpoints = round(cuts, 1),
    cox_os  = cox_uni("os_time", "os_event"),
    cox_pfs = cox_uni("pfs_time", "pfs_event"),
    os_by_tertile  = km_summary("os_time", "os_event", d$interval_tertile, d),
    pfs_by_tertile = km_summary("pfs_time", "pfs_event", d$interval_tertile, d),
    cutpoints = cutpoints
  )
}
