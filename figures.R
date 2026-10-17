## Figures 1-2.

## Figure 1. Kaplan-Meier curves for OS (A) and PFS (B).
figure1 <- function(data = load_septralu()) {
  load_dependencies()

  os <- data[!is.na(data$os_time) & !is.na(data$os_event), ]
  fit_os <- survival::survfit(survival::Surv(os_time, os_event) ~ 1, data = os)
  p_os <- survminer::ggsurvplot(
    fit_os, data = os, conf.int = TRUE, risk.table = TRUE,
    xlab = "Months from first R-PRRT cycle", ylab = "Overall survival",
    title = "A. Overall survival", censor.shape = "|",
    palette = "#2E5A87", ggtheme = ggplot2::theme_minimal())

  pfs <- data[!is.na(data$pfs_time) & !is.na(data$pfs_event), ]
  fit_pfs <- survival::survfit(survival::Surv(pfs_time, pfs_event) ~ 1, data = pfs)
  p_pfs <- survminer::ggsurvplot(
    fit_pfs, data = pfs, conf.int = TRUE, risk.table = TRUE,
    xlab = "Months from first R-PRRT cycle", ylab = "Progression-free survival",
    title = "B. Progression-free survival", censor.shape = "|",
    palette = "#8C3B4A", ggtheme = ggplot2::theme_minimal())

  survminer::arrange_ggsurvplots(list(p_os, p_pfs), ncol = 2, nrow = 1,
                                 print = TRUE)
}

## Figure 2. Forest plot of multivariable Cox models for PFS and OS.
figure2 <- function(data = load_septralu(), hr = multivariable_cox(data)) {
  load_dependencies()
  labels <- c(
    ki67 = "Ki-67 index (IQR, 15 vs 2%)",
    ecog = "ECOG PS (per 1-point increase)",
    n_metastatic_sites = "No. of metastatic sites (IQR)",
    pancreas = "Primary site (pancreas vs others)")
  hr$label <- factor(labels[hr$term], levels = rev(labels))
  hr$hrtext <- sprintf("%.2f (%.2f–%.2f)", hr$HR, hr$low, hr$high)

  ggplot2::ggplot(hr, ggplot2::aes(x = HR, y = label)) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_pointrange(ggplot2::aes(xmin = low, xmax = high),
                             colour = "#2E5A87") +
    ggplot2::geom_text(ggplot2::aes(x = high, label = hrtext),
                       hjust = -0.12, size = 2.7, colour = "grey20") +
    ggplot2::facet_wrap(~ Endpoint) +
    ggplot2::scale_x_log10(expand = ggplot2::expansion(mult = c(0.05, 0.6))) +
    ggplot2::labs(x = "Hazard ratio (95% CI)", y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
}
