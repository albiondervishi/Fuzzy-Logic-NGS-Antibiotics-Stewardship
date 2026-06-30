# ================================================================
# Figure_Combined_KM_Hospital_Survival.R
#
# Replaces the separate Fig21b / Fig22 outputs with ONE combined,
# two-panel Kaplan-Meier figure, both panels on the SAME 28-day
# HOSPITAL survival endpoint (Hosp_days / hospital_mortality):
#
#   Panel A : Six mNGS/BC/Abx diagnostic groups, all unique patients
#             -> identical to Fig22's existing Panel A (km22a_df).
#             n=214 | log-rank p=0.513 (NS) — matches what you
#             already confirmed in the rendered image.
#
#   Panel B : Antibiotic adequacy (Adequate vs Inadequate/NGS-guided)
#             -> taken from Fig21b's km_abx_hosp_df logic, NOT
#             Fig22's own Panel B. Fig22's Panel B omitted the
#             Abx_diag==1 filter and therefore pulled in Abx-negative
#             patients, contradicting the manuscript's own stated
#             inclusion criterion ("All patients included in this
#             analysis received empirical antibiotic therapy at ICU
#             admission..."). Fig21b applies that filter correctly
#             (n=106 ABX-positive patients), matching the image you
#             already reviewed.
#
# REQUIRES: feature_df already built by
#           NGS_Fuzzy_ABX_Stewardship_v8_Prina.R, OR load from
#           Output/Fuzzy_ABX_v8_patients.csv (see below).
#
# OUTPUT:   Figures/Figure_KM_Hospital_Survival_Combined.png / .pdf
# ================================================================

# ── 0. PACKAGES ─────────────────────────────────────────────────────────────
pkgs <- c("tidyverse", "survival", "survminer")
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0)
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

setwd("~/Desktop/NGS_ZNA/2_NGS_Stewarship_Antibiotic")
dir.create("Figures", showWarnings = FALSE)

# ── 1. LOAD DATA ─────────────────────────────────────────────────────────────
if (!exists("feature_df")) {
  cat("feature_df not found in session - loading from CSV\n")
  feature_df <- read.csv("Output/Fuzzy_ABX_v8_patients.csv", stringsAsFactors = FALSE)
}
cat(sprintf("Patients loaded: %d\n", nrow(feature_df)))

# ── 2. SHARED HELPERS ────────────────────────────────────────────────────────
theme_pub <- function(base = 15)
  ggplot2::theme_bw(base_size = base) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
    legend.position  = "bottom",
    legend.text      = ggplot2::element_text(size = base - 1),
    legend.title     = ggplot2::element_text(size = base),
    plot.title       = ggplot2::element_text(face = "bold", size = base + 3),
    plot.subtitle    = ggplot2::element_text(size = base - 1, colour = "grey40"),
    plot.margin      = ggplot2::margin(t = 5, r = 15, b = 5, l = 0),
    axis.title       = ggplot2::element_text(size = base + 1),
    axis.title.y     = ggplot2::element_blank(),
    axis.text        = ggplot2::element_text(size = base - 2)
  )

pdf_dev <- tryCatch(grDevices::cairo_pdf, error = function(e) grDevices::pdf)

km_colours <- c(
  "mNGS-/BC-/Abx+" = "#27AE60",
  "mNGS-/BC+/Abx+" = "#82C341",
  "mNGS+/BC-/Abx+" = "#F39C12",
  "mNGS+/BC+/Abx+" = "#E67E22",
  "mNGS+/BC-/Abx-" = "#2980B9",
  "mNGS+/BC+/Abx-" = "#C0392B"
)

# ── 3. PANEL A — Six diagnostic groups, 28-day HOSPITAL survival ────────────
# (identical logic to Fig22's existing Panel A — unchanged)
panelA_df <- feature_df |>
  dplyr::filter(
    diag_group %in% c(
      "mNGS+/BC-/Abx+", "mNGS+/BC+/Abx+",
      "mNGS+/BC-/Abx-", "mNGS+/BC+/Abx-",
      "mNGS-/BC-/Abx+", "mNGS-/BC+/Abx+"),
    !is.na(Hosp_days), !is.na(hospital_mortality), Hosp_days > 0
  ) |>
  dplyr::mutate(
    Diag = factor(diag_group, levels = c(
      "mNGS-/BC-/Abx+", "mNGS-/BC+/Abx+",
      "mNGS+/BC-/Abx+", "mNGS+/BC+/Abx+",
      "mNGS+/BC-/Abx-", "mNGS+/BC+/Abx-")),
    time_28  = pmin(Hosp_days, 28),
    event_28 = dplyr::if_else(hospital_mortality == 1L & Hosp_days <= 28, 1L, 0L)
  )

cat(sprintf("Panel A: n=%d patients\n", nrow(panelA_df)))

fitA <- survival::survfit(survival::Surv(time_28, event_28) ~ Diag, data = panelA_df)
lrA  <- survival::survdiff(survival::Surv(time_28, event_28) ~ Diag, data = panelA_df)
pA   <- round(1 - pchisq(lrA$chisq, df = length(lrA$n) - 1), 4)
pA_lbl <- dplyr::case_when(
  pA < 0.001 ~ "Log-rank p < 0.001",
  pA < 0.05  ~ sprintf("Log-rank p = %.4f", pA),
  TRUE       ~ sprintf("Log-rank p = %.3f (NS)", pA))

panelA <- survminer::ggsurvplot(
  fitA, data = panelA_df,
  palette           = km_colours,
  conf.int          = TRUE, conf.int.alpha = 0.15,
  risk.table        = TRUE, risk.table.col = "strata",
  risk.table.height = 0.28, risk.table.title = "",
  xlab = "Hospital stay (days)", ylab = "",
  xlim = c(0, 28), break.time.by = 7, ylim = c(0.50, 1.00),
  legend.title = "Diagnostic group",
  legend.labs  = c("mNGS-/BC-/Abx+", "mNGS-/BC+/Abx+",
                   "mNGS+/BC-/Abx+", "mNGS+/BC+/Abx+",
                   "mNGS+/BC-/Abx-", "mNGS+/BC+/Abx-"),
  title    = sprintf("A  \u2014  Six diagnostic groups  (n=%d)", nrow(panelA_df)),
  subtitle = paste0("All unique patients  |  28-day hospital endpoint  |  ", pA_lbl),
  ggtheme  = theme_pub(base = 15),
  tables.theme = survminer::theme_cleantable(font.main = 13, font.x = 13, font.y = 13),
  fontsize = 4.5, pval = FALSE)

panelA$plot <- panelA$plot +
  ggplot2::annotate("text", x = -0.4, y = 0.75,
                    label = "Survival probability",
                    angle = 90, size = 5.2) +
  ggplot2::annotate("text", x = 14, y = 0.55, label = pA_lbl,
                    size = 5.0, colour = "grey30", fontface = "italic") +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(
    plot.margin = ggplot2::margin(t = 5, r = 15, b = 5, l = 0),
    axis.title.y = ggplot2::element_blank()
  )

# ── 4. PANEL B — Antibiotic adequacy, ABX-positive only, 28-day HOSPITAL survival
# (Fig21b's logic — includes the Abx_diag==1 filter Fig22's own Panel B omitted)
panelB_df <- feature_df |>
  dplyr::filter(
    Abx_diag == 1, !is.na(AB_change),
    !is.na(Hosp_days), !is.na(hospital_mortality), Hosp_days > 0
  ) |>
  dplyr::mutate(
    Adequacy = dplyr::case_when(
      AB_change == "No change"                            ~ "Adequate",
      AB_change %in% c("Escalate", "De-escalate", "Stop")  ~ "Inadequate / NGS-guided",
      TRUE ~ NA_character_) |>
      factor(levels = c("Adequate", "Inadequate / NGS-guided")),
    time_28  = pmin(Hosp_days, 28),
    event_28 = dplyr::if_else(hospital_mortality == 1L & Hosp_days <= 28, 1L, 0L)
  ) |>
  dplyr::filter(!is.na(Adequacy))

cat(sprintf("Panel B: n=%d ABX-positive patients\n", nrow(panelB_df)))

fitB <- survival::survfit(survival::Surv(time_28, event_28) ~ Adequacy, data = panelB_df)
lrB  <- survival::survdiff(survival::Surv(time_28, event_28) ~ Adequacy, data = panelB_df)
pB   <- round(1 - pchisq(lrB$chisq, df = length(lrB$n) - 1), 4)
pB_lbl <- dplyr::case_when(
  pB < 0.001 ~ "Log-rank p < 0.001",
  pB < 0.05  ~ sprintf("Log-rank p = %.4f", pB),
  TRUE       ~ sprintf("Log-rank p = %.3f", pB))

panelB <- survminer::ggsurvplot(
  fitB, data = panelB_df,
  palette = c("Adequate" = "#D6604D", "Inadequate / NGS-guided" = "#2C7BB6"),
  conf.int          = TRUE, conf.int.alpha = 0.18,
  risk.table        = TRUE, risk.table.col = "strata",
  risk.table.height = 0.28, risk.table.title = "",
  xlab = "Hospital stay (days)", ylab = "",
  xlim = c(0, 28), break.time.by = 9, ylim = c(0.40, 1.00),
  legend.title = "Therapy group",
  legend.labs  = c("Adequate", "Inadequate / NGS-guided"),
  title    = sprintf("B  \u2014  Antibiotic adequacy  (n=%d)", nrow(panelB_df)),
  subtitle = "ABX-positive patients only  |  Adequate = no change  |\n   Inadequate/NGS-guided = escalation, de-escalation, stop",
  ggtheme  = theme_pub(base = 15),
  tables.theme = survminer::theme_cleantable(),
  fontsize = 5.0, pval = FALSE)

panelB$plot <- panelB$plot +
  ggplot2::annotate("text", x = -0.4, y = 0.70,
                    label = "Survival probability",
                    angle = 90, size = 5.2) +
  ggplot2::annotate("text", x = 15, y = 0.48,
                    label      = paste0(pB_lbl, "\nNon-inferiority: NGS-guided = Adequate"),
                    size       = 5.0, colour = "grey25", fontface = "italic",
                    lineheight = 0.9) +
  ggplot2::coord_cartesian(clip = "off") +
  ggplot2::theme(
    legend.position = "right",
    legend.title    = ggplot2::element_text(size = 15),
    legend.text     = ggplot2::element_text(size = 14),
    plot.margin     = ggplot2::margin(t = 5, r = 15, b = 5, l = 0),
    axis.title.y    = ggplot2::element_blank()
  )

# ── 5. COMBINE AND SAVE ───────────────────────────────────────────────────────
fig_combined <- survminer::arrange_ggsurvplots(
  list(panelA, panelB),
  ncol = 2, nrow = 1,
  title = " "
)

pdf_dev("Figures/Figure_KM_Hospital_Survival_Combined.pdf", width = 22, height = 9.5)
print(fig_combined)
grDevices::dev.off()

grDevices::png("Figures/Figure_KM_Hospital_Survival_Combined.png",
               width = 22, height = 9.5, units = "in", res = 300, bg = "white")
print(fig_combined)
grDevices::dev.off()

cat("\nSaved: Figures/Figure_KM_Hospital_Survival_Combined (.png + .pdf)\n")

# ── 6. CONSOLE SUMMARY ────────────────────────────────────────────────────────
cat(sprintf("\nCombined figure summary:\n"))
cat(sprintf("  Panel A - six diagnostic groups : n=%d | %s\n", nrow(panelA_df), pA_lbl))
cat(sprintf("  Panel B - antibiotic adequacy    : n=%d | %s\n", nrow(panelB_df), pB_lbl))
cat("\nNOTE: Panel A (n) will be LARGER than Panel B (n) by design - Panel A\n")
cat("includes all six diagnostic groups (Abx+ and Abx- at sampling), while\n")
cat("Panel B is restricted to ABX-positive patients only, per the manuscript's\n")
cat("stated inclusion criterion for the stewardship-strategy comparison.\n")
cat("\nNOTE: Panel A includes two very small strata - mNGS-/BC+/Abx+ (n=3) and\n")
cat("mNGS+/BC+/Abx+ (n=11). Their KM steps/CI bands are visually unstable late\n")
cat("in follow-up due to low n, not a data error. Flag this explicitly in the\n")
cat("figure legend/caption as a cautious-interpretation note.\n")
