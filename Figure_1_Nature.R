# ── Fig 2 · Nature-style marker-specific thresholds ──────────────────────────

nature_theme <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "grey25"),
      axis.ticks = element_line(linewidth = 0.3, colour = "grey25"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base),
      plot.title = element_text(face = "bold", size = base + 3),
      plot.subtitle = element_text(size = base, colour = "grey35"),
      legend.position = "bottom",
      legend.title = element_blank(),
      panel.spacing = unit(1.2, "lines")
    )
}

col_response <- "#2E7D32"
col_alarm    <- "#C62828"
col_neutral  <- "grey35"
col_auc      <- "black"

# A — severity membership curves
fig2A <- sev_curves |>
  tidyr::pivot_longer(
    cols = c(Normal, Alarm_high, Alarm_low),
    names_to = "Membership",
    values_to = "mu"
  ) |>
  dplyr::filter(!(Membership == "Alarm_low" & mu == 0)) |>
  dplyr::mutate(
    Membership = dplyr::recode(
      Membership,
      Normal = "Physiological / low-risk",
      Alarm_high = "High-severity alarm",
      Alarm_low = "Low-severity alarm"
    )
  ) |>
  ggplot(aes(x = x, y = mu, colour = Membership)) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~Marker, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(
    "Physiological / low-risk" = col_response,
    "High-severity alarm" = col_alarm,
    "Low-severity alarm" = col_alarm
  )) +
  scale_y_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.5, 1),
    labels = scales::percent_format()
  ) +
  labs(
    title = "A  Absolute severity membership",
    x = "Biomarker value",
    y = expression(Membership~mu)
  ) +
  nature_theme(8)

# B — kinetic threshold map
kinetic_bars2 <- threshold_validation |>
  dplyr::mutate(
    Marker = factor(
      Marker,
      levels = rev(c("CRP", "PCT", "Leukocytes", "Temperature", "SAPS II"))
    ),
    mild_resp_val    = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_resp),
    strong_resp_val  = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_resp),
    mild_alarm_val   = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_alarm),
    strong_alarm_val = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_alarm),
    auc_lab = dplyr::if_else(
      is.na(AUC),
      "",
      sprintf("AUC %.2f", AUC)
    )
  )

fig2B <- ggplot(kinetic_bars2, aes(y = Marker)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "grey45") +
  
  geom_segment(
    aes(x = strong_resp_val, xend = mild_resp_val, yend = Marker),
    linewidth = 7, alpha = 0.65, colour = col_response,
    lineend = "round"
  ) +
  
  geom_segment(
    aes(x = mild_alarm_val, xend = strong_alarm_val, yend = Marker),
    linewidth = 7, alpha = 0.65, colour = col_alarm,
    lineend = "round"
  ) +
  
  geom_segment(
    data = tibble::tibble(
      Marker = factor("Leukocytes", levels = levels(kinetic_bars2$Marker)),
      x = -100,
      xend = -60
    ),
    aes(x = x, xend = xend, y = Marker, yend = Marker),
    linewidth = 7, alpha = 0.65, colour = col_alarm,
    lineend = "round",
    inherit.aes = FALSE
  ) +
  
  geom_point(
    aes(x = optimal_cut),
    shape = 23, size = 3.5, fill = "white", colour = col_auc, stroke = 0.8
  ) +
  
  geom_text(
    aes(x = optimal_cut, label = auc_lab),
    nudge_y = 0.28, size = 2.7, colour = "grey20"
  ) +
  
  annotate("text", x = -42, y = 0.55, label = "Response zone",
           colour = col_response, fontface = "bold", size = 3) +
  annotate("text", x = 42, y = 0.55, label = "Alarm zone",
           colour = col_alarm, fontface = "bold", size = 3) +
  
  labs(
    title = "B  Marker-specific kinetic response and alarm thresholds",
    subtitle = "Bars indicate rule-defined zones; diamonds show cohort-derived Youden cut-points",
    x = expression(Delta~threshold~"(% or °C for temperature)"),
    y = NULL
  ) +
  nature_theme(9)

# C — compact validation table
fig2C_data <- threshold_validation |>
  dplyr::transmute(
    Marker,
    Metric = dplyr::if_else(metric == "delta_abs", "Δ°C", "Δ%"),
    `Response` = paste0(
      purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_resp),
      " to ",
      purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_resp)
    ),
    `Alarm` = paste0(
      purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_alarm),
      " to ",
      purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_alarm)
    ),
    `AUC [95% CI]` = dplyr::if_else(
      is.na(AUC),
      "n/a",
      sprintf("%.2f [%.2f–%.2f]", AUC, AUC_lower, AUC_upper)
    )
  )

fig2C <- gridExtra::tableGrob(
  fig2C_data,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    core = list(
      fg_params = list(cex = 0.72),
      padding = unit(c(3, 4), "mm")
    ),
    colhead = list(
      fg_params = list(cex = 0.75, fontface = "bold"),
      padding = unit(c(3, 4), "mm")
    )
  )
)

fig2_final <- (fig2A / fig2B / patchwork::wrap_elements(fig2C)) +
  plot_layout(heights = c(1.05, 1.15, 0.65)) +
  plot_annotation(
    title = "Figure 2 | Marker-specific kinetic thresholds",
    subtitle = "Severity membership, kinetic response/alarm zones, and cohort-derived validation",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15),
      plot.subtitle = element_text(size = 10, colour = "grey35")
    )
  )

save_fig(fig2_final, "Fig2_Nature_marker_thresholds", 14, 11)












# Packages needed
library(tidyverse)
library(patchwork)
library(gridExtra)
library(cowplot)
library(grid)

# ── Figure 2 design settings ────────────────────────────────────────────────
blue  <- "#003B73"
green <- "#1B7F32"
red   <- "#D7191C"
dark  <- "#111111"

theme_big <- function(base = 15) {
  theme_classic(base_size = base) +
    theme(
      plot.title = element_text(face = "bold", size = base + 4),
      plot.subtitle = element_text(size = base, colour = "grey35"),
      axis.title = element_text(face = "bold", size = base),
      axis.text = element_text(size = base - 1),
      strip.text = element_text(face = "bold", size = base),
      legend.text = element_text(size = base - 1),
      legend.position = "bottom"
    )
}

# ── Panel A: absolute severity curves ───────────────────────────────────────
figA_new <- sev_curves |>
  pivot_longer(
    c(Normal, Alarm_high, Alarm_low),
    names_to = "Zone",
    values_to = "mu"
  ) |>
  filter(!(Zone == "Alarm_low" & mu == 0)) |>
  mutate(
    Zone = recode(
      Zone,
      Normal = "Physiological / low-risk",
      Alarm_high = "High-severity alarm",
      Alarm_low = "Low-severity alarm"
    )
  ) |>
  ggplot(aes(x = x, y = mu, colour = Zone)) +
  geom_line(linewidth = 1.3) +
  facet_wrap(~Marker, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(
    "Physiological / low-risk" = green,
    "High-severity alarm" = red,
    "Low-severity alarm" = red
  )) +
  scale_y_continuous(
    labels = scales::percent_format(),
    breaks = c(0, 0.5, 1),
    limits = c(0, 1.03)
  ) +
  labs(x = "Biomarker value", y = expression(Membership~mu)) +
  theme_big(14)


figA_new <- figA_new +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    plot.margin = margin(
      t = 8,
      r = 20,
      b = 12,
      l = 20
    )
  )

# ── Panel B: kinetic threshold bars ─────────────────────────────────────────
dash_df <- threshold_validation |>
  mutate(
    Marker = factor(
      Marker,
      levels = rev(c("CRP", "PCT", "Leukocytes", "Temperature", "SAPS II"))
    ),
    response_start = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_resp),
    response_end   = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_resp),
    alarm_start    = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_alarm),
    alarm_end      = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_alarm),
    auc_label = if_else(is.na(AUC), "AUC n/a", sprintf("AUC %.2f", AUC))
  )

leuk_low_alarm <- tibble(
  Marker = factor("Leukocytes", levels = levels(dash_df$Marker)),
  x = -100,
  xend = -60
)

figB_new <- ggplot(dash_df, aes(y = Marker)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.7, colour = "grey45") +
  
  geom_segment(
    aes(x = response_start, xend = response_end, yend = Marker),
    linewidth = 13, colour = green, alpha = 0.65, lineend = "round"
  ) +
  geom_segment(
    aes(x = alarm_start, xend = alarm_end, yend = Marker),
    linewidth = 13, colour = red, alpha = 0.65, lineend = "round"
  ) +
  geom_segment(
    data = leuk_low_alarm,
    aes(x = x, xend = xend, y = Marker, yend = Marker),
    inherit.aes = FALSE,
    linewidth = 13, colour = red, alpha = 0.65, lineend = "round"
  ) +
  geom_point(
    aes(x = optimal_cut),
    shape = 23, size = 5.5, fill = "white", colour = dark, stroke = 1.1
  ) +
  geom_text(
    aes(x = 125, label = auc_label),
    size = 5,
    hjust = 0,
    colour = "grey15"
  ) +
  annotate("text", x = -45, y = 0.35, label = "RESPONSE ZONE",
           colour = green, fontface = "bold", size = 5) +
  annotate("text", x = 55, y = 0.35, label = "ALARM ZONE",
           colour = red, fontface = "bold", size = 5) +
  scale_x_continuous(limits = c(-110, 155), breaks = c(-100, -50, 0, 50, 100, 150)) +
  labs(
    x = expression(Delta~threshold~"(% or °C for temperature)"),
    y = NULL
  ) +
  theme_big(15)

figB_new <- figB_new +
  coord_cartesian(clip = "off") +
  theme(
    plot.margin = margin(
      t = 10,
      r = 35,
      b = 35,
      l = 35
    )
  )

# Larger AUC labels
figB_new <- figB_new +
  geom_text(
    data = dash_df,
    aes(x = 125, label = auc_label),
    size = 5.8,
    fontface = "bold",
    hjust = 0,
    inherit.aes = FALSE
  )

# Larger diamonds
figB_new <- figB_new +
  geom_point(
    data = dash_df,
    aes(x = optimal_cut, y = Marker),
    shape = 23,
    size = 7,
    stroke = 1.4,
    fill = "white",
    colour = "black",
    inherit.aes = FALSE
  )


# ── Panel C: clean validation table ─────────────────────────────────────────
table_df <- threshold_validation |>
  transmute(
    Marker,
    Metric = if_else(metric == "delta_abs", "Δ°C", "Δ%"),
    `Response` = paste0(
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_resp),
      " to ",
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_resp)
    ),
    `Alarm` = paste0(
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_alarm),
      " to ",
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_alarm)
    ),
    `AUC [95% CI]` = if_else(
      is.na(AUC),
      "n/a",
      sprintf("%.2f [%.2f–%.2f]", AUC, AUC_lower, AUC_upper)
    ),
    Source = case_when(
      str_detect(source_resp, "Literature") & str_detect(source_alarm, "Literature") ~ "Literature",
      str_detect(source_resp, "Data") | str_detect(source_alarm, "Data") ~ "Cohort-derived / hybrid",
      TRUE ~ "Hybrid"
    )
  )

table_grob <- tableGrob(
  table_df,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(
      fg_params = list(cex = 1.05),
      padding = unit(c(5, 6), "mm")
    ),
    colhead = list(
      fg_params = list(cex = 1.1, fontface = "bold", col = "white"),
      bg_params = list(fill = blue, col = "white"),
      padding = unit(c(5, 6), "mm")
    )
  )
)

figC_new <- cowplot::ggdraw(table_grob)

# ── Helper: boxed panel with blue header ─────────────────────────────────────
boxed_panel <- function(plot, label, title) {
  cowplot::ggdraw() +
    cowplot::draw_grob(
      grid::rectGrob(
        x = 0.5, y = 0.93,
        width = 1, height = 0.14,
        gp = grid::gpar(fill = blue, col = blue)
      )
    ) +
    cowplot::draw_plot(plot, x = 0.02, y = 0.02, width = 0.96, height = 0.82) +
    cowplot::draw_label(
      label, x = 0.025, y = 0.93,
      colour = "white", fontface = "bold",
      size = 22, hjust = 0
    ) +
    cowplot::draw_label(
      title, x = 0.075, y = 0.93,
      colour = "white", fontface = "bold",
      size = 18, hjust = 0
    )
}

panelA <- boxed_panel(figA_new, "A", "ABSOLUTE SEVERITY MEMBERSHIP")
panelB <- boxed_panel(figB_new, "B", "MARKER-SPECIFIC KINETIC RESPONSE AND ALARM THRESHOLDS")
panelC <- boxed_panel(figC_new, "C", "THRESHOLD SUMMARY AND VALIDATION")

# ── Final composite ─────────────────────────────────────────────────────────
panelA <- boxed_panel(
  figA_new,
  "A",
  "ABSOLUTE SEVERITY MEMBERSHIP"
)

panelB <- boxed_panel(
  figB_new,
  "B",
  "MARKER-SPECIFIC KINETIC RESPONSE AND ALARM THRESHOLDS"
)

panelC <- boxed_panel(
  figC_new,
  "C",
  "THRESHOLD SUMMARY AND VALIDATION"
)

# slightly taller Panel C
fig2_new_design <- plot_grid(
  panelA,
  panelB,
  panelC,
  ncol = 1,
  rel_heights = c(
    1.25,
    1.25,
    0.95
  )
)

# ──────────────────────────────────────────────────────────────
# FINAL FIGURE
# ──────────────────────────────────────────────────────────────

final_fig2 <- plot_grid(
  title,
  fig2_new_design,
  ncol = 1,
  rel_heights = c(
    0.12,
    1
  )
)

# ──────────────────────────────────────────────────────────────
# EXPORT
# ──────────────────────────────────────────────────────────────

ggsave(
  "Figures/Fig2_new_infographic_design.png",
  final_fig2,
  width = 16,
  height = 11.5,
  dpi = 400,
  bg = "white"
)

ggsave(
  "Figures/Fig2_new_infographic_design.pdf",
  final_fig2,
  width = 16,
  height = 11.5,
  device = cairo_pdf,
  bg = "white"
)






# Packages needed
library(tidyverse)
library(patchwork)
library(gridExtra)
library(cowplot)
library(grid)

# ── Figure 2 design settings ────────────────────────────────────────────────
blue  <- "#003B73"
green <- "#1B7F32"
red   <- "#D7191C"
dark  <- "#111111"

theme_big <- function(base = 15) {
  theme_classic(base_size = base) +
    theme(
      plot.title = element_text(face = "bold", size = base + 4),
      plot.subtitle = element_text(size = base, colour = "grey35"),
      axis.title = element_text(face = "bold", size = base),
      axis.text = element_text(size = base - 1),
      strip.text = element_text(face = "bold", size = base),
      legend.text = element_text(size = base - 1),
      legend.title = element_blank(),
      legend.position = "bottom",
      plot.margin = margin(8, 20, 12, 20)
    )
}

# ── Panel A: absolute severity curves, 2-zone model ─────────────────────────
marker_units <- c(
  "CRP"         = "CRP (mg/L)",
  "Leukocytes" = "Leukocytes (10^9/L)",
  "PCT"         = "PCT (ng/mL)",
  "SAPS II"     = "SAPS II (points)",
  "Temperature" = "Temperature (°C)"
)

sev_curves_A <- sev_curves |>
  transmute(
    Marker,
    x,
    `Physiological / low-risk` = Normal,
    `High-risk / alarm` = pmax(Alarm_high, Alarm_low, na.rm = TRUE)
  ) |>
  mutate(
    Marker = factor(
      Marker,
      levels = c("CRP", "Leukocytes", "PCT", "SAPS II", "Temperature"),
      labels = marker_units[c("CRP", "Leukocytes", "PCT", "SAPS II", "Temperature")]
    )
  ) |>
  pivot_longer(
    cols = c(`Physiological / low-risk`, `High-risk / alarm`),
    names_to = "Zone",
    values_to = "mu"
  )

figA_new <- ggplot(sev_curves_A, aes(x = x, y = mu, colour = Zone)) +
  geom_line(linewidth = 1.35, na.rm = TRUE) +
  facet_wrap(~Marker, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = c(
    "Physiological / low-risk" = green,
    "High-risk / alarm" = red
  )) +
  scale_y_continuous(
    labels = scales::percent_format(),
    breaks = c(0, 0.5, 1),
    limits = c(0, 1.03)
  ) +
  labs(
    x = NULL,
    y = expression(Membership~mu),
    colour = "Zone"
  ) +
  theme_big(14) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.key.width = unit(1.4, "cm"),
    plot.margin = margin(8, 20, 12, 20)
  )

# ── Panel B ─────────────────────────────────────────────────────────────────
dash_df <- threshold_validation |>
  mutate(
    Marker = factor(
      Marker,
      levels = rev(c("CRP", "PCT", "Leukocytes", "Temperature", "SAPS II"))
    ),
    response_start = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_resp),
    response_end   = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_resp),
    alarm_start    = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_alarm),
    alarm_end      = map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_alarm),
    auc_label = if_else(is.na(AUC), "AUC n/a", sprintf("AUC %.2f", AUC))
  )

leuk_low_alarm <- tibble(
  Marker = factor("Leukocytes", levels = levels(dash_df$Marker)),
  x = -100,
  xend = -60
)

figB_new <- ggplot(dash_df, aes(y = Marker)) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    colour = "grey45"
  ) +
  geom_segment(
    aes(x = response_start, xend = response_end, yend = Marker),
    linewidth = 13,
    colour = green,
    alpha = 0.65,
    lineend = "round"
  ) +
  geom_segment(
    aes(x = alarm_start, xend = alarm_end, yend = Marker),
    linewidth = 13,
    colour = red,
    alpha = 0.65,
    lineend = "round"
  ) +
  geom_segment(
    data = leuk_low_alarm,
    aes(x = x, xend = xend, y = Marker, yend = Marker),
    inherit.aes = FALSE,
    linewidth = 13,
    colour = red,
    alpha = 0.65,
    lineend = "round"
  ) +
  geom_point(
    aes(x = optimal_cut),
    shape = 23,
    size = 7,
    stroke = 1.4,
    fill = "white",
    colour = dark
  ) +
  geom_text(
    aes(x = 125, label = auc_label),
    size = 5.8,
    fontface = "bold",
    hjust = 0,
    colour = "grey15"
  ) +
  annotate(
    "text",
    x = -45,
    y = 0.30,
    label = "RESPONSE ZONE",
    colour = green,
    fontface = "bold",
    size = 5
  ) +
  annotate(
    "text",
    x = 55,
    y = 0.30,
    label = "ALARM ZONE",
    colour = red,
    fontface = "bold",
    size = 5
  ) +
  scale_x_continuous(
    limits = c(-110, 155),
    breaks = c(-100, -50, 0, 50, 100, 150)
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = expression(Delta~threshold~"(% or °C for temperature)"),
    y = NULL
  ) +
  theme_big(15) +
  theme(
    plot.margin = margin(10, 35, 45, 35)
  )

# ── Panel C ─────────────────────────────────────────────────────────────────
table_df <- threshold_validation |>
  transmute(
    Marker,
    Metric = if_else(metric == "delta_abs", "Δ°C", "Δ%"),
    Response = paste0(
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_resp),
      " to ",
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_resp)
    ),
    Alarm = paste0(
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_alarm),
      " to ",
      map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_alarm)
    ),
    `AUC [95% CI]` = if_else(
      is.na(AUC),
      "n/a",
      sprintf("%.2f [%.2f–%.2f]", AUC, AUC_lower, AUC_upper)
    )
  )

table_grob <- tableGrob(
  table_df,
  rows = NULL,
  theme = ttheme_minimal(
    core = list(
      fg_params = list(cex = 1.08),
      padding = unit(c(5, 6), "mm")
    ),
    colhead = list(
      fg_params = list(cex = 1.12, fontface = "bold", col = "white"),
      bg_params = list(fill = blue, col = "white"),
      padding = unit(c(5, 6), "mm")
    )
  )
)

figC_new <- cowplot::ggdraw(table_grob)

# ── Boxed panel helper ──────────────────────────────────────────────────────
boxed_panel <- function(plot, label, title) {
  cowplot::ggdraw() +
    cowplot::draw_grob(
      grid::rectGrob(
        x = 0.5,
        y = 0.94,
        width = 1,
        height = 0.13,
        gp = grid::gpar(fill = blue, col = blue)
      )
    ) +
    cowplot::draw_plot(
      plot,
      x = 0.02,
      y = 0.03,
      width = 0.96,
      height = 0.82
    ) +
    cowplot::draw_label(
      label,
      x = 0.025,
      y = 0.94,
      colour = "white",
      fontface = "bold",
      size = 22,
      hjust = 0
    ) +
    cowplot::draw_label(
      title,
      x = 0.075,
      y = 0.94,
      colour = "white",
      fontface = "bold",
      size = 18,
      hjust = 0
    )
}

panelA <- boxed_panel(figA_new, "A", "ABSOLUTE SEVERITY MEMBERSHIP")
panelB <- boxed_panel(figB_new, "B", "MARKER-SPECIFIC KINETIC RESPONSE AND ALARM THRESHOLDS")
panelC <- boxed_panel(figC_new, "C", "THRESHOLD SUMMARY AND VALIDATION")

# ── Main title ──────────────────────────────────────────────────────────────
title <- ggdraw() +
  draw_label(
    "Figure 2 | Marker-specific kinetic thresholds",
    x = 0,
    y = 0.72,
    hjust = 0,
    fontface = "bold",
    size = 26
  ) +
  draw_label(
    "Severity membership, kinetic response/alarm zones, and cohort-derived validation",
    x = 0,
    y = 0.25,
    hjust = 0,
    size = 16,
    colour = "grey35"
  )

# ── Final composite ─────────────────────────────────────────────────────────
fig2_new_design <- plot_grid(
  panelA,
  panelB,
  panelC,
  ncol = 1,
  rel_heights = c(1.25, 1.25, 0.85)
)

final_fig2 <- plot_grid(
  title,
  fig2_new_design,
  ncol = 1,
  rel_heights = c(0.12, 1)
)

# ── Export ──────────────────────────────────────────────────────────────────
ggsave(
  "Figures/Fig2_new_infographic_design.png",
  final_fig2,
  width = 16,
  height = 11.5,
  dpi = 400,
  bg = "white"
)

ggsave(
  "Figures/Fig2_new_infographic_design.pdf",
  final_fig2,
  width = 16,
  height = 11.5,
  device = cairo_pdf,
  bg = "white"
)





