# ================================================================
# NGS_Kinetics_Zones_v2.R
# Reproduces NGS_Kinetics_Zones_v2.xlsx in pure R
#
# OUTPUT: NGS_Kinetics_Zones_v2.xlsx  (3 sheets)
#   Sheet 1 · Zone Labels   — numbers replaced by zone category labels
#   Sheet 2 · Delta% Change — Δ% with arrows + zone + interpretation
#   Sheet 3 · Legend        — colour reference + ABX stewardship rules
#
# Place NGS_3.xlsx in working directory and source this file.
# ================================================================

rm(list = ls())
setwd("~/Desktop/NGS_ZNA/2_NGS_Stewarship_Antibiotic")

# ── 0. PACKAGES ──────────────────────────────────────────────
pkgs <- c("tidyverse", "readxl", "lubridate", "openxlsx")
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0)
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

# ── 1. DATA ──────────────────────────────────────────────────
df_raw <- read_excel("2_NGS_Main.xlsx", sheet = 1)
df_raw <- df_raw |> rename_with(str_trim)

# ── 2. HELPERS ───────────────────────────────────────────────
parse_date_flex <- function(x) {
  x <- as.character(x)
  num <- suppressWarnings(as.numeric(x))
  dplyr::if_else(
    !is.na(num),
    as.Date(num, origin = "1899-12-30"),
    as.Date(lubridate::parse_date_time(
      x, orders = c("Ymd","dmY","dmy","d.m.Y","d/m/Y","Y-m-d"), quiet = TRUE))
  )
}

to_num <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, ",", ".")
  x <- str_replace_all(x, "^-$|^\\s*$", NA_character_)
  suppressWarnings(as.numeric(x))
}

calc_age <- function(dob, exam) {
  d1 <- suppressWarnings(parse_date_flex(dob))
  d2 <- suppressWarnings(parse_date_flex(exam))
  round(as.numeric(difftime(d2, d1, units = "days")) / 365.25, 1)
}

# ── 3. REFERENCE LIMITS ──────────────────────────────────────
REFS <- list(
  Leuk = list(uln = 10,   lln = 4,    label = "Leuk (×10⁹/l)", dec = 1),
  CRP  = list(uln = 5,    lln = 0,    label = "CRP (mg/l)",     dec = 1),
  PCT  = list(uln = 0.5,  lln = 0,    label = "PCT (ng/ml)",    dec = 2),
  Temp = list(uln = 38.3, lln = 36.0, label = "Temp (°C)",      dec = 1),
  SAPS = list(uln = 29,   lln = 0,    label = "SAPS II (pts)",  dec = 0)
)

PARAMS_ORDER <- c("CRP","PCT","Leuk","SAPS","Temp")

# Column name mapping (handles trailing spaces in raw file)
PARAM_COLS <- list(
  CRP  = c("CRP Tag 1","CRP Tag 2","CRP Tag 3","CRP Tag 4","CRP  Tag 5"),
  PCT  = c("PCT Tag 1","PCT Tag 2","PCT Tag 3","PCT Tag 4","PCT  Tag 5"),
  Leuk = c("Leuk Tag 1","Leuk Tag 2","Leuk Tag 3","Leuk Tag 4","Leuk Tag 5"),
  SAPS = c("SAPS II Tag 1","SAPS II Tag 2","SAPS II Tag 3",
           "SAPS II Tag 4","SAPS II Tag 5"),
  Temp = c("Temp Tag 1","Temp Tag 2","Temp Tag 3","Temp Tag 4","Temp Tag 5")
)
# Resolve exact column names (strip whitespace for matching)
PARAM_COLS <- lapply(PARAM_COLS, function(cols) {
  sapply(cols, function(c) {
    m <- names(df_raw)[str_trim(names(df_raw)) == str_trim(c)]
    if (length(m) > 0) m[1] else c
  }, USE.NAMES = FALSE)
})

# ── 4. ZONE CLASSIFICATION ───────────────────────────────────
# Returns (zone_label, colour_key) as a list
classify_zone <- function(val, param) {
  if (is.na(val)) return(list(label = "", colour = "EMPTY"))
  r <- REFS[[param]]
  if (param == "Leuk") {
    if (val < 2)        return(list(label = "Leukopenia severe",   colour = "RED"))
    if (val < r$lln)    return(list(label = "Leukopenia",           colour = "RED"))
    if (val <= r$uln)   return(list(label = "Normal",              colour = "GREEN"))
    if (val < 20)       return(list(label = "Leukocytosis",        colour = "ORANGE"))
    return(list(label = "Leukocytosis severe", colour = "RED"))
  }
  if (param == "CRP") {
    if (val <= r$uln) return(list(label = "Normal",    colour = "GREEN"))
    return(list(label = "Elevated", colour = "ORANGE"))
  }
  if (param == "PCT") {
    if (val <= r$uln) return(list(label = "Normal",                       colour = "GREEN"))
    if (val <= 2)     return(list(label = "Possible bacterial infection",  colour = "ORANGE"))
    return(list(label = "Sepsis likely", colour = "RED"))
  }
  if (param == "Temp") {
    if (val < 35)       return(list(label = "Hypothermia severe", colour = "RED"))
    if (val < r$lln)    return(list(label = "Hypothermia",        colour = "ORANGE"))
    if (val <= r$uln)   return(list(label = "Normal",             colour = "GREEN"))
    if (val < 40)       return(list(label = "Fever",              colour = "ORANGE"))
    return(list(label = "High fever", colour = "RED"))
  }
  if (param == "SAPS") {
    if (val < r$uln)  return(list(label = "Low risk",               colour = "GREEN"))
    if (val < 60)     return(list(label = "Elevated mortality risk", colour = "ORANGE"))
    return(list(label = "Very high mortality risk", colour = "RED"))
  }
  list(label = "", colour = "EMPTY")
}

# ── 5. DELTA% INTERPRETATION ─────────────────────────────────
interpret_delta <- function(prev, curr, param) {
  if (is.na(prev) || is.na(curr) || prev == 0)
    return(list(label = "", colour = "EMPTY"))
  r     <- REFS[[param]]
  delta <- (curr - prev) / abs(prev) * 100

  if (param == "Leuk") {
    in_n <- curr >= r$lln & curr <= r$uln
    if (in_n)                                       return(list(label = "Normalising",                       colour = "GREEN"))
    if (curr < r$lln && prev > r$uln)               return(list(label = "ALARM: Overshoot → immunoparalysis", colour = "RED"))
    if (curr < r$lln && prev < r$lln && curr < prev) return(list(label = "ALARM: Worsening leukopenia",       colour = "RED"))
    if (curr < r$lln && prev < r$lln)               return(list(label = "Still leukopenic",                  colour = "ORANGE"))
    if (curr > r$uln && delta > 50)                 return(list(label = "ALARM: Rebound",                    colour = "RED"))
    if (curr > r$uln && delta < -30)                return(list(label = "Good response",                     colour = "GREEN"))
    if (curr > r$uln && delta > 20)                 return(list(label = "Rising",                            colour = "ORANGE"))
    if (curr > r$uln && delta < -20)                return(list(label = "Decreasing",                        colour = "GREEN"))
    return(list(label = "Stable", colour = "ORANGE"))
  }

  if (param == "Temp") {
    in_n <- curr >= r$lln & curr <= r$uln
    if (in_n)                                        return(list(label = "Normalising",          colour = "GREEN"))
    if (prev >= r$uln && curr < r$lln)               return(list(label = "ALARM: Hypothermia",   colour = "RED"))
    if (prev < r$lln  && curr < prev)                return(list(label = "ALARM: Worsening",     colour = "RED"))
    if (prev < r$lln  && curr >= r$lln)              return(list(label = "Recovering",            colour = "GREEN"))
    if (delta > 1.5)                                 return(list(label = "Rising fast",           colour = "ORANGE"))
    if (delta < -1.5)                                return(list(label = "Falling fast",          colour = "GREEN"))
    return(list(label = "Stable", colour = "ORANGE"))
  }

  # Unidirectional: CRP, PCT, SAPS
  if (curr <= r$uln)   return(list(label = "Normalised",    colour = "GREEN"))
  if (delta < -50)     return(list(label = "Good response", colour = "GREEN"))
  if (delta < -20)     return(list(label = "Decreasing",    colour = "GREEN"))
  if (delta > 50)      return(list(label = "ALARM: Rebound",colour = "RED"))
  if (delta > 20)      return(list(label = "Rising",        colour = "ORANGE"))
  list(label = "Stable", colour = "ORANGE")
}

# ── 6. ABX STEWARDSHIP ───────────────────────────────────────
abx_stewardship <- function(zones_list, interps_list) {
  alarms <- c(); rebounds <- c(); improving <- c(); normalised <- c()

  for (param in c("Leuk","CRP","PCT","Temp","SAPS")) {
    zi <- zones_list[[param]]
    ii <- interps_list[[param]]
    zone_labels  <- sapply(zi, `[[`, "label")
    interp_labels <- sapply(ii, `[[`, "label")

    if (any(str_detect(interp_labels, "ALARM"),  na.rm = TRUE)) alarms    <- c(alarms,    param)
    if (any(str_detect(interp_labels, "Rebound"),na.rm = TRUE)) rebounds  <- c(rebounds,  param)
    if (any(interp_labels %in% c("Good response","Normalising","Normalised","Decreasing")))
      improving <- c(improving, param)
    last_z <- tail(zone_labels[zone_labels != ""], 1)
    if (length(last_z) > 0 && last_z %in% c("Normal","Normalised","Low risk"))
      normalised <- c(normalised, param)
  }

  n_alarm <- length(unique(alarms))
  n_reb   <- length(unique(rebounds))
  n_imp   <- length(unique(improving))
  n_nor   <- length(unique(normalised))

  leuk_i    <- sapply(interps_list[["Leuk"]], `[[`, "label")
  overshoot <- any(str_detect(leuk_i, "Overshoot|Worsening"), na.rm = TRUE)

  if (overshoot)
    return(list(rec = "ESCALATE / REVIEW URGENTLY",
                reason = "Leukopenia after leukocytosis — immunoparalysis. Review spectrum, consider immunostimulation.",
                urg = "RED"))
  if (n_reb >= 2)
    return(list(rec = "ESCALATE",
                reason = sprintf("Rebound in %d parameters (%s). Possible new focus/resistance.",
                                 n_reb, paste(unique(rebounds), collapse=", ")),
                urg = "RED"))
  if (n_alarm >= 2)
    return(list(rec = "ESCALATE",
                reason = sprintf("ALARM in %d parameters. Escalate antibiotic therapy.", n_alarm),
                urg = "RED"))
  if (n_alarm == 1 && n_imp >= 2)
    return(list(rec = "WATCH & WAIT",
                reason = sprintf("ALARM in %s but %d parameters improving. Reassess in 24h.",
                                 alarms[1], n_imp),
                urg = "ORANGE"))
  if (n_nor >= 4 && n_alarm == 0)
    return(list(rec = "DE-ESCALATE / STOP",
                reason = sprintf("%d/5 parameters normalised, no alarms. Consider stopping antibiotic.", n_nor),
                urg = "GREEN"))
  if (n_nor >= 3 && n_imp >= 2 && n_alarm == 0)
    return(list(rec = "DE-ESCALATE",
                reason = sprintf("%d normal, %d improving. Step down to narrower spectrum.", n_nor, n_imp),
                urg = "GREEN"))
  if (n_imp >= 3 && n_alarm == 0)
    return(list(rec = "CONTINUE / NARROW",
                reason = sprintf("%d parameters improving. Consider narrowing if culture available.", n_imp),
                urg = "LIGHTGREEN"))
  if (n_imp >= 1 && n_alarm == 0)
    return(list(rec = "CONTINUE",
                reason = "Partial response. Monitor daily.",
                urg = "YELLOW"))
  list(rec = "CONTINUE — REASSESS",
       reason = "No clear improvement, no ALARM. Reassess at 48-72h.",
       urg = "ORANGE")
}

# ── 7. CLINICAL METADATA ─────────────────────────────────────
META_DEF <- tribble(
  ~label,           ~src,                                                              ~width, ~align,
  "JOB_ID",         "JOB_ID",                                                           9,    "center",
  "Age",            "__CALC_AGE__",                                                      7,    "center",
  "Sex",            "Geschlecht",                                                        6,    "center",
  "ICU days",       "Tage IPS",                                                          8,    "center",
  "Hosp. days",     "Tage KH",                                                           9,    "center",
  "Died",           "Verstorben KH",                                                     6,    "center",
  "Infect. type",   "INFECTION_TYPE",                                                   15,    "left",
  "BK pos/neg",     "Blutkultur pos/neg",                                                9,    "center",
  "AB change",      "Änderung angezeigt 0=nein, 1=Esc, 2= Deesc., 3=Stop",            11,    "center",
  "Diagnose",       "Hauptdiagnose (FREITEXT)",                                         26,    "left",
  "Erreger",        "DISQVER Ergebnis",                                                 24,    "left",
  "Antiinfektiva",  "ANTIINFEKTIVA // 241015",                                          24,    "left"
)
N_META <- nrow(META_DEF)

get_meta_val <- function(row, label, src) {
  if (src == "__CALC_AGE__") {
    return(calc_age(row[["Geburtsatum"]], row[["Abnahmedatum"]]))
  }
  val <- as.character(row[[src]])
  val <- if (is.na(val)) "" else str_trim(val)
  if (label == "Died")
    val <- case_when(val == "1" ~ "Yes", val == "0" ~ "No", TRUE ~ val)
  if (label == "AB change")
    val <- case_when(val == "0" ~ "None", val == "1" ~ "Escalate",
                     val == "2" ~ "De-escalate", val == "3" ~ "Stop",
                     TRUE ~ val)
  val
}

# ── 8. OPENXLSX STYLES ───────────────────────────────────────
FILLS <- list(
  RED        = createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006"),
  ORANGE     = createStyle(fgFill = "#FFEB9C", fontColour = "#9C6500"),
  YELLOW     = createStyle(fgFill = "#FFFF99", fontColour = "#7D6608"),
  GREEN      = createStyle(fgFill = "#E2EFDA", fontColour = "#375623"),
  LIGHTGREEN = createStyle(fgFill = "#CCFFCC", fontColour = "#1E6B1E"),
  EMPTY      = createStyle(fgFill = "#F5F5F5", fontColour = "#AAAAAA"),
  HEADER     = createStyle(fgFill = "#1F4E79", fontColour = "#FFFFFF",
                           textDecoration = "bold", fontSize = 10,
                           halign = "center", valign = "center", wrapText = TRUE),
  ALARM_BOLD = createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006",
                           textDecoration = "bold", halign = "center",
                           valign = "center", wrapText = TRUE),
  DIED_YES   = createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006",
                           textDecoration = "bold", halign = "center"),
  DIED_NO    = createStyle(fgFill = "#E2EFDA", fontColour = "#375623",
                           halign = "center")
)

zone_style <- function(colour, is_alarm = FALSE) {
  if (is_alarm) return(FILLS[["ALARM_BOLD"]])
  base <- FILLS[[colour]]
  if (is.null(base)) return(FILLS[["EMPTY"]])
  createStyle(fgFill   = base$fgFill,
              fontColour = base$fontColour,
              halign   = "center", valign = "center", wrapText = TRUE)
}

PARAM_HDR_FILLS <- c(
  CRP  = "#2E75B6", PCT  = "#7030A0",
  Leuk = "#C55A11", SAPS = "#375623", Temp = "#843C0C"
)

# ── 9. BUILD WORKBOOK ─────────────────────────────────────────
wb <- createWorkbook()

# Column layout (same for both data sheets)
col_start <- list()
cc <- N_META + 1L
for (p in PARAMS_ORDER) { col_start[[p]] <- cc; cc <- cc + 5L }
rec_col    <- cc
reason_col <- cc + 1L

# Helper: write shared headers (rows 1+2) to a worksheet
write_sheet_headers <- function(ws, tag_label_fn) {

  # Meta: merge rows 1+2
  for (i in seq_len(N_META)) {
    addStyle(wb, ws, FILLS[["HEADER"]], rows = 1:2, cols = i, gridExpand = TRUE)
    writeData(wb, ws, META_DEF$label[i], startRow = 1, startCol = i)
    mergeCells(wb, ws, rows = 1:2, cols = i)
    setColWidths(wb, ws, cols = i, widths = META_DEF$width[i])
  }

  # Parameter group headers (row 1) + tag sub-headers (row 2)
  for (p in PARAMS_ORDER) {
    cs  <- col_start[[p]]
    hc  <- PARAM_HDR_FILLS[[p]]
    grp_style <- createStyle(fgFill = hc, fontColour = "#FFFFFF",
                             textDecoration = "bold", fontSize = 10,
                             halign = "center", valign = "center")
    sub_style <- createStyle(fgFill = hc, fontColour = "#FFFFCC",
                             textDecoration = "bold", fontSize = 9,
                             halign = "center", valign = "center")
    writeData(wb, ws, REFS[[p]]$label, startRow = 1, startCol = cs)
    mergeCells(wb, ws, rows = 1, cols = cs:(cs+4))
    addStyle(wb, ws, grp_style, rows = 1, cols = cs:(cs+4), gridExpand = TRUE)

    for (t in 0:4) {
      lbl <- tag_label_fn(t)
      writeData(wb, ws, lbl, startRow = 2, startCol = cs + t)
      addStyle(wb, ws, sub_style, rows = 2, cols = cs + t)
      setColWidths(wb, ws, cols = cs + t,
                   widths = if (ws == "1 · Zone Labels") 22 else 20)
    }
  }

  # ABX headers
  for (col in c(rec_col, reason_col)) {
    lbl <- if (col == rec_col) "ABX Stewardship" else "Clinical Reasoning"
    writeData(wb, ws, lbl, startRow = 1, startCol = col)
    mergeCells(wb, ws, rows = 1:2, cols = col)
    addStyle(wb, ws, FILLS[["HEADER"]], rows = 1:2, cols = col, gridExpand = TRUE)
    setColWidths(wb, ws, cols = col,
                 widths = if (col == rec_col) 22 else 52)
  }

  setRowHeights(wb, ws, rows = 1, heights = 22)
  setRowHeights(wb, ws, rows = 2, heights = 28)
  freezePane(wb, ws, firstRow = FALSE,
             firstActiveRow = 3, firstActiveCol = N_META + 1)
}

# ── SHEET 1: ZONE LABELS ─────────────────────────────────────
addWorksheet(wb, "1 · Zone Labels")
write_sheet_headers("1 · Zone Labels", function(t) paste0("T", t + 1))

for (row_idx in seq_len(nrow(df_raw))) {
  row <- df_raw[row_idx, ]
  r   <- row_idx + 2L   # Excel row (row 1+2 = headers)

  setRowHeights(wb, "1 · Zone Labels", rows = r, heights = 40)

  # Meta columns
  for (i in seq_len(N_META)) {
    val <- get_meta_val(row, META_DEF$label[i], META_DEF$src[i])
    writeData(wb, "1 · Zone Labels", val, startRow = r, startCol = i)
    sty <- createStyle(
      halign = META_DEF$align[i], valign = "center",
      fontName = "Arial", fontSize = 9
    )
    if (META_DEF$label[i] == "Died" && val == "Yes")
      sty <- FILLS[["DIED_YES"]]
    else if (META_DEF$label[i] == "Died" && val == "No")
      sty <- FILLS[["DIED_NO"]]
    addStyle(wb, "1 · Zone Labels", sty, rows = r, cols = i)
  }

  # Kinetics cells — zone label only
  zones_list  <- list()
  interps_list <- list()

  for (p in PARAMS_ORDER) {
    vals   <- sapply(PARAM_COLS[[p]], function(c) to_num(row[[c]]))
    zones  <- lapply(vals, classify_zone, param = p)
    interps <- lapply(seq_along(vals), function(t) {
      prev <- if (t > 1) vals[t-1] else NA_real_
      interpret_delta(prev, vals[t], p)
    })
    zones_list[[p]]   <- zones
    interps_list[[p]] <- interps

    cs <- col_start[[p]]
    for (t in seq_along(vals)) {
      v  <- vals[t]
      zl <- zones[[t]]$label
      zc <- zones[[t]]$colour
      il <- interps[[t]]$label
      ic <- interps[[t]]$colour

      if (is.na(v)) {
        writeData(wb, "1 · Zone Labels", "—", startRow = r, startCol = cs + t - 1)
        addStyle(wb, "1 · Zone Labels", FILLS[["EMPTY"]],
                 rows = r, cols = cs + t - 1)
      } else {
        # Sheet 1: zone label + interpretation (no raw number)
        cell_text <- if (nchar(il) > 0) paste0(zl, "\n", il) else zl
        writeData(wb, "1 · Zone Labels", cell_text,
                  startRow = r, startCol = cs + t - 1)
        is_alarm <- str_detect(il, "ALARM")
        final_col <- if (is_alarm) "RED" else zc
        sty <- createStyle(
          fgFill     = switch(final_col,
                              RED="#FFC7CE", ORANGE="#FFEB9C",
                              GREEN="#E2EFDA", "#F5F5F5"),
          fontColour = switch(final_col,
                              RED="#9C0006", ORANGE="#9C6500",
                              GREEN="#375623", "#AAAAAA"),
          textDecoration = if (is_alarm) "bold" else NULL,
          halign = "center", valign = "center",
          wrapText = TRUE, fontName = "Arial", fontSize = 9
        )
        addStyle(wb, "1 · Zone Labels", sty, rows = r, cols = cs + t - 1)
      }
    }
  }

  # ABX Stewardship
  abx <- abx_stewardship(zones_list, interps_list)
  urg_fill <- switch(abx$urg,
    RED="#FFC7CE", ORANGE="#FFEB9C", YELLOW="#FFFF99",
    GREEN="#E2EFDA", LIGHTGREEN="#CCFFCC", "#FFFF99")
  urg_font <- switch(abx$urg,
    RED="#9C0006", ORANGE="#9C6500", YELLOW="#7D6608",
    GREEN="#375623", LIGHTGREEN="#1E6B1E", "#7D6608")

  writeData(wb, "1 · Zone Labels", abx$rec,    startRow = r, startCol = rec_col)
  writeData(wb, "1 · Zone Labels", abx$reason, startRow = r, startCol = reason_col)

  for (col in c(rec_col, reason_col)) {
    addStyle(wb, "1 · Zone Labels",
             createStyle(fgFill = urg_fill, fontColour = urg_font,
                         textDecoration = if (col == rec_col) "bold" else NULL,
                         fontSize = if (col == rec_col) 10 else 9,
                         fontName = "Arial",
                         halign = if (col == rec_col) "center" else "left",
                         valign = "center", wrapText = TRUE),
             rows = r, cols = col)
  }
}

# ── SHEET 2: DELTA% CHANGE ───────────────────────────────────
addWorksheet(wb, "2 · Delta% Change")
write_sheet_headers("2 · Delta% Change",
                    function(t) if (t == 0) "T1" else paste0("T", t, "\u2192T", t+1))

for (row_idx in seq_len(nrow(df_raw))) {
  row <- df_raw[row_idx, ]
  r   <- row_idx + 2L

  setRowHeights(wb, "2 · Delta% Change", rows = r, heights = 44)

  # Meta
  for (i in seq_len(N_META)) {
    val <- get_meta_val(row, META_DEF$label[i], META_DEF$src[i])
    writeData(wb, "2 · Delta% Change", val, startRow = r, startCol = i)
    sty <- createStyle(halign = META_DEF$align[i], valign = "center",
                       fontName = "Arial", fontSize = 9)
    if (META_DEF$label[i] == "Died" && val == "Yes") sty <- FILLS[["DIED_YES"]]
    else if (META_DEF$label[i] == "Died" && val == "No") sty <- FILLS[["DIED_NO"]]
    addStyle(wb, "2 · Delta% Change", sty, rows = r, cols = i)
  }

  zones_list  <- list()
  interps_list <- list()

  for (p in PARAMS_ORDER) {
    vals    <- sapply(PARAM_COLS[[p]], function(c) to_num(row[[c]]))
    zones   <- lapply(vals, classify_zone, param = p)
    interps <- lapply(seq_along(vals), function(t) {
      prev <- if (t > 1) vals[t-1] else NA_real_
      interpret_delta(prev, vals[t], p)
    })
    zones_list[[p]]   <- zones
    interps_list[[p]] <- interps

    cs  <- col_start[[p]]
    dec <- REFS[[p]]$dec
    fmt <- paste0("%.", dec, "f")

    for (t in seq_along(vals)) {
      v   <- vals[t]
      zl  <- zones[[t]]$label
      zc  <- zones[[t]]$colour
      il  <- interps[[t]]$label

      col_idx <- cs + t - 1

      if (is.na(v)) {
        writeData(wb, "2 · Delta% Change", "—", startRow = r, startCol = col_idx)
        addStyle(wb, "2 · Delta% Change", FILLS[["EMPTY"]], rows = r, cols = col_idx)
        next
      }

      val_str <- sprintf(fmt, v)

      if (t == 1) {
        # T1: baseline — raw value + zone label, no delta
        cell_text <- paste0(val_str, "\n[", zl, "]")
        dcol <- zc
        is_alarm <- FALSE
      } else {
        prev <- vals[t-1]
        if (is.na(prev) || prev == 0) {
          cell_text <- paste0(val_str, "\n[", zl, "]\n\u0394% = n/a")
          dcol <- zc
          is_alarm <- FALSE
        } else {
          delta    <- (v - prev) / abs(prev) * 100
          sign_str <- if (delta >= 0) "+" else ""
          arrow    <- if (delta > 0) "\u25b2" else if (delta < 0) "\u25bc" else "\u25ba"
          is_alarm <- str_detect(il, "ALARM")
          dcol     <- if (is_alarm)    "RED"
                      else if (delta > 20) "ORANGE"
                      else if (delta < -20) "GREEN"
                      else "NEUTRAL"
          cell_text <- paste0(arrow, " ", sign_str, sprintf("%.1f", delta), "%",
                              "\n", val_str, " [", zl, "]",
                              "\n", il)
        }
      }

      writeData(wb, "2 · Delta% Change", cell_text, startRow = r, startCol = col_idx)

      fill_hex <- switch(dcol,
        RED="#FFC7CE", ORANGE="#FFEB9C", GREEN="#E2EFDA",
        NEUTRAL="#F5F5F5", "#F5F5F5")
      font_hex <- switch(dcol,
        RED="#9C0006", ORANGE="#9C6500", GREEN="#375623",
        NEUTRAL="#555555", "#555555")

      addStyle(wb, "2 · Delta% Change",
               createStyle(fgFill = fill_hex, fontColour = font_hex,
                           textDecoration = if (is_alarm) "bold" else NULL,
                           fontName = "Arial", fontSize = 9,
                           halign = "center", valign = "center", wrapText = TRUE),
               rows = r, cols = col_idx)
    }
  }

  # ABX
  abx <- abx_stewardship(zones_list, interps_list)
  urg_fill <- switch(abx$urg,
    RED="#FFC7CE", ORANGE="#FFEB9C", YELLOW="#FFFF99",
    GREEN="#E2EFDA", LIGHTGREEN="#CCFFCC", "#FFFF99")
  urg_font <- switch(abx$urg,
    RED="#9C0006", ORANGE="#9C6500", YELLOW="#7D6608",
    GREEN="#375623", LIGHTGREEN="#1E6B1E", "#7D6608")

  writeData(wb, "2 · Delta% Change", abx$rec,    startRow = r, startCol = rec_col)
  writeData(wb, "2 · Delta% Change", abx$reason, startRow = r, startCol = reason_col)

  for (col in c(rec_col, reason_col)) {
    addStyle(wb, "2 · Delta% Change",
             createStyle(fgFill = urg_fill, fontColour = urg_font,
                         textDecoration = if (col == rec_col) "bold" else NULL,
                         fontSize = if (col == rec_col) 10 else 9,
                         fontName = "Arial",
                         halign = if (col == rec_col) "center" else "left",
                         valign = "center", wrapText = TRUE),
             rows = r, cols = col)
  }
}

# ── SHEET 3: LEGEND ──────────────────────────────────────────
addWorksheet(wb, "3 · Legend")
setColWidths(wb, "3 · Legend", cols = 1:3, widths = c(28, 32, 48))

legend_rows <- list(
  list("SHEET 1 — ZONE LABELS",  "Cell content",         "Zone label + interpretation (no raw number)", "HEADER"),
  list("",  "GREEN cell",         "Normal / good response / normalising",         "GREEN"),
  list("",  "ORANGE cell",        "Warning — elevated, leukocytosis, fever",      "ORANGE"),
  list("",  "RED cell",           "ALARM — severe, rebound, immunoparalysis",     "RED"),
  list("",  "T1 cell",            "Zone label only (no delta at baseline)",        ""),
  list("",  "T2-T5 cell",         "Zone label + interpretation",                  ""),
  list("",  "—",                  "No measurement available",                     ""),
  list("",  "",                   "",                                              ""),
  list("SHEET 2 — DELTA% CHANGE","Cell content",         "Δ% + raw value + zone + interpretation",      "HEADER"),
  list("",  "▲ +Δ% RED",          "Rising > +20% or ALARM",                       "RED"),
  list("",  "▲ +Δ% ORANGE",       "Rising warning",                               "ORANGE"),
  list("",  "▼ -Δ% GREEN",        "Falling > -20% (good response)",               "GREEN"),
  list("",  "► ±Δ% GREY",         "Stable change",                                ""),
  list("",  "T1 baseline",        "Raw value + [zone] — no delta",                ""),
  list("",  "",                   "",                                              ""),
  list("ZONE REFERENCE",          "Parameter",            "Range → Zone",          "HEADER"),
  list("",  "Leuk < 2",           "Leukopenia severe",                            "RED"),
  list("",  "Leuk 2-4",           "Leukopenia",                                   "RED"),
  list("",  "Leuk 4-10",          "Normal",                                       "GREEN"),
  list("",  "Leuk 10-20",         "Leukocytosis",                                 "ORANGE"),
  list("",  "Leuk > 20",          "Leukocytosis severe",                          "RED"),
  list("",  "CRP < 5 mg/l",       "Normal",                                       "GREEN"),
  list("",  "CRP > 5 mg/l",       "Elevated",                                     "ORANGE"),
  list("",  "PCT < 0.5",          "Normal",                                       "GREEN"),
  list("",  "PCT 0.5-2",          "Possible bacterial infection",                 "ORANGE"),
  list("",  "PCT > 2",            "Sepsis likely",                                "RED"),
  list("",  "Temp < 35°C",        "Hypothermia severe",                           "RED"),
  list("",  "Temp 35-36°C",       "Hypothermia",                                  "ORANGE"),
  list("",  "Temp 36-38.3°C",     "Normal",                                       "GREEN"),
  list("",  "Temp 38.3-40°C",     "Fever",                                        "ORANGE"),
  list("",  "Temp > 40°C",        "High fever",                                   "RED"),
  list("",  "SAPS < 29",          "Low risk",                                     "GREEN"),
  list("",  "SAPS 29-59",         "Elevated mortality risk",                      "ORANGE"),
  list("",  "SAPS >= 60",         "Very high mortality risk",                     "RED"),
  list("",  "",                   "",                                              ""),
  list("ABX STEWARDSHIP",         "Recommendation",       "Trigger",              "HEADER"),
  list("",  "ESCALATE / REVIEW URGENTLY", "Leukopenia after leukocytosis (immunoparalysis)", "RED"),
  list("",  "ESCALATE",           "Rebound >=2 params OR ALARM >=2 params",       "RED"),
  list("",  "WATCH & WAIT",       "ALARM 1 param, >=2 improving",                 "ORANGE"),
  list("",  "CONTINUE — REASSESS","No improvement, no ALARM at 48-72h",           "ORANGE"),
  list("",  "CONTINUE",           "Partial response",                             "YELLOW"),
  list("",  "CONTINUE / NARROW",  ">=3 improving, no alarms",                     "LIGHTGREEN"),
  list("",  "DE-ESCALATE",        ">=3 normal + >=2 improving",                   "GREEN"),
  list("",  "DE-ESCALATE / STOP", ">=4 normalised, no alarms",                    "GREEN")
)

for (ri in seq_along(legend_rows)) {
  entry  <- legend_rows[[ri]]
  fk     <- entry[[4]]
  setRowHeights(wb, "3 · Legend", rows = ri, heights = 16)

  for (ci in 1:3) {
    writeData(wb, "3 · Legend", entry[[ci]], startRow = ri, startCol = ci)
    if (fk == "HEADER") {
      addStyle(wb, "3 · Legend", FILLS[["HEADER"]], rows = ri, cols = ci)
    } else if (fk != "" && ci %in% c(2,3)) {
      fill_hex <- switch(fk,
        RED="#FFC7CE", ORANGE="#FFEB9C", YELLOW="#FFFF99",
        GREEN="#E2EFDA", LIGHTGREEN="#CCFFCC", "#F5F5F5")
      font_hex <- switch(fk,
        RED="#9C0006", ORANGE="#9C6500", YELLOW="#7D6608",
        GREEN="#375623", LIGHTGREEN="#1E6B1E", "#333333")
      addStyle(wb, "3 · Legend",
               createStyle(fgFill = fill_hex, fontColour = font_hex,
                           textDecoration = if (ci == 2) "bold" else NULL,
                           fontName = "Arial", fontSize = 10,
                           halign = "left", valign = "center"),
               rows = ri, cols = ci)
    } else {
      addStyle(wb, "3 · Legend",
               createStyle(textDecoration = if (ci == 1 && entry[[1]] != "") "bold" else NULL,
                           fontName = "Arial", fontSize = 10,
                           halign = "left", valign = "center"),
               rows = ri, cols = ci)
    }
  }
}

# ── 10. SAVE ─────────────────────────────────────────────────
saveWorkbook(wb, "NGS_Kinetics_Zones_v2.xlsx", overwrite = TRUE)

cat(sprintf("\nSaved: NGS_Kinetics_Zones_v2.xlsx\n"))
cat(sprintf("Patients : %d\n", nrow(df_raw)))
cat(sprintf("Columns  : %d meta + %d kinetics + 2 ABX = %d total\n",
            N_META, length(PARAMS_ORDER) * 5, N_META + length(PARAMS_ORDER)*5 + 2))
cat("\nSheets:\n")
cat("  1 · Zone Labels   — numbers replaced by zone category labels\n")
cat("  2 · Delta% Change — delta% with arrows, raw value, zone, interpretation\n")
cat("  3 · Legend        — colour reference + ABX stewardship rules\n")
