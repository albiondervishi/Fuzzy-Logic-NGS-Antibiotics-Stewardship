# ================================================================
# NGS_Fuzzy_ABX_Stewardship_v7_Prina.R
#
# KEY UPGRADES FROM v6:
#   Marker-specific kinetic thresholds replacing uniform 20%
#
#   RATIONALE:
#   A 20% delta is biologically meaningless for Temperature (°C)
#   and too lenient for PCT (literature standard = 50% fall).
#   Each marker has different kinetics, units, and clinical meaning.
#
#   THRESHOLD STRATEGY (hybrid — most defensible for publication):
#
#   PCT   → Literature: 50% fall = response (PRORATA trial, 2010)
#                        80% fall = strong response
#                        rising   = alarm (any rise >20%)
#
#   CRP   → Literature: 25% fall = mild response (Póvoa 2005)
#                        50% fall = good response
#                        Data-derived for alarm threshold
#
#   Leuk  → Literature: bidirectional (both rise AND fall alarm)
#                        >50% rise  = alarm (rebound)
#                        >50% fall from high = possible immunoparalysis
#                        Data-derived optimal cut
#
#   Temp  → ABSOLUTE °C change (not %) — % change is meaningless
#                        fall >1.0°C  = mild response
#                        fall >1.5°C  = strong response
#                        rise >1.0°C  = warning
#                        rise >1.5°C  = alarm
#
#   SAPS  → No literature standard for kinetics
#                        Data-derived from your cohort vs ICU mortality
#
#   All thresholds reported in a methods table (Fig_Thresholds)
#   with source (literature vs data-derived) and 95% CI where applicable
#
# NEW IN v7:
#   [1] Diagnostic group variable (mNGS/BC/Abx combination)
#   [2] AUC of full model (danger_score) vs ICU mortality
#   [3] Fig14 — Diagnostic group alluvial + ICU mortality
#   [4] Fig15 — Decision by diagnostic group (stacked bar, matches publication table)
#   [5] Fig16 — Danger score vs ICU mortality by diagnostic group
#   [6] Statistics Excel sheet for publication (Table 1-8)
#   [7] Hospital mortality analysis (more events, better power)
#   [8] Stop vs Rest binary safety comparison (Fisher exact + danger score)
#   [9] Fig17 — Stop vs Rest safety figure
#
# OUTPUT
#   ./Figures/   — Fig1–Fig16 (PNG + PDF)
#   ./Output/    — Results Excel + Publication_Statistics Excel
# ================================================================

rm(list = ls())
#setwd("~/Desktop/NGS_ZNA")
setwd("~/Desktop/NGS_ZNA/2_NGS_Stewarship_Antibiotic")
# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
pkgs <- c(
  "tidyverse", "readxl", "janitor", "nnet", "pROC",
  "caret", "writexl", "scales", "patchwork",
  "RColorBrewer", "broom", "ggalluvial", "gridExtra"
)
new_pkgs <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(new_pkgs) > 0)
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

dir.create("Figures", showWarnings = FALSE)
dir.create("Output",  showWarnings = FALSE)

# ── 1. MARKER-SPECIFIC THRESHOLD DEFINITIONS ──────────────────────────────────
# Each marker has: strong_response, mild_response, mild_alarm, strong_alarm
# Temperature uses ABSOLUTE °C; all others use delta%
# Source: L = Literature reference; D = Data-derived from this cohort
#
# These are used in delta_score_mk() below.
# After running the data-derived optimisation (Section 7),
# D thresholds are updated automatically from your data.

MARKER_CUTS <- list(
  
  pct = list(
    unit          = "ng/mL",
    metric        = "delta_pct",       # delta% used
    strong_resp   = -50,               # L: PRORATA trial (de Jong 2016)
    mild_resp     = -25,               # L: Schuetz 2018 meta-analysis
    mild_alarm    =  20,               # L: any rise concerning
    strong_alarm  =  50,               # L: rebound definition
    source_resp   = "Literature (PRORATA/Schuetz)",
    source_alarm  = "Literature",
    bidirectional = FALSE              # only high PCT is dangerous
  ),
  
  crp = list(
    unit          = "mg/L",
    metric        = "delta_pct",
    strong_resp   = -50,               # L: Póvoa 2005 (ICM)
    mild_resp     = -25,               # L: Povoa 2005
    mild_alarm    =  25,               # D: data-derived (updated below)
    strong_alarm  =  50,               # D: data-derived (updated below)
    source_resp   = "Literature (Póvoa 2005)",
    source_alarm  = "Data-derived",
    bidirectional = FALSE
  ),
  
  leuk = list(
    unit          = "×10⁹/L",
    metric        = "delta_pct",
    strong_resp   = -30,               # D: data-derived
    mild_resp     = -15,               # D: data-derived
    mild_alarm    =  30,               # L: leukocytosis rebound concern
    strong_alarm  =  60,               # L: severe rebound / >20×10⁹
    source_resp   = "Data-derived",
    source_alarm  = "Literature (immunoparalysis threshold)",
    bidirectional = TRUE               # both rise AND fall can be alarm
  ),
  
  temp = list(
    unit          = "°C",
    metric        = "delta_abs",       # ABSOLUTE change used, NOT %
    strong_resp   = -1.5,              # L: clinical standard (defervescence)
    mild_resp     = -1.0,              # L: partial defervescence
    mild_alarm    =  1.0,              # L: fever spike
    strong_alarm  =  1.5,             # L: high fever alarm
    source_resp   = "Literature (clinical defervescence standard)",
    source_alarm  = "Literature",
    bidirectional = FALSE
  ),
  
  saps = list(
    unit          = "points",
    metric        = "delta_pct",
    strong_resp   = -20,               # D: data-derived
    mild_resp     = -10,               # D: data-derived
    mild_alarm    =  15,               # D: data-derived
    strong_alarm  =  30,               # D: data-derived
    source_resp   = "Data-derived",
    source_alarm  = "Data-derived",
    bidirectional = FALSE
  )
)

# Convenience: PARAM block for decision rule cuts (unchanged from v5)
PARAM <- list(
  escalate_cut       = 8.0,
  deesc_response_cut = 5.0,
  stop_response_cut  = 6.0,
  stop_danger_cut    = 3.0,
  w_pct  = 1.60,  w_saps = 1.50,  w_temp = 1.10,
  w_leuk = 1.00,  w_crp  = 0.60,
  w_bk   = 0.00,  w_mngs = 1.00,  w_age  = 0.50,
  w_complexity = 1.00
)

# ── 2. READ WORKBOOK ──────────────────────────────────────────────────────────
INPUT_FILE <- "~/Desktop/NGS_ZNA/NGS_Stewarship_Antibiotics/NGS_Kinetics_Zones_v2.xlsx"
if (!file.exists(INPUT_FILE)) stop("Cannot find ", INPUT_FILE)

sheets <- readxl::excel_sheets(INPUT_FILE)
cat("Sheets found:", paste(sheets, collapse = ", "), "\n")

read_sheet_clean <- function(sn)
  readxl::read_excel(INPUT_FILE, sheet = sn, guess_max = 10000) |>
  janitor::clean_names()

zone_df_raw  <- read_sheet_clean(sheets[1])
delta_df_raw <- if (length(sheets) >= 2) read_sheet_clean(sheets[2]) else zone_df_raw

# Deduplicate on JOB_ID — keep first occurrence
id_raw <- names(zone_df_raw)[stringr::str_detect(
  names(zone_df_raw), stringr::regex("^job_id$|patient|case|^id$", ignore_case=TRUE))][1]

cat(sprintf("Raw rows: %d\n", nrow(zone_df_raw)))
zone_df_raw  <- zone_df_raw  |> dplyr::distinct(.data[[id_raw]], .keep_all=TRUE)
delta_df_raw <- delta_df_raw |> dplyr::distinct(.data[[id_raw]], .keep_all=TRUE)
cat(sprintf("After deduplication: %d unique patients\n", nrow(zone_df_raw)))

# ── 3. COLUMN FINDERS ─────────────────────────────────────────────────────────
find_col <- function(df, patterns, required=FALSE, label="column") {
  hit <- names(df)[stringr::str_detect(
    names(df), stringr::regex(paste(patterns, collapse="|"), ignore_case=TRUE))]
  if (!length(hit)) { if (required) stop("Cannot find: ", label); return(NA_character_) }
  hit[1]
}

id_col        <- find_col(zone_df_raw, c("^job_id$","patient","case","^id$"),   required=TRUE)
ab_col        <- find_col(zone_df_raw, c("^ab_change$","antibiotic.*change","abx.*change","steward","decision"), required=TRUE)
age_col       <- find_col(zone_df_raw, c("^age$","alter","years"))
infect_col    <- find_col(zone_df_raw, c("infect.*type","infection.*type","micro.*type"))
bk_col        <- find_col(zone_df_raw, c("^bk","blood.*culture","blut","bc_pos"))
mngs_col      <- find_col(zone_df_raw, c("^erreger$","erreger","disqver","mngs","ngs"))
abx_col       <- find_col(zone_df_raw, c("^abx_at_sampling_raw$","abx.*sampling.*raw","abx_bin","antibiotic.*on"))
died_col      <- find_col(zone_df_raw, c("died","verstorben","death","mortalit"))
icu_days_col  <- find_col(zone_df_raw, c("icu.*day","tage.*ips","ips.*tag"))
hosp_days_col <- find_col(zone_df_raw, c("hosp.*day","tage.*kh","kh.*tag"))
icu_mort_col  <- find_col(zone_df_raw, c("icu.*mort","icu_mort"))

# ── 4. NORMALISATION HELPERS ──────────────────────────────────────────────────
clean_decision <- function(x) {
  y <- stringr::str_squish(as.character(x))
  dplyr::case_when(
    stringr::str_detect(y, stringr::regex("de.?escal",  ignore_case=TRUE)) ~ "De-escalate",
    stringr::str_detect(y, stringr::regex("escal",      ignore_case=TRUE)) ~ "Escalate",
    stringr::str_detect(y, stringr::regex("^stop$|discontinue|cessation|^off$", ignore_case=TRUE)) ~ "Stop",
    stringr::str_detect(y, stringr::regex("no.?change|none|continue|watch|wait|unchanged", ignore_case=TRUE)) ~ "No change",
    is.na(y) | y %in% c("","—","-","NA","na") ~ NA_character_,
    TRUE ~ y
  )
}

bin01 <- function(x) {
  y <- tolower(stringr::str_squish(as.character(x)))
  dplyr::case_when(
    y %in% c("1","pos","positive","yes","ja","true")       ~ 1L,
    stringr::str_detect(y, "pos|detected|growth|wachstum") ~ 1L,
    y %in% c("0","neg","negative","no","nein","false")     ~ 0L,
    stringr::str_detect(y, "neg|not detected|kein")        ~ 0L,
    is.na(x) ~ 0L,
    TRUE ~ suppressWarnings(as.integer(as.numeric(y)))
  ) |> tidyr::replace_na(0L)
}

extract_number <- function(x)
  readr::parse_number(as.character(x),
                      locale = readr::locale(decimal_mark="."),
                      na = c("","NA","—","-"))

last_nonmissing  <- function(x) { x <- x[!is.na(x) & x != ""]; if (!length(x)) NA else x[length(x)] }
first_nonmissing <- function(x) { x <- x[!is.na(x) & x != ""]; if (!length(x)) NA else x[1] }

# ── 5. ZONE SEVERITY ─────────────────────────────────────────────────────────
zone_severity <- function(zone) {
  z <- tolower(as.character(zone))
  dplyr::case_when(
    is.na(zone) | z == ""                                                      ~ 0,
    stringr::str_detect(z, "normal|low risk|normalised")                       ~ 0,
    stringr::str_detect(z, "leukopenia$|possible bacterial|elevated$|fever$|hypothermia$|elevated mortality") ~ 1.5,
    stringr::str_detect(z, "leukopenia severe")                                ~ 2.5,
    stringr::str_detect(z, "sepsis likely|high fever|very high mortality|leukocytosis severe|severe") ~ 3.0,
    TRUE ~ 1.0
  )
}

extract_zone_text <- function(x) {
  s <- stringr::str_squish(as.character(x))
  dplyr::case_when(
    stringr::str_detect(s, stringr::regex("leukopenia severe",   ignore_case=TRUE)) ~ "Leukopenia severe",
    stringr::str_detect(s, stringr::regex("leukopenia",          ignore_case=TRUE)) ~ "Leukopenia",
    stringr::str_detect(s, stringr::regex("leukocytosis severe", ignore_case=TRUE)) ~ "Leukocytosis severe",
    stringr::str_detect(s, stringr::regex("leukocytosis",        ignore_case=TRUE)) ~ "Leukocytosis",
    stringr::str_detect(s, stringr::regex("possible bacterial",  ignore_case=TRUE)) ~ "Possible bacterial infection",
    stringr::str_detect(s, stringr::regex("sepsis likely",       ignore_case=TRUE)) ~ "Sepsis likely",
    stringr::str_detect(s, stringr::regex("hypothermia severe",  ignore_case=TRUE)) ~ "Hypothermia severe",
    stringr::str_detect(s, stringr::regex("hypothermia",         ignore_case=TRUE)) ~ "Hypothermia",
    stringr::str_detect(s, stringr::regex("high fever",          ignore_case=TRUE)) ~ "High fever",
    stringr::str_detect(s, stringr::regex("fever",               ignore_case=TRUE)) ~ "Fever",
    stringr::str_detect(s, stringr::regex("very high mortality", ignore_case=TRUE)) ~ "Very high mortality risk",
    stringr::str_detect(s, stringr::regex("elevated mortality",  ignore_case=TRUE)) ~ "Elevated mortality risk",
    stringr::str_detect(s, stringr::regex("low risk",            ignore_case=TRUE)) ~ "Low risk",
    stringr::str_detect(s, stringr::regex("normalised|normal",   ignore_case=TRUE)) ~ "Normal",
    stringr::str_detect(s, stringr::regex("elevated",            ignore_case=TRUE)) ~ "Elevated",
    TRUE ~ NA_character_
  )
}

# ── 6. MARKER-SPECIFIC DELTA SCORE FUNCTION ───────────────────────────────────
# Replaces the uniform delta_score() from v5
# Temperature uses absolute °C change; all others use delta%
# Returns -2 (strong improvement) to +2 (strong alarm)

delta_score_mk <- function(delta_val, mk) {
  cuts <- MARKER_CUTS[[mk]]
  if (is.null(cuts)) return(rep(0, length(delta_val)))
  
  d <- as.numeric(delta_val)
  
  if (cuts$metric == "delta_abs") {
    # Temperature: absolute °C
    dplyr::case_when(
      is.na(d)                      ~  0,
      d <= cuts$strong_resp         ~ -2,
      d <  cuts$mild_resp           ~ -1,
      d >= cuts$strong_alarm        ~  2,
      d >  cuts$mild_alarm          ~  1,
      TRUE                          ~  0
    )
  } else {
    # All others: delta%
    if (cuts$bidirectional) {
      # Leukocytes: large fall also alarming (immunoparalysis)
      dplyr::case_when(
        is.na(d)                      ~  0,
        d <= cuts$strong_resp         ~ -2,   # large fall = good (from high)
        d <  cuts$mild_resp           ~ -1,
        d >= cuts$strong_alarm        ~  2,   # large rise = alarm
        d >  cuts$mild_alarm          ~  1,
        d <  -60                      ~  2,   # extreme fall = immunoparalysis alarm
        TRUE                          ~  0
      )
    } else {
      dplyr::case_when(
        is.na(d)                      ~  0,
        d <= cuts$strong_resp         ~ -2,
        d <  cuts$mild_resp           ~ -1,
        d >= cuts$strong_alarm        ~  2,
        d >  cuts$mild_alarm          ~  1,
        TRUE                          ~  0
      )
    }
  }
}

# ── 7. FEATURE CONSTRUCTION ───────────────────────────────────────────────────
markers       <- c("leuk","crp","pct","temp","saps")
marker_labels <- c(leuk="Leukocytes", crp="CRP", pct="PCT",
                   temp="Temperature", saps="SAPS II")

get_marker_cols <- function(df, marker)
  names(df)[stringr::str_detect(names(df),
                                stringr::regex(marker, ignore_case=TRUE))]

# ── 7a. Base dataframe ────────────────────────────────────────────────────────
base_df <- zone_df_raw |>
  dplyr::transmute(
    JOB_ID      = .data[[id_col]],
    AB_change   = clean_decision(.data[[ab_col]]),
    Age         = if (!is.na(age_col))    extract_number(.data[[age_col]])  else NA_real_,
    Infect_type = if (!is.na(infect_col)) as.character(.data[[infect_col]]) else NA_character_,
    BK_bin      = if (!is.na(bk_col))    bin01(.data[[bk_col]])             else 0L,
    mNGS_bin    = if (!is.na(mngs_col))  bin01(.data[[mngs_col]])           else 0L,
    Abx_bin     = if (!is.na(abx_col))   bin01(.data[[abx_col]])            else 1L,
    
    ICU_mortality = {
      if (!is.na(icu_mort_col)) {
        bin01(.data[[icu_mort_col]])
      } else if (!is.na(died_col) & !is.na(icu_days_col) & !is.na(hosp_days_col)) {
        icu_d  <- suppressWarnings(as.numeric(as.character(.data[[icu_days_col]])))
        hosp_d <- suppressWarnings(as.numeric(as.character(.data[[hosp_days_col]])))
        died_v <- bin01(.data[[died_col]])
        as.integer(!is.na(icu_d) & !is.na(hosp_d) & died_v==1L & icu_d==hosp_d)
      } else { NA_integer_ }
    },
    
    ICU_days  = if (!is.na(icu_days_col))  suppressWarnings(as.numeric(as.character(.data[[icu_days_col]])))  else NA_real_,
    Hosp_days = if (!is.na(hosp_days_col)) suppressWarnings(as.numeric(as.character(.data[[hosp_days_col]]))) else NA_real_,
    Died      = if (!is.na(died_col))      bin01(.data[[died_col]])                                           else NA_integer_
  )


# Find the "durchgeführt" column
durchgefuehrt_col <- find_col(zone_df_raw, 
                              c("durchgef", "executed", "carried.out", "done", "anderu.*durch"))

# Use it as a gate: if durchgeführt == 0, override AB_change to "No change"
if (!is.na(durchgefuehrt_col)) {
  executed_flag <- suppressWarnings(
    as.integer(as.character(zone_df_raw[[durchgefuehrt_col]]))
  )
  base_df$AB_change <- dplyr::if_else(
    executed_flag == 0 & !is.na(executed_flag),
    factor("No change", levels = levels(base_df$AB_change)),
    base_df$AB_change
  )
  cat(sprintf("Execution gate applied: %d rows overridden to 'No change'\n",
              sum(executed_flag == 0, na.rm = TRUE)))
}

cat(sprintf("ICU mortality: %d Yes / %d No / %d missing\n",
            sum(base_df$ICU_mortality==1L, na.rm=TRUE),
            sum(base_df$ICU_mortality==0L, na.rm=TRUE),
            sum(is.na(base_df$ICU_mortality))))

# ── 7b. Phase A — T1 baseline severity ───────────────────────────────────────
for (mk in markers) {
  zcols <- get_marker_cols(zone_df_raw, mk)
  if (!length(zcols)) {
    base_df[[paste0(mk,"_t1_zone")]] <- NA_character_
    base_df[[paste0(mk,"_t1_sev")]]  <- 0
    next
  }
  t1_zone <- purrr::map_chr(zone_df_raw[[zcols[1]]], extract_zone_text)
  base_df[[paste0(mk,"_t1_zone")]] <- t1_zone
  base_df[[paste0(mk,"_t1_sev")]]  <- zone_severity(t1_zone)
}

base_df <- base_df |>
  dplyr::mutate(
    admission_severity =
      PARAM$w_pct  * pct_t1_sev  + PARAM$w_saps * saps_t1_sev +
      PARAM$w_temp * temp_t1_sev + PARAM$w_leuk * leuk_t1_sev +
      PARAM$w_crp  * crp_t1_sev
  )

# ── 7c. Phase B — Zone severity at last timepoint ────────────────────────────
for (mk in markers) {
  zcols <- get_marker_cols(zone_df_raw, mk)
  if (!length(zcols)) {
    base_df[[paste0(mk,"_last_zone")]] <- NA_character_
    base_df[[paste0(mk,"_last_sev")]]  <- 0
    next
  }
  zmat <- zone_df_raw[zcols] |>
    dplyr::mutate(dplyr::across(dplyr::everything(),
                                ~purrr::map_chr(.x, extract_zone_text)))
  base_df[[paste0(mk,"_last_zone")]] <- apply(zmat, 1, last_nonmissing)
  base_df[[paste0(mk,"_last_sev")]]  <- zone_severity(base_df[[paste0(mk,"_last_zone")]])
}

# ── 7d. Phase B — Delta kinetics with MARKER-SPECIFIC scoring ────────────────
# Also store raw absolute delta for Temperature
for (mk in markers) {
  all_dcols <- get_marker_cols(delta_df_raw, mk)
  dcols <- all_dcols[stringr::str_detect(all_dcols,
                                         stringr::regex("t[12345].*t[12345]|delta|change|percent|t2|t3|t4|t5",
                                                        ignore_case=TRUE))]
  if (!length(dcols)) dcols <- all_dcols
  
  if (!length(dcols)) {
    for (sfx in c("_last_delta","_last_delta_score","_max_delta_score",
                  "_min_delta_score","_n_rising","_n_falling",
                  "_trajectory","_response_flag"))
      base_df[[paste0(mk,sfx)]] <- 0
    next
  }
  
  dmat      <- delta_df_raw[dcols] |>
    dplyr::mutate(dplyr::across(dplyr::everything(), extract_number))
  score_mat <- as.data.frame(lapply(dmat, function(col) delta_score_mk(col, mk)))
  
  base_df[[paste0(mk,"_last_delta")]]       <- apply(dmat, 1, last_nonmissing) |> as.numeric()
  base_df[[paste0(mk,"_last_delta_score")]] <- delta_score_mk(base_df[[paste0(mk,"_last_delta")]], mk)
  base_df[[paste0(mk,"_max_delta_score")]]  <- apply(score_mat, 1, \(x) max(x, na.rm=TRUE))
  base_df[[paste0(mk,"_min_delta_score")]]  <- apply(score_mat, 1, \(x) min(x, na.rm=TRUE))
  base_df[[paste0(mk,"_n_rising")]]         <- apply(score_mat, 1, \(x) sum(x>0, na.rm=TRUE))
  base_df[[paste0(mk,"_n_falling")]]        <- apply(score_mat, 1, \(x) sum(x<0, na.rm=TRUE))
  
  cuts <- MARKER_CUTS[[mk]]
  base_df[[paste0(mk,"_trajectory")]] <- dplyr::case_when(
    base_df[[paste0(mk,"_last_sev")]] == 0                                      ~ "normalised",
    base_df[[paste0(mk,"_last_delta_score")]] <= -1 &
      base_df[[paste0(mk,"_max_delta_score")]] <= 1                             ~ "improving",
    base_df[[paste0(mk,"_max_delta_score")]] >= 2                               ~ "worsening",
    TRUE                                                                         ~ "unclear"
  )
  base_df[[paste0(mk,"_response_flag")]] <- as.integer(
    base_df[[paste0(mk,"_trajectory")]] %in% c("improving","normalised")
  )
}

# ── 7e. Clinical modifiers ────────────────────────────────────────────────────
feature_df <- base_df |>
  dplyr::mutate(
    AB_change = factor(AB_change, levels=c("Stop","De-escalate","No change","Escalate")),
    Infect_clean      = tolower(dplyr::coalesce(Infect_type,"unknown")),
    infect_complexity = dplyr::case_when(
      stringr::str_detect(Infect_clean,"poly|mixed|polymicrobial|fung") ~ 2.0,
      stringr::str_detect(Infect_clean,"mono|single")                   ~ 1.0,
      stringr::str_detect(Infect_clean,"neg|none|unknown|no growth")    ~ 0.5,
      TRUE ~ 1.5),
    age_risk = dplyr::case_when(
      is.na(Age) ~ 0.3, Age>=80 ~ 1.0, Age>=70 ~ 0.5, TRUE ~ 0.0),
    n_markers_responding = pct_response_flag + crp_response_flag +
      leuk_response_flag + temp_response_flag + saps_response_flag,
    n_markers_normal =
      as.integer(pct_last_sev==0) + as.integer(crp_last_sev==0) +
      as.integer(leuk_last_sev==0) + as.integer(temp_last_sev==0) +
      as.integer(saps_last_sev==0)
  )

# ── 8. DATA-DERIVED THRESHOLD OPTIMISATION ───────────────────────────────────
# For markers where source = "Data-derived", find optimal delta threshold
# against ICU mortality using Youden index
# Updates MARKER_CUTS in place so all downstream scoring is consistent

cat("\n── Optimising data-derived thresholds vs ICU mortality ──\n")

threshold_validation <- purrr::map_dfr(markers, function(mk) {
  cuts    <- MARKER_CUTS[[mk]]
  delta_col <- paste0(mk, "_last_delta")
  
  sub <- feature_df |>
    dplyr::filter(!is.na(ICU_mortality), !is.na(.data[[delta_col]]))
  
  n_events <- sum(sub$ICU_mortality == 1L, na.rm=TRUE)
  
  # Need at least 10 events for stable ROC
  if (nrow(sub) < 30 || n_events < 10) {
    cat(sprintf("  %s: insufficient events (n=%d, events=%d) — keeping prior thresholds\n",
                mk, nrow(sub), n_events))
    return(tibble::tibble(
      Marker         = marker_labels[[mk]],
      marker_key     = mk,
      metric         = cuts$metric,
      unit           = cuts$unit,
      n_patients     = nrow(sub),
      n_icu_deaths   = n_events,
      AUC            = NA_real_,
      AUC_lower      = NA_real_,
      AUC_upper      = NA_real_,
      optimal_cut    = NA_real_,
      sensitivity    = NA_real_,
      specificity    = NA_real_,
      youden         = NA_real_,
      prior_mild_resp   = cuts$mild_resp,
      prior_strong_resp = cuts$strong_resp,
      prior_mild_alarm  = cuts$mild_alarm,
      prior_strong_alarm= cuts$strong_alarm,
      source_resp    = cuts$source_resp,
      source_alarm   = cuts$source_alarm,
      updated        = FALSE
    ))
  }
  
  ro <- pROC::roc(sub$ICU_mortality, sub[[delta_col]],
                  quiet=TRUE, levels=c(0,1), direction="auto")
  
  auc_ci  <- tryCatch(
    as.numeric(pROC::ci.auc(ro, conf.level=0.95, method="delong")),
    error = \(e) c(NA,NA,NA)
  )
  
  best <- tryCatch(
    pROC::coords(ro, "best",
                 ret=c("threshold","sensitivity","specificity","youden"),
                 best.method="youden", transpose=FALSE),
    error = \(e) data.frame(threshold=NA,sensitivity=NA,specificity=NA,youden=NA)
  )
  
  opt_cut <- if (!is.na(best$threshold[1])) round(best$threshold[1], 1) else NA_real_
  
  # Only update MARKER_CUTS for data-derived markers with meaningful AUC
  auc_val <- as.numeric(pROC::auc(ro))
  updated <- FALSE
  
  if (!is.na(opt_cut) && auc_val > 0.55) {
    if (stringr::str_detect(cuts$source_alarm, "Data")) {
      MARKER_CUTS[[mk]]$mild_alarm   <<- abs(opt_cut) * 0.6
      MARKER_CUTS[[mk]]$strong_alarm <<- abs(opt_cut)
      updated <- TRUE
    }
    if (stringr::str_detect(cuts$source_resp, "Data")) {
      MARKER_CUTS[[mk]]$mild_resp   <<- -abs(opt_cut) * 0.5
      MARKER_CUTS[[mk]]$strong_resp <<- -abs(opt_cut)
      updated <- TRUE
    }
    if (updated)
      cat(sprintf("  %s: updated thresholds from data (AUC=%.3f, Youden cut=%.1f)\n",
                  mk, auc_val, opt_cut))
  }
  
  tibble::tibble(
    Marker          = marker_labels[[mk]],
    marker_key      = mk,
    metric          = cuts$metric,
    unit            = cuts$unit,
    n_patients      = nrow(sub),
    n_icu_deaths    = n_events,
    AUC             = round(auc_val, 3),
    AUC_lower       = if (!is.na(auc_ci[1])) round(auc_ci[1], 3) else NA_real_,
    AUC_upper       = if (!is.na(auc_ci[3])) round(auc_ci[3], 3) else NA_real_,
    optimal_cut     = opt_cut,
    sensitivity     = round(best$sensitivity[1], 3),
    specificity     = round(best$specificity[1], 3),
    youden          = round(best$youden[1], 3),
    prior_mild_resp    = cuts$mild_resp,
    prior_strong_resp  = cuts$strong_resp,
    prior_mild_alarm   = cuts$mild_alarm,
    prior_strong_alarm = cuts$strong_alarm,
    source_resp     = cuts$source_resp,
    source_alarm    = cuts$source_alarm,
    updated         = updated
  )
})

cat("\nThreshold validation summary:\n")
print(threshold_validation |>
        dplyr::select(Marker, metric, unit, AUC, optimal_cut,
                      source_resp, source_alarm, updated))

# Rebuild delta scores with updated MARKER_CUTS
cat("\nRebuilding kinetic scores with final thresholds...\n")
for (mk in markers) {
  dcol  <- paste0(mk, "_last_delta")
  if (!dcol %in% names(feature_df)) next
  feature_df[[paste0(mk,"_last_delta_score")]] <-
    delta_score_mk(feature_df[[dcol]], mk)
}

# ── 9. COMPOSITE SCORES & PHASE C DECISION ───────────────────────────────────
feature_df <- feature_df |>
  dplyr::mutate(
    danger_score =
      PARAM$w_pct  * (pct_last_sev  + pmax(pct_last_delta_score,  0)) +
      PARAM$w_saps * (saps_last_sev + pmax(saps_last_delta_score, 0)) +
      PARAM$w_temp * (temp_last_sev + pmax(temp_last_delta_score, 0)) +
      PARAM$w_leuk * (leuk_last_sev + pmax(leuk_last_delta_score, 0)) +
      PARAM$w_crp  * (crp_last_sev  + pmax(crp_last_delta_score,  0)) +
      PARAM$w_mngs * mNGS_bin + PARAM$w_age * age_risk +
      PARAM$w_complexity * infect_complexity,
    
    response_score =
      PARAM$w_pct  * abs(pmin(pct_last_delta_score,  0)) +
      PARAM$w_saps * abs(pmin(saps_last_delta_score, 0)) +
      PARAM$w_temp * abs(pmin(temp_last_delta_score, 0)) +
      PARAM$w_leuk * abs(pmin(leuk_last_delta_score, 0)) +
      PARAM$w_crp  * abs(pmin(crp_last_delta_score,  0)),
    
    normality_score = n_markers_normal,
    
    hard_alarm_escalate =
      pct_max_delta_score  >= 2 | saps_max_delta_score >= 2 |
      leuk_last_sev >= 3        | saps_last_sev >= 3,
    
    hard_alarm_stop = hard_alarm_escalate | temp_last_sev >= 3,
    
    deesc_allowed =
      !hard_alarm_escalate &
      response_score >= PARAM$deesc_response_cut &
      saps_last_delta_score <= 0 &
      pct_last_delta_score  <= 0 &
      leuk_last_delta_score <= 1,
    
    stop_allowed =
      deesc_allowed & !hard_alarm_stop &
      danger_score   < PARAM$stop_danger_cut &
      response_score >= PARAM$stop_response_cut &
      (BK_bin==0 | n_markers_normal>=4) &
      mNGS_bin==0 & infect_complexity<=1.0 & n_markers_normal>=3,
    
    Fuzzy_class = dplyr::case_when(
      hard_alarm_escalate | danger_score >= PARAM$escalate_cut ~ "Escalate",
      stop_allowed                                              ~ "Stop",
      deesc_allowed & danger_score < 6                         ~ "De-escalate",
      TRUE                                                      ~ "No change"
    ) |> factor(levels=levels(AB_change)),
    
    raw_esc  = (danger_score - response_score*0.35 - 6.5) / 1.5,
    raw_desc = (response_score - danger_score*0.45 - 2.0 - BK_bin*0.5 - 0.5*mNGS_bin) / 1.2,
    raw_stop = (response_score + normality_score - danger_score - 2*BK_bin - 2*mNGS_bin - infect_complexity - 3.5) / 1.2,
    raw_noch = 0.0,
    sm_denom = exp(raw_esc)+exp(raw_desc)+exp(raw_stop)+exp(raw_noch),
    score_Escalate   = exp(raw_esc)  / sm_denom,
    score_DeEscalate = exp(raw_desc) / sm_denom,
    score_Stop       = exp(raw_stop) / sm_denom,
    score_NoChange   = exp(raw_noch) / sm_denom
  ) |>
  dplyr::select(-raw_esc,-raw_desc,-raw_stop,-raw_noch,-sm_denom)

cat(sprintf("\nFeature matrix: %d patients, %d columns\n", nrow(feature_df), ncol(feature_df)))
cat("Decision distribution (rule engine):\n")
print(table(feature_df$Fuzzy_class, useNA="ifany"))

# ── 10. LOOCV MULTINOMIAL ─────────────────────────────────────────────────────
predictors <- c(
  "admission_severity",
  "pct_t1_sev","saps_t1_sev","leuk_t1_sev","temp_t1_sev","crp_t1_sev",
  "danger_score","response_score","normality_score",
  "n_markers_responding","n_markers_normal",
  "pct_last_sev","pct_last_delta_score",
  "saps_last_sev","saps_last_delta_score",
  "leuk_last_sev","leuk_last_delta_score",
  "temp_last_sev","temp_last_delta_score",
  "crp_last_sev","crp_last_delta_score",
  "hard_alarm_escalate","BK_bin","mNGS_bin","Abx_bin",
  "Age","infect_complexity"
)

model_df <- feature_df |>
  dplyr::filter(!is.na(AB_change)) |>
  dplyr::mutate(
    hard_alarm_escalate = as.integer(hard_alarm_escalate),
    dplyr::across(dplyr::all_of(predictors),
                  ~tidyr::replace_na(as.numeric(.x), 0))
  )

valid_cls <- names(table(model_df$AB_change))[table(model_df$AB_change) >= 5]
model_df  <- model_df |>
  dplyr::filter(as.character(AB_change) %in% valid_cls) |>
  droplevels()

cat(sprintf("\nLOOCV: %d patients, classes: %s\n",
            nrow(model_df), paste(valid_cls, collapse=", ")))

loocv_predict <- function(df, outcome, preds) {
  lev      <- levels(df[[outcome]])
  prob_mat <- matrix(NA_real_, nrow(df), length(lev), dimnames=list(NULL,lev))
  pred_cls <- rep(NA_character_, nrow(df))
  
  for (i in seq_len(nrow(df))) {
    if (i %% 25==0) cat(sprintf("  LOOCV %d / %d\n", i, nrow(df)))
    train <- df[-i,,drop=FALSE]; test <- df[i,,drop=FALSE]
    keep  <- preds[sapply(train[preds], \(x) length(unique(x[!is.na(x)]))>1)]
    fit   <- tryCatch(
      nnet::multinom(as.formula(paste(outcome,"~",paste(keep,collapse="+"))),
                     data=train, trace=FALSE, MaxNWts=10000, maxit=1000),
      error=\(e) NULL)
    if (is.null(fit)) {
      pr_v <- setNames(prop.table(table(train[[outcome]]))|>as.numeric(), lev)
    } else {
      pr   <- predict(fit, newdata=test, type="probs")
      pr_v <- setNames(rep(0,length(lev)), lev)
      if (is.null(dim(pr))) pr_v[names(pr)] <- as.numeric(pr)
      else pr_v[colnames(pr)] <- as.numeric(pr[1,])
    }
    prob_mat[i,] <- pr_v
    pred_cls[i]  <- lev[which.max(pr_v)]
  }
  dplyr::bind_cols(
    df |> dplyr::select(JOB_ID, AB_change),
    tibble::as_tibble(prob_mat) |>
      dplyr::rename_with(~paste0("prob_",make.names(.x))),
    tibble::tibble(LOOCV_class=factor(pred_cls, levels=lev))
  )
}

set.seed(42)
cat("Running LOOCV...\n")
loocv_res  <- loocv_predict(model_df, "AB_change", predictors)
loocv_acc  <- mean(loocv_res$LOOCV_class==loocv_res$AB_change, na.rm=TRUE)
loocv_cm   <- caret::confusionMatrix(loocv_res$LOOCV_class, loocv_res$AB_change)
cat(sprintf("LOOCV accuracy: %.3f  Kappa: %.3f\n",
            loocv_acc, loocv_cm$overall["Kappa"]))

# ── 11. VARIABLE IMPORTANCE ───────────────────────────────────────────────────
model_std  <- model_df |>
  dplyr::mutate(dplyr::across(dplyr::all_of(predictors), scale))
keep_preds <- predictors[sapply(model_std[predictors],
                                \(x) length(unique(x[!is.na(x)]))>1)]
full_fit   <- tryCatch(
  nnet::multinom(as.formula(paste("AB_change~",paste(keep_preds,collapse="+"))),
                 data=model_std, trace=FALSE, MaxNWts=10000, maxit=1000),
  error=\(e){message("Full model failed: ",e$message); NULL})

importance_df <- if (!is.null(full_fit)) {
  tibble::tibble(
    Variable   = colnames(coef(full_fit)),
    Importance = colMeans(abs(coef(full_fit)), na.rm=TRUE)
  ) |>
    dplyr::filter(Variable!="(Intercept)") |>
    dplyr::arrange(desc(Importance)) |>
    dplyr::mutate(Group=dplyr::case_when(
      stringr::str_detect(Variable,"admission_severity|_t1_") ~ "Phase A baseline",
      stringr::str_detect(Variable,"danger|response|normality|n_markers") ~ "Phase B composite",
      stringr::str_detect(Variable,"pct")   ~ "PCT",
      stringr::str_detect(Variable,"saps")  ~ "SAPS II",
      stringr::str_detect(Variable,"leuk")  ~ "Leukocytes",
      stringr::str_detect(Variable,"temp")  ~ "Temperature",
      stringr::str_detect(Variable,"crp")   ~ "CRP",
      stringr::str_detect(Variable,"bk|mngs|abx|infect|complex") ~ "Microbiology / Tx",
      stringr::str_detect(Variable,"age|alarm") ~ "Demographics / alarms",
      TRUE ~ "Other"))
} else tibble::tibble(Variable="(model failed)", Importance=NA_real_, Group="Other")

# ── 12. ROC ───────────────────────────────────────────────────────────────────
cls_levels <- levels(model_df$AB_change)

roc_loocv <- purrr::map_dfr(cls_levels, function(cls) {
  pcol  <- paste0("prob_",make.names(cls))
  if (!pcol %in% names(loocv_res)) return(NULL)
  truth <- as.integer(loocv_res$AB_change==cls)
  if (length(unique(truth))<2) return(NULL)
  ro <- pROC::roc(truth, loocv_res[[pcol]], quiet=TRUE, levels=c(0,1), direction="<")
  tibble::tibble(FPR=1-ro$specificities, TPR=ro$sensitivities,
                 Class=cls, AUC=round(as.numeric(pROC::auc(ro)),3),
                 Model="LOOCV Multinomial")
})

score_cols <- c(Stop="score_Stop",`De-escalate`="score_DeEscalate",
                `No change`="score_NoChange", Escalate="score_Escalate")

roc_fuzzy <- purrr::map_dfr(cls_levels, function(cls) {
  scol  <- score_cols[cls]
  if (is.na(scol)||!scol %in% names(feature_df)) return(NULL)
  sub   <- feature_df |> dplyr::filter(!is.na(AB_change))
  truth <- as.integer(sub$AB_change==cls)
  if (length(unique(truth))<2) return(NULL)
  ro <- pROC::roc(truth, sub[[scol]], quiet=TRUE, levels=c(0,1), direction="<")
  tibble::tibble(FPR=1-ro$specificities, TPR=ro$sensitivities,
                 Class=cls, AUC=round(as.numeric(pROC::auc(ro)),3),
                 Model="Fuzzy Rule Engine")
})

roc_all <- dplyr::bind_rows(roc_loocv, roc_fuzzy)
auc_tbl <- roc_all |> dplyr::distinct(Class,AUC,Model) |> dplyr::arrange(Model,desc(AUC))
print(auc_tbl)

# ================================================================
# 13. FIGURES
# ================================================================
pal4 <- c("No change"="#4393C3","Escalate"="#D6604D",
          "De-escalate"="#4DAC26","Stop"="#E08214")
pal_traj <- c("improving"="#4DAC26","normalised"="#A6D96A",
              "unclear"="#FFEDA0","worsening"="#D6604D")
pal_diag <- c(
  "mNGS+/BC-/Abx+" = "#1F4E79",
  "mNGS+/BC+/Abx+" = "#2E75B6",
  "mNGS+/BC-/Abx-" = "#9DC3E6",
  "mNGS+/BC+/Abx-" = "#BDD7EE",
  "mNGS-/BC-/Abx+" = "#843C0C",
  "mNGS-/BC+/Abx+" = "#C55A11",
  "Other"           = "#AAAAAA"
)
pal_grp  <- c("Phase A baseline"="#185FA5","Phase B composite"="#7B2D8B",
              "PCT"="#D6604D","SAPS II"="#639922","Leukocytes"="#BA7517",
              "Temperature"="#7F77DD","CRP"="#E24B4A",
              "Microbiology / Tx"="#888780",
              "Demographics / alarms"="#5DCAA5","Other"="#B4B2A9")

theme_pub <- function(base=11)
  ggplot2::theme_bw(base_size=base) +
  ggplot2::theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill="grey92",colour=NA),
    legend.position   = "bottom",
    plot.title        = element_text(face="bold",size=base+1),
    plot.subtitle     = element_text(size=base-1,colour="grey40"),
    axis.title        = element_text(size=base-1)
  )

save_fig <- function(p, name, w, h) {
  # cairo_pdf supports Unicode (arrows, Greek letters, en-dashes)
  # Falls back to standard pdf if cairo is unavailable
  pdf_dev <- tryCatch(grDevices::cairo_pdf, error = function(e) grDevices::pdf)
  ggsave(paste0("Figures/",name,".pdf"), p, width=w, height=h, device=pdf_dev)
  ggsave(paste0("Figures/",name,".png"), p, width=w, height=h, dpi=300)
  cat(sprintf("Saved %s\n", name))
}

# ── Fig 1 · Clinical pathway ──────────────────────────────────────────────────
phase_boxes <- tibble::tribble(
  ~x,   ~y,   ~label,                                        ~fill,     ~phase,
  1.0,  4.0,  "PHASE A\nDay 0 / Admission",                  "#1F4E79", "hdr",
  1.0,  3.2,  "Baseline severity\n(T1: CRP, PCT, Leuk,\nTemp, SAPS II)", "#BDD7EE","A",
  1.0,  2.2,  "Current antibiotics\n(ABX on/off)",            "#BDD7EE","A",
  1.0,  1.3,  "Suspected diagnosis\n+ cultures / NGS",        "#BDD7EE","A",
  3.0,  4.0,  "PHASE B\nDay 1\u20134 / Monitoring",           "#375623","hdr",
  3.0,  3.2,  "Biomarker kinetics\n(\u0394 T1\u2192T5)",      "#E2EFDA","B",
  3.0,  2.2,  "Microbiology results\n(BK pos/neg, NGS)",      "#E2EFDA","B",
  3.0,  1.3,  "Clinical response\n(improving / not)",         "#E2EFDA","B",
  5.0,  4.0,  "PHASE C\nDecision",                            "#7F3F00","hdr",
  5.0,  3.4,  "Escalate",                                     "#FFC7CE","C",
  5.0,  2.7,  "No change",                                    "#BDD7EE","C",
  5.0,  2.0,  "De-escalate",                                  "#E2EFDA","C",
  5.0,  1.3,  "Stop",                                         "#FFEB9C","C",
  7.0,  4.0,  "PHASE D\nOutcome",                             "#4C0022","hdr",
  7.0,  2.8,  "ICU mortality\n(ICU days = Hosp days\n& Died = Yes)", "#FFC7CE","D",
  7.0,  1.5,  "Hospital survival\n& discharge",               "#E2EFDA","D"
)
fig1 <- ggplot(phase_boxes, aes(x=x,y=y)) +
  geom_tile(aes(fill=fill),
            width=1.55,
            height=dplyr::if_else(phase_boxes$phase=="hdr",0.45,0.72),
            colour="grey55", linewidth=0.4) +
  geom_text(aes(label=label,
                colour=dplyr::if_else(phase=="hdr","#FFFFFF","#1A1A1A")),
            size=2.75, lineheight=0.92) +
  scale_fill_identity() + scale_colour_identity() +
  annotate("segment",x=1.8,xend=2.2,y=2.5,yend=2.5,
           arrow=arrow(length=unit(0.1,"inches"),type="closed"),colour="grey35",linewidth=0.8) +
  annotate("segment",x=3.8,xend=4.2,y=2.5,yend=2.5,
           arrow=arrow(length=unit(0.1,"inches"),type="closed"),colour="grey35",linewidth=0.8) +
  annotate("segment",x=5.8,xend=6.2,y=2.5,yend=2.5,
           arrow=arrow(length=unit(0.1,"inches"),type="closed"),colour="grey35",linewidth=0.8) +
  scale_x_continuous(limits=c(0.1,7.9)) + scale_y_continuous(limits=c(0.7,4.5)) +
  labs(title="Figure 1 \u2014 Clinical pathway: A\u2192B\u2192C\u2192D",
       subtitle="Admission severity \u2192 Biomarker kinetics \u2192 ABX decision \u2192 ICU mortality") +
  theme_void() +
  theme(plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=10,colour="grey40"),
        plot.margin=margin(12,12,12,12))
save_fig(fig1, "Fig1_pathway", 12, 5.5)

# ── Fig 2 · Marker-specific threshold reference ───────────────────────────────
# Shows fuzzy membership curves with marker-specific cuts
# Panel A: severity zones (absolute values)
# Panel B: kinetic thresholds (delta% or absolute)
# Panel C: threshold table with sources

# Panel A — severity membership (same across all versions)
sev_curves <- purrr::map_dfr(list(
  list(mk="crp",  x=seq(0,300,1),  uln=5,   alarm=50),
  list(mk="pct",  x=seq(0,20,0.1), uln=0.5, alarm=2),
  list(mk="leuk", x=seq(0,25,0.1), uln=10,  lln=4, alarm=20),
  list(mk="temp", x=seq(35,41,0.1),uln=38.3,lln=36,alarm=40),
  list(mk="saps", x=seq(0,80,1),   uln=29,  alarm=60)
), function(p) {
  tibble::tibble(
    Marker = marker_labels[[p$mk]],
    x = p$x,
    Normal = if (!is.null(p$lln) && !is.na(p$lln))
      pmax(0, pmin(1, (p$x - p$lln) / (p$uln - p$lln)))
    else
      pmax(0, 1 - (p$x - p$uln) / p$uln),
    # High alarm (leukocytosis, high CRP, high PCT, high temp, high SAPS)
    Alarm_high = pmax(0, pmin(1, (p$x - p$uln) / (p$alarm - p$uln))),
    # Leukopenia alarm: only for leukocytes — rises as value falls below lln=4
    # severe leukopenia threshold = 2 (from zone definitions)
    Alarm_low  = if (!is.null(p$lln) && !is.na(p$lln) && p$mk == "leuk")
      pmax(0, pmin(1, (p$lln - p$x) / (p$lln - 2)))
    else
      rep(0, length(p$x))
  )
})

figA <- ggplot(sev_curves) +
  geom_line(aes(x=x, y=Normal),     colour="#4DAC26", linewidth=1.0) +
  geom_line(aes(x=x, y=Alarm_high), colour="#D6604D", linewidth=1.0) +
  # Leukopenia alarm — red line on left side for Leukocytes only
  geom_line(aes(x=x, y=Alarm_low),  colour="#D6604D", linewidth=1.0,
            linetype="solid") +
  facet_wrap(~Marker, scales="free_x", nrow=1) +
  scale_y_continuous(limits=c(0,1), labels=scales::percent_format()) +
  labs(title="A  Severity membership (absolute values)",
       subtitle="Leukocytes: bidirectional alarm (leukocytosis right, leukopenia left)",
       x="Raw biomarker value", y="Membership (\u03bc)") +
  theme_pub(base=9) +
  theme(strip.text=element_text(face="bold"))

# Panel B — kinetic threshold visualisation
kinetic_bars <- threshold_validation |>
  dplyr::filter(!is.na(AUC)) |>
  dplyr::mutate(
    Marker = factor(Marker, levels=rev(c("CRP","PCT","Leukocytes","Temperature","SAPS II"))),
    mild_resp_val   = purrr::map2_dbl(marker_key, prior_mild_resp,
                                      ~MARKER_CUTS[[.x]]$mild_resp),
    strong_resp_val = purrr::map2_dbl(marker_key, prior_strong_resp,
                                      ~MARKER_CUTS[[.x]]$strong_resp),
    mild_alarm_val  = purrr::map2_dbl(marker_key, prior_mild_alarm,
                                      ~MARKER_CUTS[[.x]]$mild_alarm),
    strong_alarm_val= purrr::map2_dbl(marker_key, prior_strong_alarm,
                                      ~MARKER_CUTS[[.x]]$strong_alarm),
    auc_lab = sprintf("AUC=%.2f [%.2f\u2013%.2f]", AUC, AUC_lower, AUC_upper)
  )

# Leukocytes immunoparalysis alarm zone (extreme fall on negative side)
leuk_immuno <- tibble::tibble(
  Marker     = factor("Leukocytes",
                      levels = rev(c("CRP","PCT","Leukocytes","Temperature","SAPS II"))),
  x_start    = -60,    # immunoparalysis threshold starts at -60%
  x_end      = -100    # extends left (worse = more negative)
)

figB <- ggplot(kinetic_bars, aes(y=Marker)) +
  # Response zone (green — improvement signal)
  geom_segment(aes(x=strong_resp_val, xend=mild_resp_val,
                   yend=Marker), colour="#4DAC26", linewidth=8, alpha=0.5) +
  # Alarm zone — rising (red, right side)
  geom_segment(aes(x=mild_alarm_val, xend=strong_alarm_val,
                   yend=Marker), colour="#D6604D", linewidth=8, alpha=0.5) +
  # Leukocytes immunoparalysis alarm — extreme fall (red, left side)
  geom_segment(data=leuk_immuno,
               aes(x=x_start, xend=x_end, y=Marker, yend=Marker),
               colour="#D6604D", linewidth=8, alpha=0.5) +
  # Vertical red boundary line at -60% showing threshold
  annotate("segment", x=-60, xend=-60, y=2.7, yend=3.3,
           colour="#D6604D", linewidth=1.5) +
  # Immunoparalysis label
  annotate("text", x=-80, y="Leukocytes", vjust=-1.0, hjust=0.5,
           label="Immunoparalysis
(>60% fall)", size=2.4,
           colour="#D6604D", fontface="italic") +
  # Optimal cut from data
  geom_point(aes(x=optimal_cut), colour="#1A1A1A", size=4, shape=18) +
  # AUC label
  geom_text(aes(x=optimal_cut, label=auc_lab),
            vjust=-1.4, hjust=0.5, size=2.6, colour="grey30",
            check_overlap=TRUE) +
  geom_vline(xintercept=0, linetype="dashed", colour="grey50", linewidth=0.5) +
  annotate("text", x=-45, y=0.6, label="Response", colour="#4DAC26",
           fontface="bold", size=3) +
  annotate("text", x=45,  y=0.6, label="Alarm",    colour="#D6604D",
           fontface="bold", size=3) +
  labs(title="B  Kinetic thresholds by marker",
       subtitle="Green = response zone; Red = alarm zone (rise + immunoparalysis for Leukocytes); Diamond = data-derived Youden cut",
       x="Delta% threshold (or \u0394\u00b0C for Temperature)", y=NULL) +
  theme_pub(base=9)

# Panel C — threshold source table
thresh_tbl <- threshold_validation |>
  dplyr::transmute(
    Marker,
    Metric = dplyr::if_else(metric=="delta_abs","\u0394°C (absolute)","\u0394%"),
    `Mild resp` = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_resp),
    `Strong resp` = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_resp),
    `Mild alarm` = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$mild_alarm),
    `Strong alarm` = purrr::map_dbl(marker_key, ~MARKER_CUTS[[.x]]$strong_alarm),
    `Source (resp)` = source_resp,
    `Source (alarm)` = source_alarm,
    `AUC vs ICU mort` = dplyr::if_else(!is.na(AUC), sprintf("%.3f",AUC), "n/a"),
    Updated = dplyr::if_else(updated,"Yes","No")
  )

figC <- gridExtra::tableGrob(
  thresh_tbl,
  rows = NULL,
  theme = gridExtra::ttheme_minimal(
    core    = list(fg_params=list(cex=0.65)),
    colhead = list(fg_params=list(cex=0.70, fontface="bold"))
  )
)

fig2 <- (figA / figB) +
  plot_annotation(
    title    = "Figure 2 \u2014 Marker-specific kinetic thresholds",
    subtitle = "A: severity zones  |  B: kinetic response/alarm zones with data-derived validation",
    theme    = theme(plot.title=element_text(face="bold",size=13),
                     plot.subtitle=element_text(size=10,colour="grey40"))
  )
save_fig(fig2, "Fig2_thresholds", 14, 10)

# Save threshold table separately as clean PDF
if (requireNamespace("gridExtra", quietly=TRUE)) {
  grDevices::pdf("Figures/Fig2b_threshold_table.pdf", width=14, height=4)
  gridExtra::grid.arrange(figC)
  grDevices::dev.off()
  cat("Saved Fig2b_threshold_table\n")
}

# ── Fig 3 · Admission severity by decision ────────────────────────────────────
fig3 <- feature_df |>
  dplyr::filter(!is.na(AB_change)) |>
  ggplot(aes(x=AB_change,y=admission_severity,fill=AB_change,colour=AB_change)) +
  geom_violin(alpha=0.25,linewidth=0.5,trim=TRUE) +
  geom_jitter(alpha=0.40,size=1.2,width=0.15) +
  geom_boxplot(width=0.14,fill="white",alpha=0.85,outlier.shape=NA,linewidth=0.5) +
  scale_fill_manual(values=pal4,guide="none") +
  scale_colour_manual(values=pal4,guide="none") +
  scale_x_discrete(guide=guide_axis(angle=20)) +
  labs(title="Figure 3 \u2014 Phase A: Admission severity by decision class",
       subtitle="Weighted composite at Day 0 (T1)",
       x=NULL, y="Admission severity score") +
  theme_pub()
save_fig(fig3, "Fig3_admission_severity", 8, 5.5)

# ── Fig 4 · Kinetics trajectory heatmap ───────────────────────────────────────
traj_long <- feature_df |>
  dplyr::filter(!is.na(AB_change)) |>
  dplyr::select(JOB_ID, AB_change, dplyr::ends_with("_trajectory")) |>
  tidyr::pivot_longer(dplyr::ends_with("_trajectory"),
                      names_to="Marker", values_to="Trajectory") |>
  dplyr::mutate(
    Marker = stringr::str_remove(Marker,"_trajectory") |>
      stringr::str_to_upper() |>
      dplyr::recode(LEUK="Leukocytes",TEMP="Temperature",SAPS="SAPS II"),
    Marker     = factor(Marker, levels=c("CRP","PCT","Leukocytes","SAPS II","Temperature")),
    Trajectory = factor(Trajectory, levels=c("normalised","improving","unclear","worsening"))
  )

pt_order4 <- feature_df |>
  dplyr::filter(!is.na(AB_change)) |>
  dplyr::arrange(AB_change, desc(n_markers_responding)) |>
  dplyr::pull(JOB_ID)

fig4 <- traj_long |>
  dplyr::mutate(JOB_ID=factor(JOB_ID, levels=rev(pt_order4))) |>
  ggplot(aes(x=Marker, y=JOB_ID, fill=Trajectory)) +
  geom_tile(colour="white",linewidth=0.15) +
  facet_grid(AB_change~., scales="free_y", space="free_y") +
  scale_fill_manual(values=pal_traj, na.value="grey85", name="Trajectory") +
  scale_x_discrete(guide=guide_axis(angle=30)) +
  labs(title="Figure 4 \u2014 Phase B: Biomarker kinetics trajectory (T1\u2192T5)",
       subtitle="Marker-specific thresholds applied; sorted by decision class + response count",
       x=NULL, y="Patient (JOB_ID)") +
  theme_pub() +
  theme(axis.text.y=element_text(size=4.5),
        strip.text.y=element_text(angle=0,face="bold"),
        legend.position="right",
        panel.spacing.y=unit(0.3,"lines"))
save_fig(fig4, "Fig4_kinetics_heatmap", 9, 12)

# ── Fig 5 · Delta% distribution per marker vs ICU mortality ───────────────────
delta_mort_long <- feature_df |>
  dplyr::filter(!is.na(ICU_mortality)) |>
  dplyr::select(JOB_ID, ICU_mortality,
                pct_last_delta, crp_last_delta,
                leuk_last_delta, temp_last_delta, saps_last_delta) |>
  tidyr::pivot_longer(dplyr::ends_with("_last_delta"),
                      names_to="Marker_raw", values_to="Delta") |>
  dplyr::mutate(
    Marker = stringr::str_remove(Marker_raw,"_last_delta") |>
      stringr::str_to_upper() |>
      dplyr::recode(LEUK="Leukocytes",TEMP="Temperature",SAPS="SAPS II"),
    Marker = factor(Marker, levels=c("CRP","PCT","Leukocytes","SAPS II","Temperature")),
    Outcome = factor(ICU_mortality, levels=c(0,1), labels=c("Survived","ICU death"))
  ) |>
  dplyr::filter(!is.na(Delta))

# Add marker-specific threshold lines
thresh_lines <- tibble::tribble(
  ~Marker,        ~thresh,  ~label,
  "CRP",          -25,      "Resp (-25%)",
  "CRP",          -50,      "Strong resp (-50%)",
  "PCT",          -25,      "Resp (-25%)",
  "PCT",          -50,      "Strong resp (-50%)",
  "Leukocytes",   -30,      "Resp (-30%)",
  "Leukocytes",    30,      "Alarm (+30%)",
  "Temperature",   -1.0,   "Resp (-1.0\u00b0C)",
  "Temperature",    1.0,   "Alarm (+1.0\u00b0C)",
  "SAPS II",      -20,      "Resp (-20%)",
  "SAPS II",       15,      "Alarm (+15%)"
) |>
  dplyr::mutate(Marker=factor(Marker, levels=c("CRP","PCT","Leukocytes","SAPS II","Temperature")))

fig5_pvals <- delta_mort_long |>
  dplyr::filter(!is.na(Delta), !is.na(ICU_mortality)) |>
  dplyr::group_by(Marker) |>
  dplyr::summarise(
    p_val = tryCatch(wilcox.test(Delta ~ ICU_mortality)$p.value, error=\(e) NA_real_),
    y_pos = quantile(Delta, 0.97, na.rm=TRUE),
    .groups="drop"
  ) |>
  dplyr::mutate(
    p_lab = dplyr::case_when(
      is.na(p_val)  ~ "n/a",
      p_val < 0.001 ~ "p<0.001",
      p_val < 0.05  ~ sprintf("p=%.3f", p_val),
      TRUE          ~ sprintf("p=%.2f NS", p_val)
    )
  )

fig5 <- ggplot(delta_mort_long,
               aes(x=Outcome, y=Delta, fill=Outcome, colour=Outcome)) +
  geom_violin(alpha=0.25, linewidth=0.5, trim=TRUE) +
  geom_jitter(alpha=0.35, size=0.8, width=0.15) +
  geom_boxplot(width=0.14, fill="white", alpha=0.85,
               outlier.shape=NA, linewidth=0.5) +
  geom_hline(data=thresh_lines, aes(yintercept=thresh),
             linetype="dashed", colour="grey30", linewidth=0.45) +
  geom_text(data=thresh_lines, aes(y=thresh, label=label),
            x=2.48, hjust=0, size=2.0, colour="grey30",
            inherit.aes=FALSE) +
  geom_text(data=fig5_pvals,
            aes(x=1.5, y=y_pos, label=p_lab, fill=NULL, colour=NULL),
            size=2.7, colour="grey25", fontface="italic",
            inherit.aes=FALSE) +
  facet_wrap(~Marker, scales="free_y", nrow=1) +
  scale_fill_manual(values=c("Survived"="#4DAC26","ICU death"="#D6604D"),
                    guide="none") +
  scale_colour_manual(values=c("Survived"="#4DAC26","ICU death"="#D6604D"),
                      guide="none") +
  scale_x_discrete(guide=guide_axis(angle=20)) +
  labs(title="Figure 5 \u2014 Delta% at last timepoint vs ICU mortality, by marker",
       subtitle="Dashed lines = marker-specific thresholds  |  Wilcoxon p-values shown",
       x=NULL, y="Delta at last timepoint (% or \u0394\u00b0C)") +
  theme_pub() +
  theme(panel.spacing.x=unit(1.2,"lines"))
save_fig(fig5, "Fig5_delta_vs_mortality", 14, 7)

# ── Fig 6 · ROC: each marker delta vs ICU mortality ───────────────────────────
roc_mk_df <- purrr::map_dfr(markers, function(mk) {
  dcol  <- paste0(mk,"_last_delta")
  if (!dcol %in% names(feature_df)) return(NULL)
  sub   <- feature_df |> dplyr::filter(!is.na(ICU_mortality),!is.na(.data[[dcol]]))
  if (nrow(sub)<20||length(unique(sub$ICU_mortality))<2) return(NULL)
  ro <- pROC::roc(sub$ICU_mortality, sub[[dcol]], quiet=TRUE,
                  levels=c(0,1), direction="auto")
  auc_ci <- tryCatch(as.numeric(pROC::ci.auc(ro,conf.level=0.95)),
                     error=\(e) c(NA,NA,NA))
  tibble::tibble(
    FPR   = 1-ro$specificities,
    TPR   = ro$sensitivities,
    Marker = marker_labels[[mk]],
    AUC   = round(as.numeric(pROC::auc(ro)),3),
    AUC_lo= if(!is.na(auc_ci[1])) round(auc_ci[1],3) else NA_real_,
    AUC_hi= if(!is.na(auc_ci[3])) round(auc_ci[3],3) else NA_real_
  )
})

auc_labs6 <- roc_mk_df |>
  dplyr::distinct(Marker,AUC,AUC_lo,AUC_hi) |>
  dplyr::mutate(lab=sprintf("%s\nAUC=%.3f [%.3f\u2013%.3f]",Marker,AUC,AUC_lo,AUC_hi))

fig6 <- ggplot(roc_mk_df, aes(FPR,TPR,colour=Marker)) +
  geom_line(linewidth=1.0) +
  geom_abline(slope=1,intercept=0,linetype="dashed",colour="grey65",linewidth=0.4) +
  geom_label(data=auc_labs6, aes(x=0.65,y=0.25,label=lab,colour=Marker),
             size=2.5, fill="white", label.padding=unit(0.2,"lines"),
             show.legend=FALSE) +
  scale_colour_brewer(palette="Set1", name="Marker") +
  scale_x_continuous(labels=percent_format(),limits=c(0,1)) +
  scale_y_continuous(labels=percent_format(),limits=c(0,1)) +
  coord_equal() +
  facet_wrap(~Marker, nrow=1) +
  labs(title="Figure 6 \u2014 ROC: individual marker kinetics vs ICU mortality",
       subtitle="Delta at last timepoint as predictor; validates marker-specific thresholds",
       x="1 \u2212 Specificity", y="Sensitivity") +
  theme_pub() + theme(legend.position="none")
save_fig(fig6, "Fig6_roc_markers_mortality", 14, 5)

# ── Fig 7–12: standard figures from v5 (decision scatter, LOOCV CM, etc.) ────
# Fig 7 · Danger vs response scatter
fig7 <- feature_df |>
  dplyr::filter(!is.na(AB_change)) |>
  ggplot(aes(x=danger_score,y=response_score,colour=AB_change,shape=AB_change)) +
  geom_jitter(alpha=0.75,size=2.0,width=0.15,height=0.10) +
  geom_vline(xintercept=PARAM$escalate_cut,linetype="dashed",colour="grey35",linewidth=0.55) +
  geom_hline(yintercept=PARAM$deesc_response_cut,linetype="dotted",colour="grey35",linewidth=0.55) +
  scale_colour_manual(values=pal4,name="True decision") +
  scale_shape_manual(values=c(16,17,15,18),name="True decision") +
  facet_wrap(~AB_change,nrow=2) +
  labs(title="Figure 7 \u2014 Phase C: Danger vs response score",
       x="Danger score",y="Response score") +
  theme_pub() + theme(legend.position="none")
save_fig(fig7, "Fig7_scatter", 9, 7)

# Fig 8 · LOOCV ROC
auc_labels8 <- roc_all |>
  dplyr::distinct(Class, AUC, Model) |>
  dplyr::mutate(
    lab   = sprintf("AUC=%.3f", AUC),
    x_pos = 0.68,
    y_pos = dplyr::if_else(Model=="LOOCV Multinomial", 0.16, 0.06)
  )

fig8 <- ggplot(roc_all, aes(FPR, TPR, colour=Class, linetype=Model)) +
  geom_line(linewidth=0.9) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              colour="grey65", linewidth=0.4) +
  geom_text(data=auc_labels8,
            aes(x=x_pos, y=y_pos, label=lab, colour=Class),
            size=2.7, fontface="bold", show.legend=FALSE,
            inherit.aes=FALSE) +
  scale_colour_manual(values=pal4, name="Decision class") +
  scale_linetype_manual(
    values=c("LOOCV Multinomial"="solid","Fuzzy Rule Engine"="dashed"),
    name="Model") +
  scale_x_continuous(labels=percent_format(), limits=c(0,1)) +
  scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
  coord_equal() +
  facet_wrap(~Class, nrow=2) +
  labs(title="Figure 8 \u2014 Phase C: One-vs-rest ROC curves",
       subtitle="Solid = LOOCV Multinomial  |  Dashed = Fuzzy Rule Engine  |  AUC values labelled",
       x="1 \u2212 Specificity", y="Sensitivity") +
  theme_pub() +
  theme(legend.position="right")
save_fig(fig8, "Fig8_roc_decision", 11, 8)

# Fig 9 · Confusion matrix
cm_df <- as.data.frame(loocv_cm$table) |>
  setNames(c("Predicted","Reference","Freq")) |>
  dplyr::group_by(Reference) |>
  dplyr::mutate(RowPct=Freq/sum(Freq)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    Predicted=factor(Predicted,levels=rev(levels(model_df$AB_change))),
    Reference=factor(Reference,    levels=    levels(model_df$AB_change))
  )
fig9 <- ggplot(cm_df,aes(x=Reference,y=Predicted,fill=RowPct)) +
  geom_tile(colour="white",linewidth=0.5) +
  geom_text(aes(label=paste0(Freq,"\n(",percent(RowPct,accuracy=1),")")),
            size=3.5,colour="grey10") +
  scale_fill_gradient(low="#F7F7F7",high="#2166AC",labels=percent_format(),name="Row %") +
  scale_x_discrete(guide=guide_axis(angle=25)) +
  labs(title="Figure 9 \u2014 LOOCV confusion matrix",
       subtitle=sprintf("Accuracy: %.1f%%  |  Kappa: %.3f",loocv_acc*100,loocv_cm$overall["Kappa"]),
       x="True label",y="Predicted") +
  theme_pub()
save_fig(fig9, "Fig9_confusion", 7, 6)

# Fig 10 · Variable importance
fig10 <- importance_df |>
  dplyr::slice_head(n=18) |>
  dplyr::mutate(VarLabel=stringr::str_replace_all(Variable,"_"," ")|>stringr::str_to_sentence()) |>
  ggplot(aes(x=reorder(VarLabel,Importance),y=Importance,fill=Group)) +
  geom_col(width=0.72) +
  geom_text(aes(label=round(Importance,2)),hjust=-0.15,size=2.9) +
  scale_fill_manual(values=pal_grp,name="Feature group") +
  scale_y_continuous(expand=expansion(mult=c(0,0.15))) +
  coord_flip() +
  labs(title="Figure 10 \u2014 Variable importance (standardised)",
       x=NULL,y="Mean |standardised coefficient|") +
  theme_pub() + theme(legend.position="right",panel.grid.major.y=element_blank())
save_fig(fig10, "Fig10_importance", 11, 7)

# Fig 11 · ICU mortality by decision
mort_df <- feature_df |>
  dplyr::filter(!is.na(AB_change),!is.na(ICU_mortality)) |>
  dplyr::group_by(AB_change) |>
  dplyr::summarise(
    n_total=dplyr::n(), n_icu_death=sum(ICU_mortality==1L),
    pct_mort=n_icu_death/n_total,
    se=sqrt(pct_mort*(1-pct_mort)/n_total),
    ci_lo=pmax(pct_mort-1.96*se,0), ci_hi=pmin(pct_mort+1.96*se,1),
    .groups="drop")

fig11 <- ggplot(mort_df, aes(x=AB_change, y=pct_mort, fill=AB_change)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi),
                width=0.22, linewidth=0.7) +
  geom_text(aes(label=sprintf("%d/%d\n(%.0f%%)", n_icu_death, n_total, pct_mort*100),
                y=ci_hi + 0.004),
            vjust=0, size=3.0, colour="grey20", lineheight=0.85) +
  scale_fill_manual(values=pal4, guide="none") +
  scale_y_continuous(labels=percent_format(),
                     expand=expansion(mult=c(0, 0.28))) +
  labs(title="Figure 11 \u2014 Phase D: ICU mortality by ABX decision class",
       subtitle="Error bars = 95% CI  |  Labels positioned above error bars",
       x=NULL, y="ICU mortality rate") +
  theme_pub() +
  theme(panel.grid.major.x=element_blank())
save_fig(fig11, "Fig11_icu_mortality", 7, 6.5)

# Fig 12 · Danger score vs mortality
fig12 <- feature_df |>
  dplyr::filter(!is.na(ICU_mortality)) |>
  dplyr::mutate(Outcome=factor(ICU_mortality,levels=c(0,1),
                               labels=c("Survived","ICU death"))) |>
  ggplot(aes(x=Outcome,y=danger_score,fill=Outcome,colour=Outcome)) +
  geom_violin(alpha=0.25,linewidth=0.5,trim=TRUE) +
  geom_jitter(alpha=0.4,size=1.5,width=0.15) +
  geom_boxplot(width=0.14,fill="white",alpha=0.85,outlier.shape=NA,linewidth=0.5) +
  scale_fill_manual(values=c("Survived"="#4DAC26","ICU death"="#D6604D"),guide="none") +
  scale_colour_manual(values=c("Survived"="#4DAC26","ICU death"="#D6604D"),guide="none") +
  labs(title="Figure 12 \u2014 Danger score vs ICU mortality",
       x=NULL,y="Danger score") +
  theme_pub()
save_fig(fig12, "Fig12_danger_mortality", 6.5, 5.5)

# Fig 13 · Alluvial
alluv_df <- feature_df |>
  dplyr::filter(!is.na(AB_change),!is.na(ICU_mortality)) |>
  dplyr::mutate(Outcome=factor(ICU_mortality,levels=c(0,1),
                               labels=c("Survived","ICU death"))) |>
  dplyr::count(AB_change,Fuzzy_class,Outcome)

fig13 <- ggplot(alluv_df,
                aes(axis1=AB_change,axis2=Fuzzy_class,axis3=Outcome,y=n)) +
  geom_alluvium(aes(fill=AB_change),alpha=0.65,width=0.2) +
  geom_stratum(width=0.28,fill="grey92",colour="grey55") +
  geom_label(stat="stratum",aes(label=after_stat(stratum)),size=2.9) +
  scale_x_discrete(limits=c("Workbook decision","Rule engine","Outcome"),
                   expand=c(0.08,0.05)) +
  scale_fill_manual(values=pal4,guide="none") +
  labs(title="Figure 13 \u2014 Alluvial: workbook \u2192 rule engine \u2192 ICU outcome",
       y="Number of patients") +
  theme_pub() + theme(axis.text.x=element_text(face="bold",size=11))
save_fig(fig13, "Fig13_alluvial", 10, 7)

# ── 13a. FINAL THRESHOLD TABLE (needed by publication stats below) ───────────
final_thresh_tbl <- tibble::tibble(
  Marker        = names(marker_labels),
  Marker_label  = unname(marker_labels),
  Metric        = purrr::map_chr(names(marker_labels), ~MARKER_CUTS[[.x]]$metric),
  Unit          = purrr::map_chr(names(marker_labels), ~MARKER_CUTS[[.x]]$unit),
  Mild_response = purrr::map_dbl(names(marker_labels), ~MARKER_CUTS[[.x]]$mild_resp),
  Strong_response= purrr::map_dbl(names(marker_labels),~MARKER_CUTS[[.x]]$strong_resp),
  Mild_alarm    = purrr::map_dbl(names(marker_labels), ~MARKER_CUTS[[.x]]$mild_alarm),
  Strong_alarm  = purrr::map_dbl(names(marker_labels), ~MARKER_CUTS[[.x]]$strong_alarm),
  Source_resp   = purrr::map_chr(names(marker_labels), ~MARKER_CUTS[[.x]]$source_resp),
  Source_alarm  = purrr::map_chr(names(marker_labels), ~MARKER_CUTS[[.x]]$source_alarm)
)

# ── 13b. DIAGNOSTIC GROUP VARIABLE ──────────────────────────────────────────
# mNGS / Blood Culture / Antibiotic combination → 6 clinical groups
# Mirrors the publication table exactly

# Diagnostic: check what values mNGS, BK, Abx columns actually contain
cat("\nDiagnostic group column check:\n")
cat(sprintf("  mNGS_bin unique values: %s\n",
            paste(sort(unique(feature_df$mNGS_bin)), collapse=", ")))
cat(sprintf("  BK_bin   unique values: %s\n",
            paste(sort(unique(feature_df$BK_bin)),   collapse=", ")))
cat(sprintf("  Abx_bin  unique values: %s\n",
            paste(sort(unique(feature_df$Abx_bin)),  collapse=", ")))

# Re-read raw columns directly for diagnostic grouping
# bin01() may not have fired correctly if column names differ
# Use a fresh, tolerant re-extraction from zone_df_raw
get_bin_fresh <- function(df, patterns) {
  col <- find_col(df, patterns)
  if (is.na(col)) return(rep(NA_integer_, nrow(df)))
  bin01(df[[col]])
}

# Actual column names found in your Excel:
# mNGS result  → "erreger"          (pathogen/organism found by NGS)
# BK pos/neg   → "bk_pos_neg"       (blood culture result — already in base_df as BK_bin)
# ABX on/off   → "abx_at_sampling_raw" (antibiotic at sampling, numeric 0/1)

# "erreger" = pathogen name (free text) or empty → positive if non-empty
# Convert: any non-empty, non-NA, non-"neg"/"none"/"kein" string = mNGS positive
erreger_col <- find_col(zone_df_raw, c("^erreger$","erreger"))
if (!is.na(erreger_col)) {
  mNGS_fresh <- as.integer(
    !is.na(zone_df_raw[[erreger_col]]) &
      stringr::str_squish(as.character(zone_df_raw[[erreger_col]])) != "" &
      !stringr::str_detect(
        tolower(stringr::str_squish(as.character(zone_df_raw[[erreger_col]]))),
        "^neg|^no|^kein|^none|^-$|^na$|^0$"
      )
  )
} else {
  mNGS_fresh <- get_bin_fresh(zone_df_raw, c("mngs","ngs","disqver"))
}

# BK pos/neg — already correctly handled by bk_col
BK_fresh <- get_bin_fresh(zone_df_raw, c("^bk","bk_pos","blutkultur","blood.*cult"))

# ABX at sampling raw — numeric 0/1
Abx_fresh <- get_bin_fresh(zone_df_raw, c("^abx_at_sampling_raw$","abx.*sampling.*raw"))

cat(sprintf("  mNGS (erreger): %d positive / %d negative\n",
            sum(mNGS_fresh==1, na.rm=TRUE), sum(mNGS_fresh==0, na.rm=TRUE)))
cat(sprintf("  BK:             %d positive / %d negative\n",
            sum(BK_fresh==1,   na.rm=TRUE), sum(BK_fresh==0,   na.rm=TRUE)))
cat(sprintf("  Abx:            %d on / %d off\n",
            sum(Abx_fresh==1,  na.rm=TRUE), sum(Abx_fresh==0,  na.rm=TRUE)))

cat(sprintf("  mNGS_fresh unique: %s\n", paste(sort(unique(mNGS_fresh)), collapse=", ")))
cat(sprintf("  BK_fresh   unique: %s\n", paste(sort(unique(BK_fresh)),   collapse=", ")))
cat(sprintf("  Abx_fresh  unique: %s\n", paste(sort(unique(Abx_fresh)),  collapse=", ")))

feature_df <- feature_df |>
  dplyr::mutate(
    mNGS_diag = mNGS_fresh,
    BK_diag   = BK_fresh,
    Abx_diag  = Abx_fresh,
    diag_group = dplyr::case_when(
      mNGS_diag==1 & BK_diag==0 & Abx_diag==1 ~ "mNGS+/BC-/Abx+",
      mNGS_diag==1 & BK_diag==1 & Abx_diag==1 ~ "mNGS+/BC+/Abx+",
      mNGS_diag==1 & BK_diag==0 & Abx_diag==0 ~ "mNGS+/BC-/Abx-",
      mNGS_diag==1 & BK_diag==1 & Abx_diag==0 ~ "mNGS+/BC+/Abx-",
      mNGS_diag==0 & BK_diag==0 & Abx_diag==1 ~ "mNGS-/BC-/Abx+",
      mNGS_diag==0 & BK_diag==1 & Abx_diag==1 ~ "mNGS-/BC+/Abx+",
      TRUE ~ "Other"
    ),
    diag_group = factor(diag_group,
                        levels = c("mNGS+/BC-/Abx+","mNGS+/BC+/Abx+",
                                   "mNGS+/BC-/Abx-","mNGS+/BC+/Abx-",
                                   "mNGS-/BC-/Abx+","mNGS-/BC+/Abx+","Other"))
  )

cat("\nDiagnostic group distribution:\n")
print(table(feature_df$diag_group, useNA="ifany"))
cat("\nIf most patients are still 'Other', run this to find correct column names:\n")
cat("  names(zone_df_raw)[grep('ngs|blut|abx|antibiotic|antinfekt', names(zone_df_raw), ignore.case=TRUE)]\n")

# ── 13c. ICU MORTALITY AUC — FULL MODEL ──────────────────────────────────────
# AUC of danger_score as continuous predictor of ICU mortality
# + logistic regression model: danger_score + response_score + normality_score
# This is the key "does the model predict outcome?" analysis

cat("\n── ICU Mortality AUC analysis ──\n")

mort_sub <- feature_df |>
  dplyr::filter(!is.na(ICU_mortality)) |>
  dplyr::mutate(dplyr::across(
    c(danger_score, response_score, normality_score,
      admission_severity, pct_last_delta_score, saps_last_delta_score),
    ~tidyr::replace_na(as.numeric(.x), 0)
  ))

# 1. AUC: danger_score alone vs ICU mortality
roc_danger <- pROC::roc(mort_sub$ICU_mortality, mort_sub$danger_score,
                        quiet=TRUE, levels=c(0,1), direction="<")
auc_danger_ci <- tryCatch(
  as.numeric(pROC::ci.auc(roc_danger, conf.level=0.95, method="delong")),
  error=\(e) c(NA,NA,NA)
)

# 2. AUC: admission_severity (Phase A) vs ICU mortality
roc_admission <- pROC::roc(mort_sub$ICU_mortality, mort_sub$admission_severity,
                           quiet=TRUE, levels=c(0,1), direction="<")
auc_admission_ci <- tryCatch(
  as.numeric(pROC::ci.auc(roc_admission, conf.level=0.95, method="delong")),
  error=\(e) c(NA,NA,NA)
)

# 3. Full logistic model: danger + response + normality + admission_severity
logit_mort <- tryCatch(
  glm(ICU_mortality ~ danger_score + response_score + normality_score +
        admission_severity + pct_last_delta_score + saps_last_delta_score,
      data=mort_sub, family=binomial()),
  error=\(e){message("Logistic model failed: ",e$message); NULL}
)

if (!is.null(logit_mort)) {
  mort_sub$pred_logit <- predict(logit_mort, type="response")
  roc_logit <- pROC::roc(mort_sub$ICU_mortality, mort_sub$pred_logit,
                         quiet=TRUE, levels=c(0,1), direction="<")
  auc_logit_ci <- tryCatch(
    as.numeric(pROC::ci.auc(roc_logit, conf.level=0.95, method="delong")),
    error=\(e) c(NA,NA,NA)
  )
} else {
  roc_logit <- NULL
  auc_logit_ci <- c(NA,NA,NA)
}

# 4. DeLong test: does full model outperform danger_score alone?
delong_test <- if (!is.null(roc_logit)) {
  tryCatch(
    pROC::roc.test(roc_logit, roc_danger, method="delong"),
    error=\(e) NULL
  )
} else NULL

icu_auc_summary <- tibble::tibble(
  Model = c("Admission severity (Phase A)",
            "Danger score (Phase C)",
            "Full logistic model (Phase A+B+C)"),
  AUC   = c(round(as.numeric(pROC::auc(roc_admission)),3),
            round(as.numeric(pROC::auc(roc_danger)),3),
            if(!is.null(roc_logit)) round(as.numeric(pROC::auc(roc_logit)),3) else NA_real_),
  AUC_lower = c(round(auc_admission_ci[1],3),
                round(auc_danger_ci[1],3),
                round(auc_logit_ci[1],3)),
  AUC_upper = c(round(auc_admission_ci[3],3),
                round(auc_danger_ci[3],3),
                round(auc_logit_ci[3],3)),
  DeLong_p  = c(NA_real_, NA_real_,
                if(!is.null(delong_test)) round(delong_test$p.value,4) else NA_real_)
)

cat("\nICU mortality AUC summary:\n")
print(icu_auc_summary)

# Logistic model summary for publication
if (!is.null(logit_mort)) {
  logit_tbl <- broom::tidy(logit_mort, exponentiate=TRUE, conf.int=TRUE) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x,3))) |>
    dplyr::rename(OR=estimate, CI_lower=conf.low, CI_upper=conf.high)
  cat("\nLogistic model (OR for ICU mortality):\n")
  print(logit_tbl)
} else {
  logit_tbl <- tibble::tibble(term="(model failed)")
}

# ── 13d. HOSPITAL MORTALITY + STOP vs REST BINARY ANALYSIS ──────────────────
# Hospital mortality = Died in hospital (broader endpoint, more events)
# Stop vs Rest = key safety question: is stopping ABX safe?

cat("\n── Hospital mortality + Stop vs Rest analysis ──\n")

# Add hospital_mortality to feature_df
feature_df <- feature_df |>
  dplyr::mutate(
    hospital_mortality = as.integer(Died == 1L)
  )

cat(sprintf("Hospital mortality: %d / %d (%.1f%%)\n",
            sum(feature_df$hospital_mortality==1L, na.rm=TRUE),
            sum(!is.na(feature_df$hospital_mortality)),
            mean(feature_df$hospital_mortality, na.rm=TRUE)*100))

# ── Hospital mortality by decision class ──────────────────────────────────────
hosp_mort_df <- feature_df |>
  dplyr::filter(!is.na(AB_change), !is.na(hospital_mortality)) |>
  dplyr::group_by(AB_change) |>
  dplyr::summarise(
    n_total       = dplyr::n(),
    n_hosp_death  = sum(hospital_mortality==1L),
    pct_mort      = n_hosp_death / n_total,
    se            = sqrt(pct_mort*(1-pct_mort)/n_total),
    ci_lo         = pmax(pct_mort - 1.96*se, 0),
    ci_hi         = pmin(pct_mort + 1.96*se, 1),
    .groups       = "drop"
  )
cat("\nHospital mortality by decision class:\n")
print(hosp_mort_df)

# Fisher exact — hospital mortality across 4 decision classes
hosp_fisher4 <- tryCatch({
  cont <- table(
    feature_df$AB_change[!is.na(feature_df$hospital_mortality)],
    feature_df$hospital_mortality[!is.na(feature_df$hospital_mortality)]
  )
  fisher.test(cont, simulate.p.value=TRUE, B=10000)
}, error=\(e) NULL)
hosp_p4 <- if (!is.null(hosp_fisher4)) round(hosp_fisher4$p.value, 4) else NA_real_
cat(sprintf("Fisher exact (4 groups) — hospital mortality: p = %.4f\n", hosp_p4))

# ── STOP vs REST binary comparison ───────────────────────────────────────────
# Primary safety question: is stopping antibiotics safe?
# Stop (n~21) vs Escalate+NoChange+DeEscalate

feature_df <- feature_df |>
  dplyr::mutate(
    stop_vs_rest = dplyr::case_when(
      AB_change == "Stop"                          ~ "Stop",
      AB_change %in% c("Escalate","No change",
                       "De-escalate")             ~ "Continue/Change",
      TRUE                                         ~ NA_character_
    ) |> factor(levels=c("Continue/Change","Stop"))
  )

# 2×2 Fisher: ICU mortality
stop_icu_fisher <- tryCatch({
  cont <- table(
    feature_df$stop_vs_rest[!is.na(feature_df$ICU_mortality)],
    feature_df$ICU_mortality[!is.na(feature_df$ICU_mortality)]
  )
  fisher.test(cont)
}, error=\(e) NULL)

# 2×2 Fisher: hospital mortality
stop_hosp_fisher <- tryCatch({
  cont <- table(
    feature_df$stop_vs_rest[!is.na(feature_df$hospital_mortality)],
    feature_df$hospital_mortality[!is.na(feature_df$hospital_mortality)]
  )
  fisher.test(cont)
}, error=\(e) NULL)

stop_summary <- feature_df |>
  dplyr::filter(!is.na(stop_vs_rest)) |>
  dplyr::group_by(stop_vs_rest) |>
  dplyr::summarise(
    N                  = dplyr::n(),
    ICU_deaths         = sum(ICU_mortality==1L, na.rm=TRUE),
    ICU_mort_pct       = round(ICU_deaths/N*100, 1),
    Hospital_deaths    = sum(hospital_mortality==1L, na.rm=TRUE),
    Hospital_mort_pct  = round(Hospital_deaths/N*100, 1),
    Danger_score_med   = round(median(danger_score, na.rm=TRUE), 1),
    Admission_sev_med  = round(median(admission_severity, na.rm=TRUE), 1),
    .groups            = "drop"
  )

cat("\nStop vs Rest summary:\n")
print(stop_summary)
cat(sprintf("Fisher exact ICU mortality (Stop vs Rest):      p = %.4f\n",
            if(!is.null(stop_icu_fisher))  round(stop_icu_fisher$p.value, 4)  else NA_real_))
cat(sprintf("Fisher exact Hospital mortality (Stop vs Rest): p = %.4f\n",
            if(!is.null(stop_hosp_fisher)) round(stop_hosp_fisher$p.value, 4) else NA_real_))

# ── AUC: full model vs hospital mortality ─────────────────────────────────────
hosp_sub <- feature_df |>
  dplyr::filter(!is.na(hospital_mortality)) |>
  dplyr::mutate(dplyr::across(
    c(danger_score, response_score, normality_score,
      admission_severity, pct_last_delta_score, saps_last_delta_score),
    ~tidyr::replace_na(as.numeric(.x), 0)
  ))

roc_danger_hosp <- pROC::roc(hosp_sub$hospital_mortality, hosp_sub$danger_score,
                             quiet=TRUE, levels=c(0,1), direction="<")
auc_danger_hosp_ci <- tryCatch(
  as.numeric(pROC::ci.auc(roc_danger_hosp, conf.level=0.95, method="delong")),
  error=\(e) c(NA,NA,NA)
)

logit_hosp <- tryCatch(
  glm(hospital_mortality ~ danger_score + response_score + normality_score +
        admission_severity + pct_last_delta_score + saps_last_delta_score,
      data=hosp_sub, family=binomial()),
  error=\(e){message("Hospital logistic failed: ",e$message); NULL}
)

if (!is.null(logit_hosp)) {
  hosp_sub$pred_logit_hosp <- predict(logit_hosp, type="response")
  roc_logit_hosp <- pROC::roc(hosp_sub$hospital_mortality, hosp_sub$pred_logit_hosp,
                              quiet=TRUE, levels=c(0,1), direction="<")
  auc_logit_hosp_ci <- tryCatch(
    as.numeric(pROC::ci.auc(roc_logit_hosp, conf.level=0.95, method="delong")),
    error=\(e) c(NA,NA,NA)
  )
  logit_hosp_tbl <- broom::tidy(logit_hosp, exponentiate=TRUE, conf.int=TRUE) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x,3))) |>
    dplyr::rename(OR=estimate, CI_lower=conf.low, CI_upper=conf.high)
  cat("\nHospital mortality logistic model (OR):\n")
  print(logit_hosp_tbl)
} else {
  roc_logit_hosp   <- NULL
  auc_logit_hosp_ci <- c(NA,NA,NA)
  logit_hosp_tbl   <- tibble::tibble(term="(model failed)")
}

hosp_auc_summary <- tibble::tibble(
  Outcome = "Hospital mortality",
  Model   = c("Danger score (Phase C)",
              "Full logistic model (Phase A+B+C)"),
  AUC     = c(round(as.numeric(pROC::auc(roc_danger_hosp)),3),
              if(!is.null(roc_logit_hosp))
                round(as.numeric(pROC::auc(roc_logit_hosp)),3) else NA_real_),
  AUC_lower = c(round(auc_danger_hosp_ci[1],3),
                round(auc_logit_hosp_ci[1],3)),
  AUC_upper = c(round(auc_danger_hosp_ci[3],3),
                round(auc_logit_hosp_ci[3],3))
)
cat("\nHospital mortality AUC:\n")
print(hosp_auc_summary)

# ── Fig 17 · ICU mortality by diagnostic group and decision class ─────────────
# Panel A: ICU mortality by diagnostic group (mNGS/BC/Abx combination)
#           with Fisher exact p-value
# Panel B: ICU mortality by decision class with Fisher exact p-value
# Panel C: Danger score by diagnostic group
#           (validates that high-mortality groups were truly sicker)

# Fisher exact — ICU mortality across diagnostic groups
diag_icu_fisher <- tryCatch({
  sub  <- feature_df |>
    dplyr::filter(!is.na(ICU_mortality), diag_group != "Other")
  cont <- table(sub$diag_group, sub$ICU_mortality)
  cont <- cont[rowSums(cont) > 0, , drop=FALSE]
  fisher.test(cont, simulate.p.value=TRUE, B=10000)
}, error=\(e) NULL)
diag_icu_p <- if (!is.null(diag_icu_fisher))
  round(diag_icu_fisher$p.value, 4) else NA_real_
cat(sprintf("Fisher exact ICU mortality by diagnostic group: p = %.4f\n", diag_icu_p))

# ICU mortality by diagnostic group
diag_icu_df <- feature_df |>
  dplyr::filter(!is.na(ICU_mortality), diag_group != "Other") |>
  dplyr::group_by(diag_group) |>
  dplyr::summarise(
    N          = dplyr::n(),
    ICU_deaths = sum(ICU_mortality==1L),
    pct_mort   = ICU_deaths / N,
    se         = sqrt(pct_mort*(1-pct_mort)/N),
    ci_lo      = pmax(pct_mort - 1.96*se, 0),
    ci_hi      = pmin(pct_mort + 1.96*se, 1),
    .groups    = "drop"
  )

fig17_a <- ggplot(diag_icu_df,
                  aes(x=diag_group, y=pct_mort, fill=diag_group)) +
  geom_col(width=0.68, alpha=0.88) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi),
                width=0.22, linewidth=0.7) +
  geom_text(aes(label=sprintf("%d/%d\n(%.0f%%)", ICU_deaths, N, pct_mort*100)),
            vjust=-0.4, size=2.8, colour="grey20") +
  annotate("text",
           x    = nlevels(droplevels(diag_icu_df$diag_group)) / 2 + 0.5,
           y    = max(diag_icu_df$ci_hi, na.rm=TRUE) * 1.35,
           label= sprintf("Fisher exact p = %.4f", diag_icu_p),
           size=3.2, colour="grey30", fontface="italic") +
  scale_fill_manual(values=pal_diag, guide="none") +
  scale_y_continuous(labels=percent_format(),
                     limits=c(0, max(diag_icu_df$ci_hi, na.rm=TRUE)*1.55),
                     expand=expansion(mult=c(0, 0.05))) +
  scale_x_discrete(guide=guide_axis(angle=25)) +
  labs(title="A  ICU mortality by diagnostic group (mNGS / Blood culture / Antibiotics)",
       subtitle="Highest mortality: mNGS+/BC+/Abx- = pathogen detected, culture positive, not yet on antibiotics",
       x=NULL, y="ICU mortality rate") +
  theme_pub(base=10) +
  theme(panel.grid.major.x=element_blank())

# Fisher exact — ICU mortality across 4 decision classes
icu_decision_fisher <- tryCatch({
  cont <- table(
    feature_df$AB_change[!is.na(feature_df$ICU_mortality)],
    feature_df$ICU_mortality[!is.na(feature_df$ICU_mortality)]
  )
  fisher.test(cont, simulate.p.value=TRUE, B=10000)
}, error=\(e) NULL)
icu_decision_p <- if (!is.null(icu_decision_fisher))
  round(icu_decision_fisher$p.value, 4) else NA_real_

# ICU mortality by decision class
icu_mort_decision <- feature_df |>
  dplyr::filter(!is.na(AB_change), !is.na(ICU_mortality)) |>
  dplyr::group_by(AB_change) |>
  dplyr::summarise(
    N          = dplyr::n(),
    ICU_deaths = sum(ICU_mortality==1L),
    pct_mort   = ICU_deaths / N,
    se         = sqrt(pct_mort*(1-pct_mort)/N),
    ci_lo      = pmax(pct_mort - 1.96*se, 0),
    ci_hi      = pmin(pct_mort + 1.96*se, 1),
    .groups    = "drop"
  )

fig17_b <- ggplot(icu_mort_decision,
                  aes(x=AB_change, y=pct_mort, fill=AB_change)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=0.22, linewidth=0.7) +
  geom_text(aes(label=sprintf("%d/%d\n(%.0f%%)", ICU_deaths, N, pct_mort*100)),
            vjust=-0.4, size=3.0, colour="grey20") +
  annotate("text", x=2.5,
           y    = max(icu_mort_decision$ci_hi, na.rm=TRUE)*1.5,
           label= sprintf("Fisher p = %.4f", icu_decision_p),
           size=3.0, colour="grey40", fontface="italic") +
  scale_fill_manual(values=pal4, guide="none") +
  scale_y_continuous(labels=percent_format(),
                     limits=c(0, max(icu_mort_decision$ci_hi, na.rm=TRUE)*1.7),
                     expand=expansion(mult=c(0,0.05))) +
  scale_x_discrete(guide=guide_axis(angle=15)) +
  labs(title="B  ICU mortality by ABX decision class",
       subtitle="Stop = 0 deaths; No change = highest in continuing groups",
       x=NULL, y="ICU mortality rate") +
  theme_pub(base=10) +
  theme(panel.grid.major.x=element_blank())

# Panel C: Danger score by diagnostic group
fig17_c <- feature_df |>
  dplyr::filter(!is.na(AB_change), diag_group != "Other") |>
  ggplot(aes(x=diag_group, y=danger_score, fill=diag_group)) +
  geom_violin(alpha=0.30, linewidth=0.5, trim=TRUE) +
  geom_boxplot(width=0.18, fill="white", outlier.shape=NA, linewidth=0.5) +
  scale_fill_manual(values=pal_diag, guide="none") +
  scale_x_discrete(guide=guide_axis(angle=25)) +
  labs(title="C  Danger score by diagnostic group",
       subtitle="Were high-mortality groups truly sicker at the decision point?",
       x=NULL, y="Danger score") +
  theme_pub(base=10)

fig17 <- (fig17_a / (fig17_b | fig17_c)) +
  plot_annotation(
    title    = "Figure 17 \u2014 ICU mortality: diagnostic group and decision class",
    subtitle = paste0(
      "A: ICU mortality by mNGS/BC/Abx group (n=",
      sum(feature_df$diag_group != "Other", na.rm=TRUE),
      ")  |  B: by decision class  |  C: danger score by group"
    ),
    theme = theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=10, colour="grey40")
    )
  )
save_fig(fig17, "Fig17_icu_mortality_diag_decision", 16, 12)

# ── Fig 18 · Kaplan-Meier: 28-day ICU survival by empirical antibiotic adequacy
#
# CLINICAL CONTEXT:
#   All patients: mNGS + BC taken Day 1 → empirical broad-spectrum ABX started Day 1
#   mNGS result returns Day 2-4 → therapy adjusted if needed
#
# ADEQUACY DEFINITION (within mNGS+/Abx+ patients only):
#
#   ADEQUATE              = mNGS+ AND on ABX AND AB_change == "No change"
#                           Empirical antibiotic was correct — mNGS confirmed pathogen
#                           already covered → no adjustment needed
#
#   INADEQUATE /          = mNGS+ AND on ABX AND AB_change == Stop/Escalate/De-escalate
#   NGS-guided              Empirical antibiotic was WRONG or unnecessary →
#                           mNGS guided correction of therapy
#
# HYPOTHESIS:
#   Empirical therapy mismatch (inadequate) → worse 28-day ICU survival
#   even though both groups started on ABX from Day 1
#   Divergence expected from Day 3 when mNGS results return
#
# Time variable : ICU_days (censored at 28 days)
# Event         : ICU_mortality (1 = died in ICU)

cat("\n── Kaplan-Meier: 28-day ICU survival by antibiotic adequacy ──\n")

# SIX-GROUP KM COMPARISON — all diagnostic groups:
#   mNGS+/BC-/Abx+  Pathogen found (NGS), culture neg, on ABX
#   mNGS+/BC+/Abx+  Pathogen found (NGS+culture), on ABX
#   mNGS+/BC-/Abx-  Pathogen found (NGS), culture neg, NOT on ABX
#   mNGS+/BC+/Abx-  Pathogen found (NGS+culture), NOT on ABX — highest risk
#   mNGS-/BC-/Abx+  No pathogen, culture neg, on ABX
#   mNGS-/BC+/Abx+  No pathogen by NGS, culture pos, on ABX
#
# Expected survival order (based on Fig 17):
#   mNGS-/BC-/Abx+ > mNGS-/BC+/Abx+ > mNGS+/BC-/Abx+ >
#   mNGS+/BC+/Abx+ > mNGS+/BC-/Abx- > mNGS+/BC+/Abx-

km_df <- feature_df |>
  dplyr::filter(
    diag_group %in% c("mNGS+/BC-/Abx+","mNGS+/BC+/Abx+",
                      "mNGS+/BC-/Abx-","mNGS+/BC+/Abx-",
                      "mNGS-/BC-/Abx+","mNGS-/BC+/Abx+"),
    !is.na(ICU_days),
    !is.na(ICU_mortality)
  ) |>
  dplyr::mutate(
    # Use diag_group directly as KM grouping variable
    # Order from lowest to highest expected mortality
    Abx_adequacy = factor(diag_group,
                          levels = c(
                            "mNGS-/BC-/Abx+",   # reference — no pathogen, on ABX
                            "mNGS-/BC+/Abx+",   # no pathogen by NGS, culture pos, on ABX
                            "mNGS+/BC-/Abx+",   # NGS positive, culture neg, on ABX
                            "mNGS+/BC+/Abx+",   # NGS + culture positive, on ABX
                            "mNGS+/BC-/Abx-",   # NGS positive, not on ABX
                            "mNGS+/BC+/Abx-"    # NGS + culture positive, NOT on ABX — highest risk
                          )
    ),
    # 28-day censoring
    time_28  = pmin(ICU_days, 28),
    event_28 = dplyr::if_else(ICU_mortality==1L & ICU_days<=28, 1L, 0L)
  ) |>
  dplyr::filter(!is.na(Abx_adequacy))

cat(sprintf("KM analysis: %d patients across 6 diagnostic groups\n", nrow(km_df)))
cat("Group sizes and events:\n")
print(km_df |> dplyr::count(Abx_adequacy, event_28) |> tidyr::pivot_wider(
  names_from=event_28, values_from=n, names_prefix="event_"))

# Kaplan-Meier fit
km_fit <- survival::survfit(
  survival::Surv(time_28, event_28) ~ Abx_adequacy,
  data = km_df
)
print(summary(km_fit)$table)

# Log-rank test
lr_test <- survival::survdiff(
  survival::Surv(time_28, event_28) ~ Abx_adequacy,
  data = km_df
)
lr_p <- round(1 - pchisq(lr_test$chisq, df=length(lr_test$n)-1), 4)
# df = number of groups - 1 = 5 for 6 groups
cat(sprintf("Log-rank test p = %.4f\n", lr_p))

# Cox proportional hazards (adjusted for danger score)
# Cox model — mNGS-/BC-/Abx+ as reference (lowest risk group)
km_df <- km_df |>
  dplyr::mutate(
    Abx_adequacy = relevel(Abx_adequacy, ref="mNGS-/BC-/Abx+")
  )

cox_fit <- tryCatch(
  survival::coxph(
    survival::Surv(time_28, event_28) ~ Abx_adequacy + danger_score +
      admission_severity,
    data = km_df
  ),
  error = \(e) NULL
)

if (!is.null(cox_fit)) {
  cox_tbl <- broom::tidy(cox_fit, exponentiate=TRUE, conf.int=TRUE) |>
    dplyr::mutate(dplyr::across(where(is.numeric), ~round(.x, 3))) |>
    dplyr::rename(HR=estimate, CI_lower=conf.low, CI_upper=conf.high)
  cat("\nCox PH model (HR for 28-day ICU death):\n")
  print(cox_tbl)
} else {
  cox_tbl <- tibble::tibble(term="(model failed)")
}

# ── Plot ──────────────────────────────────────────────────────────────────────
km_colours <- c(
  "mNGS-/BC-/Abx+"  = "#27AE60",   # green   — no pathogen, on ABX (best)
  "mNGS-/BC+/Abx+"  = "#82C341",   # lime    — culture pos but NGS neg, on ABX
  "mNGS+/BC-/Abx+"  = "#F39C12",   # orange  — NGS pos, culture neg, on ABX
  "mNGS+/BC+/Abx+"  = "#E67E22",   # dark orange — both pos, on ABX
  "mNGS+/BC-/Abx-"  = "#2980B9",   # blue    — NGS pos, not on ABX
  "mNGS+/BC+/Abx-"  = "#C0392B"    # red     — both pos, NOT on ABX (worst)
)

# p-value label
p_label <- dplyr::case_when(
  lr_p < 0.001  ~ "Log-rank p < 0.001",
  lr_p < 0.05   ~ sprintf("Log-rank p = %.4f", lr_p),
  TRUE          ~ sprintf("Log-rank p = %.3f (NS)", lr_p)
)

fig18 <- survminer::ggsurvplot(
  km_fit,
  data            = km_df,
  palette         = km_colours,
  conf.int        = TRUE,
  conf.int.alpha  = 0.15,
  risk.table      = TRUE,
  risk.table.col  = "strata",
  risk.table.title= "",
  risk.table.height = 0.28,
  xlab            = "ICU stay (days)",
  ylab            = "Survival probability",
  xlim            = c(0, 28),
  break.time.by   = 7,
  ylim            = c(0.55, 1.0),
  legend.title    = "Therapy group",
  legend.labs     = c("mNGS-/BC-/Abx+","mNGS-/BC+/Abx+",
                      "mNGS+/BC-/Abx+","mNGS+/BC+/Abx+",
                      "mNGS+/BC-/Abx-","mNGS+/BC+/Abx-"),
  title           = "Figure 18 \u2014 Kaplan-Meier: 28-day ICU survival by diagnostic group",
  subtitle        = paste0(
    "Six mNGS/BC/ABX groups  |  n=", nrow(km_df),
    "  |  ", p_label
  ),
  ggtheme         = theme_pub(base=11),
  tables.theme    = survminer::theme_cleantable(),
  fontsize        = 3.5,
  font.main       = c(13, "bold"),
  font.submain    = c(10, "plain", "grey40"),
  font.x          = c(10, "plain"),
  font.y          = c(10, "plain"),
  font.tickslab   = c(9, "plain"),
  surv.median.line= "hv",
  pval            = FALSE   # we annotate manually for control
)

# Add p-value annotation
fig18$plot <- fig18$plot +
  annotate("text",
           x     = 14, y = 0.62,
           label = p_label,
           size  = 3.8, colour = "grey30", fontface = "italic")

# Save — survminer objects must be printed inside an open device
# PDF
pdf_dev <- tryCatch(grDevices::cairo_pdf, error=function(e) grDevices::pdf)
pdf_dev("Figures/Fig18_km_abx_adequacy.pdf", width=9, height=8)
print(fig18)
grDevices::dev.off()

# PNG — open PNG device, print, close
grDevices::png("Figures/Fig18_km_abx_adequacy.png",
               width=9, height=8, units="in", res=300,
               bg="white")
print(fig18)
grDevices::dev.off()
cat("Saved Fig18_km_abx_adequacy\n")

# ── Fig 20 · Incremental AUC + Calibration ───────────────────────────────────
# Panel A: Incremental AUC — Phase A → Phase A+B → Phase A+B+C
# Panel B: Calibration plot — observed vs predicted mortality (quintiles)
#           for both ICU and hospital mortality endpoints

cat("\n── Building Fig 20: Incremental AUC + Calibration ──\n")

# ── Panel A: Incremental AUC ──────────────────────────────────────────────────
# Phase A = admission_severity alone
# Phase A+B = admission_severity + pct_last_delta_score + saps_last_delta_score
# Phase A+B+C = full logistic model

mort_inc <- feature_df |>
  dplyr::filter(!is.na(ICU_mortality)) |>
  dplyr::mutate(dplyr::across(
    c(admission_severity, danger_score, response_score, normality_score,
      pct_last_delta_score, saps_last_delta_score),
    ~tidyr::replace_na(as.numeric(.x), 0)
  ))

# Phase A
fit_A <- glm(ICU_mortality ~ admission_severity,
             data=mort_inc, family=binomial())
roc_A <- pROC::roc(mort_inc$ICU_mortality,
                   predict(fit_A, type="response"),
                   quiet=TRUE, levels=c(0,1), direction="<")
ci_A  <- as.numeric(pROC::ci.auc(roc_A, conf.level=0.95))

# Phase A+B
fit_AB <- glm(ICU_mortality ~ admission_severity +
                pct_last_delta_score + saps_last_delta_score +
                normality_score,
              data=mort_inc, family=binomial())
roc_AB <- pROC::roc(mort_inc$ICU_mortality,
                    predict(fit_AB, type="response"),
                    quiet=TRUE, levels=c(0,1), direction="<")
ci_AB  <- as.numeric(pROC::ci.auc(roc_AB, conf.level=0.95))

# Phase A+B+C (full model — already computed as roc_logit)
ci_ABC <- auc_logit_ci

inc_auc <- tibble::tibble(
  Phase     = factor(c("Phase A","Phase A+B","Phase A+B+C"),
                     levels=c("Phase A","Phase A+B","Phase A+B+C")),
  AUC       = c(as.numeric(pROC::auc(roc_A)),
                as.numeric(pROC::auc(roc_AB)),
                as.numeric(pROC::auc(roc_logit))),
  AUC_lo    = c(ci_A[1],  ci_AB[1],  ci_ABC[1]),
  AUC_hi    = c(ci_A[3],  ci_AB[3],  ci_ABC[3])
)
cat("\nIncremental AUC:\n"); print(inc_auc)

fig20_a <- ggplot(inc_auc, aes(x=Phase, y=AUC)) +
  geom_pointrange(aes(ymin=AUC_lo, ymax=AUC_hi),
                  size=0.9, linewidth=1.2, colour="grey25") +
  geom_text(aes(label=round(AUC,3)),
            vjust=-1.2, size=3.8, fontface="bold", colour="grey20") +
  geom_line(aes(group=1), colour="grey50", linewidth=0.6,
            linetype="dashed") +
  scale_y_continuous(limits=c(0.45, 0.85),
                     breaks=seq(0.5, 0.85, 0.1)) +
  labs(title="A  Incremental ICU mortality discrimination",
       subtitle="Admission severity \u2192 kinetics \u2192 full score model",
       x=NULL, y="AUC") +
  theme_pub(base=11) +
  theme(panel.grid.major.x=element_blank())

# ── Panel B: Calibration plot ─────────────────────────────────────────────────
# Quintiles of predicted risk vs observed mortality
# For both ICU and hospital mortality

calib_df <- mort_inc |>
  dplyr::mutate(
    pred_icu  = predict(logit_mort, newdata=mort_inc, type="response"),
    pred_hosp = if (!is.null(logit_hosp))
      predict(logit_hosp, newdata=mort_inc, type="response")
    else NA_real_,
    quintile_icu  = dplyr::ntile(pred_icu,  5),
    quintile_hosp = dplyr::ntile(pred_hosp, 5)
  )

calib_icu <- calib_df |>
  dplyr::group_by(quintile_icu) |>
  dplyr::summarise(
    pred_mean = mean(pred_icu,       na.rm=TRUE),
    obs_mean  = mean(ICU_mortality,  na.rm=TRUE),
    n         = dplyr::n(),
    se        = sqrt(obs_mean*(1-obs_mean)/n),
    .groups   = "drop"
  ) |>
  dplyr::mutate(Endpoint = "ICU mortality")

calib_hosp <- calib_df |>
  dplyr::filter(!is.na(pred_hosp), !is.na(hospital_mortality)) |>
  dplyr::group_by(quintile_hosp) |>
  dplyr::summarise(
    pred_mean = mean(pred_hosp,         na.rm=TRUE),
    obs_mean  = mean(hospital_mortality,na.rm=TRUE),
    n         = dplyr::n(),
    se        = sqrt(obs_mean*(1-obs_mean)/n),
    .groups   = "drop"
  ) |>
  dplyr::rename(quintile_icu=quintile_hosp) |>
  dplyr::mutate(Endpoint = "Hospital mortality")

calib_all <- dplyr::bind_rows(calib_icu, calib_hosp)

fig20_b <- ggplot(calib_all,
                  aes(x=pred_mean, y=obs_mean,
                      colour=Endpoint, group=Endpoint)) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              colour="grey55", linewidth=0.5) +
  geom_line(linewidth=0.8) +
  geom_pointrange(aes(ymin=pmax(obs_mean-1.96*se,0),
                      ymax=pmin(obs_mean+1.96*se,1)),
                  size=0.9, linewidth=0.8) +
  scale_colour_manual(
    values=c("ICU mortality"="#2C7BB6","Hospital mortality"="#D6604D"),
    name=NULL) +
  scale_x_continuous(labels=percent_format(), limits=c(0,0.5)) +
  scale_y_continuous(labels=percent_format(), limits=c(0,0.5)) +
  coord_equal() +
  labs(title="B  Calibration: ICU and hospital mortality",
       subtitle="Observed mortality across predicted-risk quintiles",
       x="Predicted risk", y="Observed mortality") +
  theme_pub(base=11) +
  theme(legend.position="bottom")

fig20 <- (fig20_a | fig20_b) +
  plot_annotation(
    title    = "Figure 20 \u2014 Incremental model performance and calibration",
    subtitle = "ICU mortality and hospital mortality endpoints; exploratory internal evaluation",
    theme    = theme(
      plot.title    = element_text(face="bold", size=13),
      plot.subtitle = element_text(size=10, colour="grey40")
    )
  )
save_fig(fig20, "Fig20_incremental_auc_calibration", 13, 6)

# ── Fig 21b · KM 28-day hospital survival by antibiotic adequacy ──────────────
# ABX-positive patients only
# Adequate          = No change (empirical ABX was correct)
# Inadequate/NGS    = Escalate + De-escalate + Stop (therapy corrected by mNGS)
# Event             = hospital death within 28 days
# KEY FINDING: p=0.991 — mNGS correction rescues survival to equal adequate therapy

cat("\n── Building Fig 21b: KM hospital survival by antibiotic adequacy ──\n")

km_abx_hosp_df <- feature_df |>
  dplyr::filter(
    Abx_diag == 1,
    !is.na(AB_change),
    !is.na(Hosp_days),
    !is.na(hospital_mortality),
    Hosp_days > 0
  ) |>
  dplyr::mutate(
    Antibiotic_Adequacy = dplyr::case_when(
      AB_change == "No change"                                       ~ "Adequate",
      AB_change %in% c("Escalate","De-escalate","Stop")             ~ "Inadequate / NGS-guided",
      TRUE ~ NA_character_
    ) |> factor(levels=c("Adequate","Inadequate / NGS-guided")),
    time_28  = pmin(Hosp_days, 28),
    event_28 = dplyr::if_else(hospital_mortality==1L & Hosp_days<=28, 1L, 0L)
  ) |>
  dplyr::filter(!is.na(Antibiotic_Adequacy))

cat("\nKM hospital mortality by antibiotic adequacy:\n")
print(table(km_abx_hosp_df$Antibiotic_Adequacy, km_abx_hosp_df$event_28))

km_abx_hosp_fit <- survival::survfit(
  survival::Surv(time_28, event_28) ~ Antibiotic_Adequacy,
  data = km_abx_hosp_df
)
lr_abx_hosp <- survival::survdiff(
  survival::Surv(time_28, event_28) ~ Antibiotic_Adequacy,
  data = km_abx_hosp_df
)
lr_abx_hosp_p <- round(
  1 - pchisq(lr_abx_hosp$chisq, df=length(lr_abx_hosp$n)-1), 4
)
p_label_abx_hosp <- dplyr::case_when(
  lr_abx_hosp_p < 0.001 ~ "Log-rank p < 0.001",
  lr_abx_hosp_p < 0.05  ~ sprintf("Log-rank p = %.4f", lr_abx_hosp_p),
  TRUE                   ~ sprintf("Log-rank p = %.3f", lr_abx_hosp_p)
)
cat(sprintf("Log-rank p = %.4f\n", lr_abx_hosp_p))

km_abx_colours <- c(
  "Adequate"                = "#D6604D",
  "Inadequate / NGS-guided" = "#2C7BB6"
)

fig21b <- survminer::ggsurvplot(
  km_abx_hosp_fit,
  data              = km_abx_hosp_df,
  palette           = km_abx_colours,
  conf.int          = TRUE,
  conf.int.alpha    = 0.18,
  risk.table        = TRUE,
  risk.table.col    = "strata",
  risk.table.height = 0.28,
  risk.table.title  = "",
  xlab              = "Hospital stay (days)",
  ylab              = "Survival probability",
  xlim              = c(0, 28),
  break.time.by     = 7,
  ylim              = c(0.40, 1.00),
  legend.title      = "Therapy group",
  legend.labs       = c("Adequate","Inadequate / NGS-guided"),
  title             = "Figure 21b \u2014 Kaplan\u2013Meier: 28-day hospital survival by antibiotic adequacy",
  subtitle          = paste0(
    "ABX-positive patients only  |  Adequate = no change  |  ",
    "Inadequate/NGS-guided = escalation, de-escalation or stop  |  n=",
    nrow(km_abx_hosp_df), "  |  ", p_label_abx_hosp
  ),
  ggtheme           = theme_pub(base=12),
  tables.theme      = survminer::theme_cleantable(),
  fontsize          = 4.0,
  pval              = FALSE
)

fig21b$plot <- fig21b$plot +
  annotate("text",
           x=15, y=0.48,
           label=paste0(p_label_abx_hosp,
                        "\nNon-inferiority: NGS-guided = Adequate"),
           size=3.8, colour="grey25", fontface="italic",
           lineheight=0.9) +
  theme(legend.position    = "right",
        legend.title       = element_text(size=11),
        legend.text        = element_text(size=10))

# Save
pdf_dev <- tryCatch(grDevices::cairo_pdf, error=function(e) grDevices::pdf)
pdf_dev("Figures/Fig21b_km_hospital_survival_antibiotic_adequacy.pdf",
        width=9.5, height=7.8)
print(fig21b)
grDevices::dev.off()

grDevices::png("Figures/Fig21b_km_hospital_survival_antibiotic_adequacy.png",
               width=9.5, height=7.8, units="in", res=300, bg="white")
print(fig21b)
grDevices::dev.off()
cat("Saved Fig21b_km_hospital_survival_antibiotic_adequacy\n")


# Store cox table for publication statistics
cox_abx_adequacy_tbl <- cox_tbl



# ── Fig 14 · Diagnostic group alluvial → decision → ICU mortality ─────────────
alluv14 <- feature_df |>
  dplyr::filter(!is.na(AB_change), !is.na(ICU_mortality),
                diag_group != "Other") |>
  dplyr::mutate(Outcome = factor(ICU_mortality, levels=c(0,1),
                                 labels=c("Survived","ICU death"))) |>
  dplyr::count(diag_group, AB_change, Outcome)

fig14 <- ggplot(alluv14,
                aes(axis1=diag_group, axis2=AB_change, axis3=Outcome, y=n)) +
  geom_alluvium(aes(fill=diag_group), alpha=0.65, width=0.18) +
  geom_stratum(width=0.25, fill="grey92", colour="grey55") +
  geom_label(stat="stratum", aes(label=after_stat(stratum)),
             size=2.6, label.padding=unit(0.15,"lines")) +
  scale_x_discrete(limits=c("Diagnostic group","Decision","Outcome"),
                   expand=c(0.08,0.05)) +
  scale_fill_manual(values=pal_diag, guide="none") +
  labs(
    title    = "Figure 14 — Diagnostic group → ABX decision → ICU outcome",
    subtitle = "mNGS / Blood culture / Antibiotic status combination drives decision pathway",
    y = "Number of patients"
  ) +
  theme_pub() +
  theme(axis.text.x=element_text(face="bold", size=11))
save_fig(fig14, "Fig14_diag_alluvial", 13, 8)

# ── Fig 15 · Decision by diagnostic group (publication table visualised) ───────
diag_decision <- feature_df |>
  dplyr::filter(!is.na(AB_change), diag_group != "Other") |>
  dplyr::group_by(diag_group) |>
  dplyr::mutate(n_group = dplyr::n()) |>
  dplyr::count(diag_group, AB_change, n_group) |>
  dplyr::mutate(
    pct      = n / n_group,
    adapted  = AB_change != "No change",
    pct_lab  = sprintf("%d\n(%.0f%%)", n, pct*100)
  ) |>
  dplyr::ungroup()

# Adaptation rate summary per group
adapt_rate <- feature_df |>
  dplyr::filter(!is.na(AB_change), diag_group != "Other") |>
  dplyr::group_by(diag_group) |>
  dplyr::summarise(
    n_total    = dplyr::n(),
    n_adapted  = sum(AB_change != "No change"),
    pct_adapted= round(n_adapted/n_total*100, 1),
    .groups    = "drop"
  ) |>
  dplyr::mutate(lab = sprintf("Total adaptations\n%d/%d (%.0f%%)",
                              n_adapted, n_total, pct_adapted))

fig15 <- ggplot(diag_decision,
                aes(x=diag_group, y=pct, fill=AB_change)) +
  geom_col(width=0.72, position="stack") +
  geom_text(aes(label=pct_lab),
            position=position_stack(vjust=0.5),
            size=2.6, colour="white", fontface="bold", lineheight=0.85) +
  # Adaptation rate annotation
  geom_text(data=adapt_rate,
            aes(x=diag_group, y=1.07, label=lab, fill=NULL),
            size=2.4, colour="grey30", lineheight=0.85, inherit.aes=FALSE) +
  scale_fill_manual(values=pal4, name="Decision") +
  scale_y_continuous(labels=percent_format(), limits=c(0,1.25),
                     breaks=seq(0,1,0.25)) +
  scale_x_discrete(guide=guide_axis(angle=25)) +
  labs(
    title    = "Figure 15 — ABX decision by diagnostic group",
    subtitle = "Stacked bars = proportion per decision class; annotation = total adaptation rate",
    x = NULL, y = "Proportion of patients"
  ) +
  theme_pub() +
  theme(legend.position="right", panel.grid.major.x=element_blank())
save_fig(fig15, "Fig15_decision_by_diag_group", 12, 7)

# ── Fig 16 · Danger score + ICU mortality by diagnostic group ─────────────────
fig16_a <- feature_df |>
  dplyr::filter(!is.na(AB_change), diag_group != "Other") |>
  ggplot(aes(x=diag_group, y=danger_score, fill=diag_group)) +
  geom_violin(alpha=0.3, linewidth=0.5, trim=TRUE) +
  geom_boxplot(width=0.18, fill="white", outlier.shape=NA, linewidth=0.5) +
  scale_fill_manual(values=pal_diag, guide="none") +
  scale_x_discrete(guide=guide_axis(angle=35)) +
  labs(title="A  Danger score by diagnostic group",
       x=NULL, y="Danger score") +
  theme_pub(base=9) +
  theme(axis.text.x=element_text(size=7.5))

fig16_b <- feature_df |>
  dplyr::filter(!is.na(ICU_mortality), diag_group != "Other") |>
  dplyr::group_by(diag_group) |>
  dplyr::summarise(
    n_total    = dplyr::n(),
    n_deaths   = sum(ICU_mortality==1L),
    pct_mort   = n_deaths/n_total,
    se         = sqrt(pct_mort*(1-pct_mort)/n_total),
    ci_lo      = pmax(pct_mort-1.96*se, 0),
    ci_hi      = pmin(pct_mort+1.96*se, 1),
    .groups    = "drop"
  ) |>
  ggplot(aes(x=diag_group, y=pct_mort, fill=diag_group)) +
  geom_col(width=0.65, alpha=0.85) +
  geom_errorbar(aes(ymin=ci_lo, ymax=ci_hi), width=0.2, linewidth=0.7) +
  geom_text(aes(label=sprintf("%d/%d\n(%.0f%%)", n_deaths, n_total, pct_mort*100)),
            vjust=-0.5, size=2.8, colour="grey20") +
  scale_fill_manual(values=pal_diag, guide="none") +
  scale_y_continuous(labels=percent_format(), limits=c(0,0.35),
                     expand=expansion(mult=c(0,0.15))) +
  scale_x_discrete(guide=guide_axis(angle=30)) +
  labs(title="B  ICU mortality by diagnostic group",
       x=NULL, y="ICU mortality rate") +
  theme_pub(base=10) +
  theme(panel.grid.major.x=element_blank())

# AUC ROC panel
roc_df_icu <- if (!is.null(roc_logit)) {
  dplyr::bind_rows(
    tibble::tibble(FPR=1-roc_admission$specificities,
                   TPR=roc_admission$sensitivities,
                   Model=sprintf("Admission severity\nAUC=%.3f [%.3f-%.3f]",
                                 as.numeric(pROC::auc(roc_admission)),
                                 auc_admission_ci[1], auc_admission_ci[3])),
    tibble::tibble(FPR=1-roc_danger$specificities,
                   TPR=roc_danger$sensitivities,
                   Model=sprintf("Danger score\nAUC=%.3f [%.3f-%.3f]",
                                 as.numeric(pROC::auc(roc_danger)),
                                 auc_danger_ci[1], auc_danger_ci[3])),
    tibble::tibble(FPR=1-roc_logit$specificities,
                   TPR=roc_logit$sensitivities,
                   Model=sprintf("Full model\nAUC=%.3f [%.3f-%.3f]",
                                 as.numeric(pROC::auc(roc_logit)),
                                 auc_logit_ci[1], auc_logit_ci[3]))
  )
} else {
  dplyr::bind_rows(
    tibble::tibble(FPR=1-roc_admission$specificities,
                   TPR=roc_admission$sensitivities,
                   Model=sprintf("Admission severity AUC=%.3f",
                                 as.numeric(pROC::auc(roc_admission)))),
    tibble::tibble(FPR=1-roc_danger$specificities,
                   TPR=roc_danger$sensitivities,
                   Model=sprintf("Danger score AUC=%.3f",
                                 as.numeric(pROC::auc(roc_danger))))
  )
}

fig16_c <- ggplot(roc_df_icu, aes(FPR, TPR, colour=Model)) +
  geom_line(linewidth=1.1) +
  geom_abline(slope=1, intercept=0, linetype="dashed",
              colour="grey65", linewidth=0.4) +
  scale_colour_brewer(palette="Set1", name=NULL) +
  scale_x_continuous(labels=percent_format(), limits=c(0,1)) +
  scale_y_continuous(labels=percent_format(), limits=c(0,1)) +
  coord_equal() +
  labs(title="C  ROC: model scores vs ICU mortality",
       subtitle="Phase A vs Phase C vs combined logistic model",
       x="1 - Specificity", y="Sensitivity") +
  theme_pub(base=10) +
  theme(legend.position="right",
        legend.text=element_text(size=8))

fig16 <- (fig16_a | fig16_b) / fig16_c +
  plot_annotation(
    title    = "Figure 16 — Diagnostic group: severity, ICU mortality, and model AUC",
    subtitle = "A: danger score distribution  |  B: ICU mortality rate  |  C: ROC vs ICU mortality",
    theme    = theme(plot.title   = element_text(face="bold", size=13),
                     plot.subtitle= element_text(size=10, colour="grey40"))
  )
save_fig(fig16, "Fig16_diag_group_mortality_auc", 16, 11)

# ── PUBLICATION STATISTICS EXCEL ─────────────────────────────────────────────
# A separate clean Excel file with all tables ready for publication
# Table 1: Patient characteristics by decision class
# Table 2: Diagnostic group decision breakdown (mirrors your publication table)
# Table 3: ICU mortality analysis
# Table 4: Model performance (LOOCV + AUC)
# Table 5: Marker-specific thresholds with sources
# Table 6: Logistic regression (OR for ICU mortality)

# Table 1 — Patient characteristics by decision class
tbl1 <- feature_df |>
  dplyr::filter(!is.na(AB_change)) |>
  dplyr::group_by(AB_change) |>
  dplyr::summarise(
    N                  = dplyr::n(),
    Age_median         = round(median(Age, na.rm=TRUE), 1),
    Age_IQR_lo         = round(quantile(Age, 0.25, na.rm=TRUE), 1),
    Age_IQR_hi         = round(quantile(Age, 0.75, na.rm=TRUE), 1),
    ICU_days_median    = round(median(ICU_days, na.rm=TRUE), 1),
    ICU_days_IQR_lo    = round(quantile(ICU_days, 0.25, na.rm=TRUE), 1),
    ICU_days_IQR_hi    = round(quantile(ICU_days, 0.75, na.rm=TRUE), 1),
    Hosp_days_median   = round(median(Hosp_days, na.rm=TRUE), 1),
    mNGS_pos_n         = sum(mNGS_bin==1, na.rm=TRUE),
    mNGS_pos_pct       = round(sum(mNGS_bin==1, na.rm=TRUE)/dplyr::n()*100, 1),
    BK_pos_n           = sum(BK_bin==1, na.rm=TRUE),
    BK_pos_pct         = round(sum(BK_bin==1, na.rm=TRUE)/dplyr::n()*100, 1),
    Abx_on_n           = sum(Abx_bin==1, na.rm=TRUE),
    Abx_on_pct         = round(sum(Abx_bin==1, na.rm=TRUE)/dplyr::n()*100, 1),
    ICU_mortality_n    = sum(ICU_mortality==1L, na.rm=TRUE),
    ICU_mortality_pct  = round(sum(ICU_mortality==1L, na.rm=TRUE)/dplyr::n()*100, 1),
    Admission_sev_med  = round(median(admission_severity, na.rm=TRUE), 2),
    Danger_score_med   = round(median(danger_score, na.rm=TRUE), 2),
    .groups = "drop"
  )

# Fisher exact p-values for ICU mortality across groups
if (sum(!is.na(feature_df$ICU_mortality)) > 0) {
  mort_fisher <- tryCatch({
    cont_tbl <- table(
      feature_df$AB_change[!is.na(feature_df$ICU_mortality)],
      feature_df$ICU_mortality[!is.na(feature_df$ICU_mortality)]
    )
    fisher.test(cont_tbl, simulate.p.value=TRUE, B=10000)
  }, error=\(e) NULL)
  mort_p <- if (!is.null(mort_fisher)) round(mort_fisher$p.value, 4) else NA_real_
} else { mort_p <- NA_real_ }

cat(sprintf("\nFisher exact test — ICU mortality by decision: p = %.4f\n", mort_p))

# Table 2 — Diagnostic group breakdown (publication table)
tbl2 <- feature_df |>
  dplyr::filter(!is.na(AB_change), diag_group != "Other") |>
  dplyr::group_by(diag_group) |>
  dplyr::summarise(
    N_analyzed      = dplyr::n(),
    No_change_n     = sum(AB_change=="No change"),
    No_change_pct   = round(No_change_n/N_analyzed*100, 1),
    Escalation_n    = sum(AB_change=="Escalate"),
    Escalation_pct  = round(Escalation_n/N_analyzed*100, 1),
    DeEscalation_n  = sum(AB_change=="De-escalate"),
    DeEscalation_pct= round(DeEscalation_n/N_analyzed*100, 1),
    Stop_n          = sum(AB_change=="Stop"),
    Stop_pct        = round(Stop_n/N_analyzed*100, 1),
    Total_adapted_n = sum(AB_change!="No change"),
    Total_adapted_pct= round(Total_adapted_n/N_analyzed*100, 1),
    ICU_mortality_n = sum(ICU_mortality==1L, na.rm=TRUE),
    ICU_mortality_pct= round(ICU_mortality_n/N_analyzed*100, 1),
    .groups = "drop"
  )

# Table 3 — ICU mortality analysis
tbl3 <- dplyr::bind_rows(
  mort_df |>
    dplyr::transmute(
      Group         = as.character(AB_change),
      Grouping      = "Decision class",
      N             = n_total,
      ICU_deaths    = n_icu_death,
      Mortality_pct = round(pct_mort*100, 1),
      CI_95_lo      = round(ci_lo*100, 1),
      CI_95_hi      = round(ci_hi*100, 1)
    ),
  feature_df |>
    dplyr::filter(!is.na(ICU_mortality), diag_group != "Other") |>
    dplyr::group_by(diag_group) |>
    dplyr::summarise(
      N=dplyr::n(), ICU_deaths=sum(ICU_mortality==1L),
      pct=ICU_deaths/N,
      se=sqrt(pct*(1-pct)/N),
      .groups="drop"
    ) |>
    dplyr::transmute(
      Group         = as.character(diag_group),
      Grouping      = "Diagnostic group",
      N, ICU_deaths,
      Mortality_pct = round(pct*100, 1),
      CI_95_lo      = round(pmax(pct-1.96*se,0)*100, 1),
      CI_95_hi      = round(pmin(pct+1.96*se,1)*100, 1)
    )
) |>
  dplyr::mutate(Fisher_p_overall = c(rep(mort_p, nrow(mort_df)),
                                     rep(NA_real_, nrow(tbl2))))

# Table 4 — Model performance
tbl4 <- dplyr::bind_rows(
  icu_auc_summary |>
    dplyr::mutate(Analysis="ICU mortality prediction"),
  auc_tbl |>
    dplyr::transmute(
      Model   = paste0("LOOCV: ", Class, " (", Model, ")"),
      AUC, AUC_lower=NA_real_, AUC_upper=NA_real_, DeLong_p=NA_real_,
      Analysis="ABX decision classification"
    )
) |>
  dplyr::mutate(
    LOOCV_accuracy = dplyr::if_else(
      Analysis=="ABX decision classification" &
        stringr::str_detect(Model,"LOOCV Multinomial"),
      round(loocv_acc*100, 1), NA_real_),
    Kappa = dplyr::if_else(
      Analysis=="ABX decision classification" &
        stringr::str_detect(Model,"LOOCV Multinomial"),
      round(loocv_cm$overall["Kappa"], 3), NA_real_)
  )

# Table 5 — Thresholds
tbl5 <- final_thresh_tbl |>
  dplyr::left_join(
    threshold_validation |>
      dplyr::select(marker_key, n_patients, n_icu_deaths, AUC, AUC_lower, AUC_upper,
                    optimal_cut, sensitivity, specificity, updated),
    by=c("Marker"="marker_key")
  )

# Table 6 — Logistic regression
tbl6 <- logit_tbl

# Table 7 — Hospital mortality + Stop vs Rest
tbl7 <- dplyr::bind_rows(
  hosp_mort_df |>
    dplyr::transmute(
      Group           = as.character(AB_change),
      Grouping        = "Decision class",
      N               = n_total,
      Hosp_deaths     = n_hosp_death,
      Mortality_pct   = round(pct_mort*100, 1),
      CI_95_lo        = round(ci_lo*100, 1),
      CI_95_hi        = round(ci_hi*100, 1),
      Fisher_p_4groups= hosp_p4
    ),
  stop_summary |>
    dplyr::transmute(
      Group           = as.character(stop_vs_rest),
      Grouping        = "Stop vs Rest (binary)",
      N,
      Hosp_deaths     = Hospital_deaths,
      Mortality_pct   = Hospital_mort_pct,
      CI_95_lo        = NA_real_,
      CI_95_hi        = NA_real_,
      Fisher_p_4groups= if (!is.null(stop_hosp_fisher))
        round(stop_hosp_fisher$p.value, 4) else NA_real_
    )
)
tbl8 <- logit_hosp_tbl

# Write publication statistics Excel
writexl::write_xlsx(
  list(
    "Table1_Patient_Characteristics" = tbl1,
    "Table2_Diagnostic_Groups"       = tbl2,
    "Table3_ICU_Mortality"           = tbl3,
    "Table4_Model_Performance"       = tbl4,
    "Table5_Marker_Thresholds"       = tbl5,
    "Table6_Logistic_ICU"            = tbl6,
    "Table7_Hospital_Mortality"      = tbl7,
    "Table8_Logistic_Hospital"       = tbl8,
    "Table9_Cox_ABX_Adequacy"        = if(exists("cox_abx_adequacy_tbl"))
      cox_abx_adequacy_tbl else
        tibble::tibble(note="Run Fig18 first")
  ),
  path = "Output/Publication_Statistics.xlsx"
)
cat("\nSaved: Output/Publication_Statistics.xlsx\n")

# ── 14. EXCEL OUTPUT ──────────────────────────────────────────────────────────
# ── 14. EXCEL OUTPUT ──────────────────────────────────────────────────────────

patient_out <- feature_df |>
  dplyr::select(
    JOB_ID, AB_change, Fuzzy_class,
    admission_severity,
    dplyr::ends_with("_t1_sev"), dplyr::ends_with("_t1_zone"),
    n_markers_responding, n_markers_normal,
    dplyr::ends_with("_last_zone"), dplyr::ends_with("_last_sev"),
    dplyr::ends_with("_last_delta"), dplyr::ends_with("_last_delta_score"),
    dplyr::ends_with("_trajectory"), dplyr::ends_with("_response_flag"),
    danger_score, response_score, normality_score,
    hard_alarm_escalate, hard_alarm_stop, deesc_allowed, stop_allowed,
    score_Escalate, score_DeEscalate, score_NoChange, score_Stop,
    ICU_mortality, ICU_days, Hosp_days, Died,
    Age, Infect_type, infect_complexity, BK_bin, mNGS_bin, Abx_bin,
    diag_group
  ) |>
  dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~round(.x,4))) |>
  dplyr::left_join(
    loocv_res |> dplyr::select(JOB_ID,LOOCV_class,dplyr::starts_with("prob_")),
    by="JOB_ID"
  )

writexl::write_xlsx(
  list(
    Patient_Results      = patient_out,
    Threshold_Validation = threshold_validation,
    Final_Thresholds     = final_thresh_tbl,
    LOOCV_AUC            = auc_tbl,
    LOOCV_Confusion      = cm_df,
    ICU_Mortality        = mort_df,
    Variable_Importance  = importance_df,
    Parameters           = tibble::tibble(Parameter=names(PARAM),Value=unlist(PARAM))
  ),
  path = "Output/Fuzzy_ABX_v7_Results.xlsx"
)

write.csv(patient_out,          "Output/Fuzzy_ABX_v7_patients.csv",    row.names=FALSE)
write.csv(threshold_validation, "Output/Fuzzy_ABX_v7_thresholds.csv",  row.names=FALSE)
write.csv(auc_tbl,              "Output/Fuzzy_ABX_v7_auc.csv",         row.names=FALSE)

# ── 15. CONSOLE SUMMARY ───────────────────────────────────────────────────────
cat("\n",strrep("=",65),"\n")
cat("NGS FUZZY ABX STEWARDSHIP v7 — COMPLETE\n")
cat(strrep("=",65),"\n")
cat(sprintf("Patients       : %d (after deduplication)\n", nrow(feature_df)))
cat(sprintf("LOOCV accuracy : %.1f%%  |  Kappa: %.3f\n",
            loocv_acc*100, loocv_cm$overall["Kappa"]))
cat(sprintf("ICU mortality  : %d / %d (%.1f%%)\n",
            sum(feature_df$ICU_mortality==1L,na.rm=TRUE),
            sum(!is.na(feature_df$ICU_mortality)),
            mean(feature_df$ICU_mortality,na.rm=TRUE)*100))
cat("\nFinal marker-specific thresholds:\n")
print(final_thresh_tbl |> dplyr::select(Marker_label,Metric,Unit,
                                        Mild_response,Strong_response,
                                        Mild_alarm,Strong_alarm,
                                        Source_resp,Source_alarm))
cat("\nLOOCV AUC:\n"); print(auc_tbl)
cat("\nFigures saved to ./Figures/\n")
cat("  Fig1   Clinical pathway A\u2192B\u2192C\u2192D\n")
cat("  Fig2   Marker-specific thresholds (membership curves + kinetic zones)\n")
cat("  Fig2b  Threshold table (PDF)\n")
cat("  Fig3   Phase A: Admission severity by decision class\n")
cat("  Fig4   Phase B: Kinetics trajectory heatmap\n")
cat("  Fig5   Delta% at last timepoint vs ICU mortality per marker\n")
cat("  Fig6   ROC: individual marker kinetics vs ICU mortality\n")
cat("  Fig7   Phase C: Danger vs response scatter\n")
cat("  Fig8   Phase C: One-vs-rest ROC curves\n")
cat("  Fig9   LOOCV confusion matrix\n")
cat("  Fig10  Variable importance\n")
cat("  Fig11  Phase D: ICU mortality by decision class\n")
cat("  Fig12  Danger score vs ICU mortality\n")
cat("  Fig13  Alluvial: workbook \u2192 rule engine \u2192 outcome\n")
cat("  Fig14  Diagnostic group \u2192 decision \u2192 ICU outcome (alluvial)\n")
cat("  Fig15  Decision by diagnostic group (publication table)\n")
cat("  Fig16  Diagnostic group: severity, ICU mortality, model AUC\n")
cat("  Fig17  ICU mortality: diagnostic group and decision class\n")
cat("  Fig18  Kaplan-Meier: 28-day ICU survival by diagnostic group\n")
cat("  Fig20  Incremental AUC + calibration plot\n")
cat("  Fig21b Kaplan-Meier: 28-day hospital survival by antibiotic adequacy\n")
cat("\nOutput saved to ./Output/\n")
cat("  Fuzzy_ABX_v7_Results.xlsx\n")
cat("  Publication_Statistics.xlsx  <- ready for paper\n")
cat(strrep("=",65),"\n")


##################################################################################################################



# ================================================================
# Table 2 — Antibiotic decisions by diagnostic group
# Export to Excel (publication-ready)
#
# REQUIRES: feature_df already exists in your environment
#           (i.e. run NGS_Antibiotics_Stewardship_V7.R first, or
#            source() it before running this script)
#
# OUTPUT:   Output/Table2_DiagGroup.xlsx
# ================================================================

library(dplyr)
library(writexl)
library(openxlsx)   # for formatting; install if needed:
# install.packages("openxlsx")

dir.create("Output", showWarnings = FALSE)

# ── 1. BUILD THE SUMMARY TABLE ────────────────────────────────────────────────

# Row-level data (same filter as tbl2 in v7)
tbl2_data <- feature_df |>
  dplyr::filter(!is.na(AB_change), diag_group != "Other") |>
  dplyr::group_by(diag_group) |>
  dplyr::summarise(
    N_analyzed       = dplyr::n(),
    No_change_n      = sum(AB_change == "No change"),
    No_change_pct    = round(No_change_n  / N_analyzed * 100, 1),
    Escalation_n     = sum(AB_change == "Escalate"),
    Escalation_pct   = round(Escalation_n / N_analyzed * 100, 1),
    DeEscalation_n   = sum(AB_change == "De-escalate"),
    DeEscalation_pct = round(DeEscalation_n / N_analyzed * 100, 1),
    Stop_n           = sum(AB_change == "Stop"),
    Stop_pct         = round(Stop_n   / N_analyzed * 100, 1),
    Total_adapted_n  = sum(AB_change != "No change"),
    Total_adapted_pct= round(Total_adapted_n / N_analyzed * 100, 1),
    .groups = "drop"
  )

# Total row
total_row <- feature_df |>
  dplyr::filter(!is.na(AB_change), diag_group != "Other") |>
  dplyr::summarise(
    diag_group        = factor("Total"),
    N_analyzed        = dplyr::n(),
    No_change_n       = sum(AB_change == "No change"),
    No_change_pct     = round(No_change_n  / N_analyzed * 100, 1),
    Escalation_n      = sum(AB_change == "Escalate"),
    Escalation_pct    = round(Escalation_n / N_analyzed * 100, 1),
    DeEscalation_n    = sum(AB_change == "De-escalate"),
    DeEscalation_pct  = round(DeEscalation_n / N_analyzed * 100, 1),
    Stop_n            = sum(AB_change == "Stop"),
    Stop_pct          = round(Stop_n   / N_analyzed * 100, 1),
    Total_adapted_n   = sum(AB_change != "No change"),
    Total_adapted_pct = round(Total_adapted_n / N_analyzed * 100, 1)
  )

tbl2_full <- dplyr::bind_rows(tbl2_data, total_row)

# ── 2. FORMAT FOR DISPLAY (match the screenshot layout) ──────────────────────
# Combine n and % into "n (x.x%)" strings, just like the published table

fmt_cell <- function(n, p) {
  ifelse(n == 0 & p == 0, "0",
         paste0(n, " (", formatC(p, format = "f", digits = 1), "%)"))
}

# Parse mNGS / BC / Abx signs from the diag_group factor label
parse_signs <- function(label) {
  label <- as.character(label)
  if (label == "Total") return(data.frame(mNGS = "Total", BC = "", Abx = ""))
  parts <- strsplit(label, "/")[[1]]           # e.g. c("mNGS+","BC-","Abx+")
  signs <- sub(".*([+\\-])", "\\1", parts)     # extract last char: +/-
  data.frame(
    mNGS = signs[1],
    BC   = signs[2],
    Abx  = signs[3],
    stringsAsFactors = FALSE
  )
}

sign_df <- do.call(rbind, lapply(tbl2_full$diag_group, parse_signs))

pub_table <- data.frame(
  mNGS                    = sign_df$mNGS,
  BC                      = sign_df$BC,
  Abx                     = sign_df$Abx,
  `Total/excluded/analyzed (n)` = paste0(
    c(
      sapply(as.character(tbl2_data$diag_group), function(g) {
        r <- tbl2_data[as.character(tbl2_data$diag_group) == g, ]
        # total/excl not stored in tbl2 — re-derive from feature_df
        tot  <- nrow(feature_df[as.character(feature_df$diag_group) == g, ])
        excl <- tot - r$N_analyzed
        paste0(tot, "/", excl, "/", r$N_analyzed)
      }),
      paste0(
        nrow(feature_df[feature_df$diag_group != "Other", ]),
        "/",
        nrow(feature_df[feature_df$diag_group != "Other", ]) - total_row$N_analyzed,
        "/",
        total_row$N_analyzed
      )
    )
  ),
  `No change`             = fmt_cell(tbl2_full$No_change_n,      tbl2_full$No_change_pct),
  Escalation              = fmt_cell(tbl2_full$Escalation_n,     tbl2_full$Escalation_pct),
  `De-escalation`         = fmt_cell(tbl2_full$DeEscalation_n,   tbl2_full$DeEscalation_pct),
  Stop                    = fmt_cell(tbl2_full$Stop_n,           tbl2_full$Stop_pct),
  `Total adaptations n (%)` = fmt_cell(tbl2_full$Total_adapted_n, tbl2_full$Total_adapted_pct),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ── 3. EXPORT: simple writexl (no formatting) ─────────────────────────────────
writexl::write_xlsx(
  list(
    Table2_DiagGroup     = pub_table,
    Table2_Raw_Numbers   = as.data.frame(tbl2_full)
  ),
  path = "Output/Table2_DiagGroup.xlsx"
)
cat("Saved (basic): Output/Table2_DiagGroup.xlsx\n")

# ── 4. EXPORT: openxlsx (formatted, publication-ready) ────────────────────────

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Table 2")

# ── Styles ──
style_header <- openxlsx::createStyle(
  fontName  = "Calibri", fontSize = 11, fontColour = "white",
  fgFill    = "#2E4057", halign = "center", valign = "center",
  textDecoration = "bold", wrapText = TRUE,
  border    = "TopBottomLeftRight", borderColour = "#FFFFFF"
)
style_group <- openxlsx::createStyle(
  fontName  = "Calibri", fontSize = 10,
  halign    = "center", valign = "center",
  border    = "TopBottomLeftRight", borderColour = "#D0D0D0"
)
style_data <- openxlsx::createStyle(
  fontName  = "Calibri", fontSize = 10,
  halign    = "center", valign = "center",
  border    = "TopBottomLeftRight", borderColour = "#D0D0D0"
)
style_total <- openxlsx::createStyle(
  fontName  = "Calibri", fontSize = 10, textDecoration = "bold",
  halign    = "center", valign = "center", fgFill = "#F0F0F0",
  border    = "TopBottomLeftRight", borderColour = "#A0A0A0"
)
style_adapt <- openxlsx::createStyle(
  fontName  = "Calibri", fontSize = 10, textDecoration = "bold",
  halign    = "center", valign = "center",
  border    = "TopBottomLeftRight", borderColour = "#D0D0D0"
)

nrow_data <- nrow(pub_table)
ncol_data <- ncol(pub_table)

# Write header row
openxlsx::writeData(wb, "Table 2", pub_table, startRow = 2, startCol = 1,
                    colNames = TRUE, rowNames = FALSE)

# Apply header style (row 2)
openxlsx::addStyle(wb, "Table 2", style_header,
                   rows = 2, cols = 1:ncol_data, gridExpand = TRUE)

# Apply data rows (rows 3 to nrow_data+1)
for (r in 3:(nrow_data + 1)) {
  is_total <- r == (nrow_data + 2)
  s <- if (is_total) style_total else style_data
  openxlsx::addStyle(wb, "Table 2", s,
                     rows = r, cols = 1:ncol_data, gridExpand = TRUE)
}

# Bold the "Total adaptations" column
openxlsx::addStyle(wb, "Table 2", style_adapt,
                   rows = 3:(nrow_data + 2), cols = ncol_data, gridExpand = TRUE)

# Total row style (last data row)
openxlsx::addStyle(wb, "Table 2", style_total,
                   rows = nrow_data + 2, cols = 1:ncol_data, gridExpand = TRUE)

# Column widths
openxlsx::setColWidths(wb, "Table 2",
                       cols = 1:ncol_data,
                       widths = c(7, 7, 7, 22, 14, 14, 14, 10, 22))

# Row heights
openxlsx::setRowHeights(wb, "Table 2", rows = 2,             heights = 30)
openxlsx::setRowHeights(wb, "Table 2", rows = 3:(nrow_data+2), heights = 18)

# Freeze header
openxlsx::freezePane(wb, "Table 2", firstActiveRow = 3)

# Title cell above the table
openxlsx::writeData(wb, "Table 2",
                    "Table 2 — Antibiotic stewardship decisions by diagnostic group (mNGS / Blood culture / Antibiotics)",
                    startRow = 1, startCol = 1)
openxlsx::addStyle(wb, "Table 2",
                   openxlsx::createStyle(fontSize = 12, textDecoration = "bold",
                                         fontName = "Calibri"),
                   rows = 1, cols = 1)
openxlsx::mergeCells(wb, "Table 2", cols = 1:ncol_data, rows = 1)

# ── Raw numbers sheet ──
openxlsx::addWorksheet(wb, "Raw numbers")
openxlsx::writeDataTable(wb, "Raw numbers",
                         as.data.frame(tbl2_full),
                         startRow = 1, startCol = 1,
                         tableStyle = "TableStyleMedium9")

openxlsx::saveWorkbook(wb, "Output/Table2_DiagGroup_Formatted.xlsx", overwrite = TRUE)
cat("Saved (formatted): Output/Table2_DiagGroup_Formatted.xlsx\n")
cat("\nTable preview:\n")
print(pub_table, right = FALSE)

