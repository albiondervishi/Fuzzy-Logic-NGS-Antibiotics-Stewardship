# ================================================================
# Publication summary Table 1 — NGS Stewardship  (v7)
#
# CHANGE LOG vs v5:
#   [1] therapy_group is NO LONGER derived from raw "ab_change".
#       It is now derived from feature_df$AB_change, which has
#       already passed through v8's:
#         - Fallnummer (patient-level) deduplication
#         - "Anderung durchgefuhrt" execution gate
#           (indicated-but-not-executed -> "No change")
#       This is the SAME AB_change used to build Table 2, so the
#       two tables now share one source of truth.
#   [2] Cohort restricted to diag_group != "Other", matching the
#       243-patient denominator used in Table 2 (one patient with
#       indeterminate mNGS/BC/Abx classification is excluded, vs
#       244 in the old v5 cohort).
#
# CHANGE LOG vs v6 (this version):
#   [3] ICU_days / Hosp_days / Died / hospital_mortality / icu_mortality_28d
#       are NO LONGER pulled from feature_df. Verified against the raw
#       sheet: 295 unique patients (Fallnummer), Died is 100% consistent
#       across duplicate rows of the same patient, and only 1 missing
#       value total across ICU days / Hosp. days combined. feature_df's
#       versions of these columns produced 69 spurious "Unknown" LOS
#       values and undercounted deaths by ~56, almost certainly due to
#       a row-alignment bug inside v8 around the diag_group/
#       hospital_mortality mutate block (separate fix needed in v8
#       itself - tracked separately). Until that's fixed, raw_df is
#       the trustworthy source for these four outcome fields. Only
#       AB_change and diag_group are taken from feature_df (the
#       execution-gated stewardship logic, which is NOT affected by
#       the row-alignment issue since it's verified against Table 2's
#       totals).
#   [4] Sex label "w" -> "f" for publication consistency.
#   [5] Infection type "negativ" -> "Negative" for publication consistency.
#
# Expected final N: Overall = 243 | Adequate = 119 | Inadequate/NGS-guided = 124
# ================================================================

rm(list = ls())
setwd("~/Desktop/NGS_ZNA/2_NGS_Stewarship_Antibiotic")

library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(gtsummary)
library(flextable)
library(officer)

# ------------------------------------------------
# 0. Load feature_df (single source of truth for AB_change / diag_group)
# ------------------------------------------------
# Option A: run NGS_Fuzzy_ABX_Stewardship_v8_Prina.R first in this session
#           -> feature_df already exists.
# Option B: load the CSV that v8 writes out at the end of its run.

if (!exists("feature_df")) {
  cat("feature_df not found in session - loading from CSV\n")
  feature_df <- read.csv("Output/Fuzzy_ABX_v8_patients.csv", stringsAsFactors = FALSE)
}
cat(sprintf("feature_df rows loaded: %d\n", nrow(feature_df)))

# ------------------------------------------------
# 1. Load raw sheet for demographic / classification columns
#    NOT carried in feature_df (sex, infection-type free text,
#    antibiotic-at-sampling free text, pathogen free text).
#    Join key: JOB_ID == job_id (same row v8 kept after Fallnummer dedup).
# ------------------------------------------------

file <- "NGS_Kinetics_Zones_v2.xlsx"

raw_df <- read_excel(file, sheet = 1) %>%
  clean_names()

# ------------------------------------------------
# Helper functions
# ------------------------------------------------

to_num <- function(x) {
  as.numeric(str_replace_all(as.character(x), ",", "."))
}

# ------------------------------------------------
# Sex
# ------------------------------------------------

raw_df <- raw_df %>%
  mutate(
    sex = factor(
      case_when(
        str_to_lower(as.character(sex)) %in% c("m", "male", "mannlich") ~ "m",
        str_to_lower(as.character(sex)) %in% c("w", "f", "female", "weiblich") ~ "f",
        TRUE ~ NA_character_
      ),
      levels = c("m", "f")
    )
  )

# ------------------------------------------------
# Died -> binary (verified against raw sheet: 295 unique Fallnummer,
# Died 100% consistent across duplicate rows, only 1 missing value
# combined across ICU days / Hosp. days)
# ------------------------------------------------

died_to_binary <- function(x) {
  case_when(
    str_to_lower(as.character(x)) %in% c("1", "yes", "ja", "dead", "death", "died", "true", "1.0") ~ 1L,
    str_to_lower(as.character(x)) %in% c("0", "no", "nein", "alive", "survived", "false", "0.0") ~ 0L,
    TRUE ~ NA_integer_
  )
}

raw_df <- raw_df %>%
  mutate(
    icu_days_raw  = to_num(icu_days),
    hosp_days_raw = to_num(hosp_days),
    died_raw      = died_to_binary(died),

    # Hospital mortality - straight from raw "Died" column
    hospital_mortality_raw = factor(
      died_raw,
      levels = c(0L, 1L),
      labels = c("Survived", "Died")
    ),

    # 28-day ICU mortality (same definition v5 used):
    # died AND never left ICU (ICU days == Hosp. days) -> ICU death
    icu_mortality_28d_raw = factor(
      case_when(
        died_raw == 1L &
          !is.na(icu_days_raw) & !is.na(hosp_days_raw) &
          icu_days_raw == hosp_days_raw                ~ 1L,
        died_raw == 0L                                  ~ 0L,
        died_raw == 1L                                  ~ 0L,  # died after ICU transfer out
        TRUE                                             ~ NA_integer_
      ),
      levels = c(0L, 1L),
      labels = c("Survived", "Died")
    )
  )

# ------------------------------------------------
# Antibiotic at sampling - classify free-text antiinfektiva
# into the 11 publication categories  (unchanged from v5)
# ------------------------------------------------

classify_abx <- function(x) {
  x_lo <- str_to_lower(as.character(x))
  case_when(
    is.na(x) | str_trim(x_lo) == ""
      ~ "No antibiotic therapy",

    str_detect(x_lo, "caspo|caspofungin|anidulafungin|amphotericin|fluconazol|voriconazol") &
      !str_detect(x_lo, "pip|taz|mero|ceftri|ampi|amp|levo|cipro|vanco|linezolid|clinda|penicillin|cefaz|cefuro|ceftaz|cotrim|doxy")
      ~ "Antifungal therapy",

    str_detect(x_lo, "aciclovir|acilovir|acyclovir") &
      !str_detect(x_lo, "pip|taz|mero|ceftri|ampi|amp|levo|cipro|vanco|linezolid|clinda|penicillin|cefaz|cefuro|ceftaz|cotrim|doxy")
      ~ "Antiviral therapy",

    str_detect(x_lo, "meropenem|menopenem|imipenem|ertapenem|ceftazidim.*avibactam|ceftazidim/avibactam")
      ~ "Carbapenem / reserve broad-spectrum therapy",

    str_detect(x_lo, "linezolid|vancomycin|vanco|daptomycin|staphylex|flucloxacillin") &
      !str_detect(x_lo, "meropenem|menopenem|imipenem")
      ~ "Gram-positive reserve therapy",

    str_detect(x_lo, "clinda") &
      !str_detect(x_lo, "meropenem|pip|taz|ceftri|ceftaz|levo|cipro|ampi|amp|vanco|linezolid")
      ~ "Clindamycin / toxin-suppression coverage",

    str_detect(x_lo, "pip.*taz|piperacillin.*tazobactam|pip/taz|piptaz") &
      !str_detect(x_lo, "meropenem|menopenem|imipenem")
      ~ "Broad-spectrum beta-lactam therapy",

    (str_detect(x_lo, "ceftriax|ceftr|cetfr|cetriax") &
       str_detect(x_lo, "metro|metronidazol|metronidaxol|metromidazle|clont|met\\b")) |
    str_detect(x_lo, "cef/clont|cef/met|cetfr.*clont")
      ~ "Cephalosporin / anaerobic coverage",

    str_detect(x_lo, "levo|cipro|ciproflox|roxi|roxithromycin|clarithromycin|erythromycin|doxycyclin|cotrimoxazol|cotrim|azithromycin") &
      !str_detect(x_lo, "meropenem|pip.*taz|piperacillin") &
      !str_detect(x_lo, "ceftriax|ceftr|amp.*sulb|ampi.*sulb|amoclav|cefazolin|cephazolin|cefuroxim")
      ~ "Quinolone / macrolide / atypical coverage",

    str_detect(x_lo, "amp.*sulb|ampi.*sulb|amoclav|ampicillin.*sulbactam|amp/sulb|ampi/sulb") |
    str_detect(x_lo, "^ceftriax|^ceftr|^cetriax|^cetfr") & !str_detect(x_lo, "metro|clont|met\\b|mero") |
    str_detect(x_lo, "cefazolin|cephazolin|cefuroxim|cefurox|cef\\b") & !str_detect(x_lo, "metro|clont") |
    str_detect(x_lo, "^penicillin|^pen g|^pen\\b") & !str_detect(x_lo, "clinda|mero") |
    str_detect(x_lo, "^flucloxacillin$|^staphylex$")
      ~ "Standard beta-lactam therapy",

    TRUE ~ "Other / unclear antimicrobial therapy"
  )
}

raw_df <- raw_df %>%
  mutate(
    abx_at_sampling = factor(
      classify_abx(antiinfektiva),
      levels = c(
        "Antifungal therapy",
        "Antiviral therapy",
        "Broad-spectrum beta-lactam therapy",
        "Carbapenem / reserve broad-spectrum therapy",
        "Cephalosporin / anaerobic coverage",
        "Clindamycin / toxin-suppression coverage",
        "Gram-positive reserve therapy",
        "No antibiotic therapy",
        "Other / unclear antimicrobial therapy",
        "Quinolone / macrolide / atypical coverage",
        "Standard beta-lactam therapy"
      )
    )
  )

other_rows <- raw_df %>%
  filter(abx_at_sampling == "Other / unclear antimicrobial therapy") %>%
  select(job_id, antiinfektiva, abx_at_sampling)
if (nrow(other_rows) > 0) {
  cat("\nRows classified as 'Other' - review:\n")
  print(other_rows)
}

# ------------------------------------------------
# Pathogen family grouping (unchanged from v5)
# ------------------------------------------------

raw_df <- raw_df %>%
  mutate(
    pathogen_family_group = case_when(
      is.na(erreger) | str_trim(as.character(erreger)) == "" ~
        "Unknown / not documented",
      str_detect(str_to_lower(erreger), "negativ|negative|no pathogen|kein erreger") ~
        "No pathogen detected",
      str_detect(str_to_lower(erreger),
                 "escherichia|klebsiella|enterobacter|citrobacter|proteus|serratia|morganella|hafnia") ~
        "Enterobacterales",
      str_detect(str_to_lower(erreger), "enterococcus") ~
        "Enterococcus spp.",
      str_detect(str_to_lower(erreger), "streptococcus") ~
        "Streptococcus spp.",
      str_detect(str_to_lower(erreger), "staphylococcus") ~
        "Staphylococcus spp.",
      str_detect(str_to_lower(erreger), "pseudomonas") ~
        "Pseudomonas spp.",
      str_detect(str_to_lower(erreger),
                 "bacteroides|phocaeicola|prevotella|clostridium|fusobacterium|alistipes|parabacteroides|finegoldia") ~
        "Anaerobes / gut-associated bacteria",
      str_detect(str_to_lower(erreger), "candida|aspergillus|saccharomyces") ~
        "Fungal pathogens",
      str_detect(str_to_lower(erreger),
                 "epstein|cytomegalovirus|herpes|hhv|hsv|varizella|torque teno|polyomavirus|influenza|sars|covid|rhinovirus|adenovirus|norovirus") ~
        "Viral pathogens",
      TRUE ~ "Other bacteria / mixed organisms"
    ),
    pathogen_family_group = factor(pathogen_family_group, levels = c(
      "Anaerobes / gut-associated bacteria",
      "Enterobacterales",
      "Enterococcus spp.",
      "Fungal pathogens",
      "No pathogen detected",
      "Other bacteria / mixed organisms",
      "Pseudomonas spp.",
      "Staphylococcus spp.",
      "Streptococcus spp.",
      "Unknown / not documented",
      "Viral pathogens"
    ))
  )

# ------------------------------------------------
# Infection type label cleanup ("negativ" -> "Negative")
# ------------------------------------------------

raw_df <- raw_df %>%
  mutate(
    infect_type_clean = case_when(
      str_to_lower(str_trim(as.character(infect_type))) == "negativ" ~ "Negative",
      TRUE ~ as.character(infect_type)
    )
  )

# ------------------------------------------------
# 2. JOIN raw demographic/classification/outcome columns onto feature_df
#    feature_df is patient-level (post Fallnummer dedup); join by job_id.
#    AB_change and diag_group come from feature_df (v8's execution-gated
#    stewardship logic - verified against Table 2's totals).
#    Age, sex, ICU/Hosp days, Died, infect_type, abx_at_sampling,
#    pathogen_family_group come from raw_df (verified complete/consistent;
#    feature_df's versions of the outcome columns are NOT used here -
#    see CHANGE LOG [3] above).
# ------------------------------------------------

# NOTE: feature_df's id column is "JOB_ID" (uppercase, from v8's transmute).
#       raw_df's id column after clean_names() is "job_id".
feature_df <- feature_df %>%
  mutate(job_id = as.character(JOB_ID))

raw_df <- raw_df %>%
  mutate(job_id = as.character(job_id))

joined <- feature_df %>%
  select(job_id, AB_change, diag_group) %>%
  left_join(
    raw_df %>%
      select(
        job_id, age, sex, icu_days_raw, hosp_days_raw,
        hospital_mortality_raw, icu_mortality_28d_raw,
        infect_type_clean, abx_at_sampling, pathogen_family_group
      ),
    by = "job_id"
  )

cat(sprintf("\nRows after join: %d (feature_df had %d)\n",
            nrow(joined), nrow(feature_df)))

# Sanity check: flag any unmatched rows (job_id in feature_df not found in raw_df)
n_unmatched <- sum(is.na(joined$hospital_mortality_raw))
if (n_unmatched > 0) {
  cat(sprintf("WARNING: %d rows did not match a raw_df job_id - review join keys.\n", n_unmatched))
}

# ------------------------------------------------
# 3. Restrict to diag_group != "Other" -> matches Table 2's N = 243
# ------------------------------------------------

n_before <- nrow(joined)
joined <- joined %>% filter(diag_group != "Other")
cat(sprintf("Excluded (diag_group == 'Other'): %d\n", n_before - nrow(joined)))

# ------------------------------------------------
# 4. therapy_group derived from AB_change (execution-gated), NOT raw ab_change
#    "No change"                              -> Adequate
#    "Escalate" / "De-escalate" / "Stop"      -> Inadequate / NGS-guided
# ------------------------------------------------

joined <- joined %>%
  mutate(
    therapy_group = factor(
      case_when(
        AB_change == "No change"                                ~ "Adequate",
        AB_change %in% c("Escalate", "De-escalate", "Stop")      ~ "Inadequate / NGS-guided",
        TRUE                                                     ~ NA_character_
      ),
      levels = c("Adequate", "Inadequate / NGS-guided")
    )
  )

cat("\nFinal N by group (should match Table 2 Total row: 119 / 124):\n")
print(table(joined$therapy_group, useNA = "ifany"))
cat("Total N:", sum(!is.na(joined$therapy_group)), "\n")

# ------------------------------------------------
# 5. Build publication dataframe
#    Age / sex / ICU & Hosp days / mortality / infect_type sourced from
#    raw_df (see [3] above); AB_change-derived therapy_group + diag_group
#    filter from feature_df.
# ------------------------------------------------

pub_df <- joined %>%
  transmute(
    job_id    = job_id,                              # unique patient ID - not in table
    therapy_group = therapy_group,
    age      = to_num(age),
    sex      = sex,
    icu_days = icu_days_raw,
    hosp_days = hosp_days_raw,

    icu_mortality_28d  = icu_mortality_28d_raw,
    hospital_mortality = hospital_mortality_raw,

    infect_type            = factor(infect_type_clean),
    abx_at_sampling        = abx_at_sampling,
    pathogen_family_group  = pathogen_family_group
  ) %>%
  filter(!is.na(therapy_group))

# ── job_id uniqueness check ──────────────────────────────
n_unique <- length(unique(pub_df$job_id))
dups     <- pub_df$job_id[duplicated(pub_df$job_id)]
if (length(dups) == 0) {
  cat(sprintf("\nUniqueness check: ALL %d job_id unique - no duplicates.\n", n_unique))
} else {
  cat(sprintf("\nWARNING: %d duplicate job_id(s) detected - review before publication:\n",
              length(dups)))
  print(pub_df %>% filter(job_id %in% dups) %>%
          select(job_id, therapy_group, icu_days, hosp_days, hospital_mortality))
}

cat("\n28-day ICU mortality distribution:\n")
print(table(pub_df$icu_mortality_28d, pub_df$therapy_group, useNA = "ifany"))

# ------------------------------------------------
# p-value formatter
# ------------------------------------------------

fmt_p <- function(x) {
  dplyr::case_when(
    is.na(x)  ~ NA_character_,
    x < 0.001 ~ "<0.001",
    x > 0.9   ~ ">0.9",
    TRUE       ~ formatC(x, digits = 3, format = "f")
  )
}

# ------------------------------------------------
# Table 1
# ------------------------------------------------

table1 <- pub_df %>%
  select(-job_id) %>%                          # ID column - not a table variable
  tbl_summary(
    by = therapy_group,
    type = list(
      age               ~ "continuous",
      icu_days          ~ "continuous",
      hosp_days         ~ "continuous",
      icu_mortality_28d ~ "categorical",
      hospital_mortality ~ "categorical"
    ),
    statistic = list(
      age       ~ "{median} [{p25}, {p75}]",
      icu_days  ~ "{median} [{p25}, {p75}]",
      hosp_days ~ "{median} [{p25}, {p75}]",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      all_continuous()  ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "ifany",
    missing_text = "Unknown",
    label = list(
      age                ~ "Age, years",
      sex                ~ "Sex",
      icu_days           ~ "ICU length of stay, days",
      hosp_days          ~ "Hospital length of stay, days",
      icu_mortality_28d  ~ "28-day ICU mortality",
      hospital_mortality ~ "Hospital mortality",
      infect_type            ~ "Infection type",
      abx_at_sampling        ~ "Antibiotic therapy at sampling",
      pathogen_family_group  ~ "Detected pathogen group"
    )
  ) %>%
  add_overall(last = FALSE) %>%
  add_p(
    test = list(
      all_continuous()  ~ "wilcox.test",
      all_categorical() ~ "chisq.test"
    ),
    pvalue_fun = fmt_p
  ) %>%
  bold_labels() %>%
  modify_header(
    label   ~ "**Characteristic**",
    stat_0  ~ "**Overall**\nN = {N}",
    stat_1  ~ "**Adequate**\nN = {n}",
    stat_2  ~ "**Inadequate / NGS-guided**\nN = {n}",
    p.value ~ "**p-value**"
  ) %>%
  modify_footnote(
    all_stat_cols() ~ "Median [Q1, Q3]; n (%)",
    p.value         ~ "Wilcoxon rank sum test; Pearson's Chi-squared test"
  )

print(table1)

# ------------------------------------------------
# Export -> Word
# ------------------------------------------------

doc <- read_docx() %>%
  body_add_par(
    paste0("Table 1. Demographic, clinical, microbiological and antibiotic ",
           "stewardship characteristics stratified by therapy adequacy."),
    style = "heading 1"
  ) %>%
  body_add_flextable(
    as_flex_table(table1) %>%
      fontsize(size = 9, part = "all") %>%
      font(fontname = "Calibri", part = "all") %>%
      set_table_properties(width = 1, layout = "autofit")
  )

print(doc, target = "Table1_Publication_Summary_v7.docx")
cat("\nDone. File written: Table1_Publication_Summary_v7.docx\n")
