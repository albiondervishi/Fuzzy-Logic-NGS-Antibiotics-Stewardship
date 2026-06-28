# ================================================================
# Publication summary Table 1 — NGS Stewardship
# Final N: Overall=244 | Adequate=131 | Inadequate/NGS-guided=113
# Transfer patient (n=1) excluded per clinical review
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

file <- "NGS_Kinetics_Zones_v2.xlsx"

df <- read_excel(file, sheet = 1) %>%
  clean_names()

# ------------------------------------------------
# Helper functions
# ------------------------------------------------

to_num <- function(x) {
  as.numeric(str_replace_all(as.character(x), ",", "."))
}

mortality_to_binary <- function(x) {
  case_when(
    str_to_lower(as.character(x)) %in%
      c("1", "yes", "ja", "dead", "death", "died", "true", "1.0") ~ 1L,
    str_to_lower(as.character(x)) %in%
      c("0", "no", "nein", "alive", "survived", "false", "0.0") ~ 0L,
    TRUE ~ NA_integer_
  )
}

# ------------------------------------------------
# Therapy group — exact match
# "None" → Adequate | "De-escalate","Escalate","Stop" → NGS-guided
# "x" → NA (transferred patient, excluded)
# ------------------------------------------------

df <- df %>%
  mutate(
    therapy_group = case_when(
      as.character(ab_change) == "None"      ~ "Adequate",
      as.character(ab_change) %in%
        c("De-escalate", "Escalate", "Stop") ~ "Inadequate / NGS-guided",
      TRUE                                    ~ NA_character_
    )
  )

# ------------------------------------------------
# Sex
# ------------------------------------------------

df <- df %>%
  mutate(
    sex = factor(
      case_when(
        str_to_lower(as.character(sex)) %in% c("m", "male", "männlich") ~ "m",
        str_to_lower(as.character(sex)) %in% c("w", "f", "female", "weiblich") ~ "w",
        TRUE ~ NA_character_
      ),
      levels = c("m", "w")
    )
  )

# ------------------------------------------------
# Antibiotic at sampling — classify free-text antiinfektiva
# into the 11 publication categories
# ------------------------------------------------

classify_abx <- function(x) {
  x_lo <- str_to_lower(as.character(x))
  case_when(
    is.na(x) | str_trim(x_lo) == ""
      ~ "No antibiotic therapy",

    # Antifungal (no antibacterial)
    str_detect(x_lo, "caspo|caspofungin|anidulafungin|amphotericin|fluconazol|voriconazol") &
      !str_detect(x_lo, "pip|taz|mero|ceftri|ampi|amp|levo|cipro|vanco|linezolid|clinda|penicillin|cefaz|cefuro|ceftaz|cotrim|doxy")
      ~ "Antifungal therapy",

    # Antiviral (no antibacterial)
    str_detect(x_lo, "aciclovir|acilovir|acyclovir") &
      !str_detect(x_lo, "pip|taz|mero|ceftri|ampi|amp|levo|cipro|vanco|linezolid|clinda|penicillin|cefaz|cefuro|ceftaz|cotrim|doxy")
      ~ "Antiviral therapy",

    # Carbapenem / reserve broad-spectrum
    str_detect(x_lo, "meropenem|menopenem|imipenem|ertapenem|ceftazidim.*avibactam|ceftazidim/avibactam")
      ~ "Carbapenem / reserve broad-spectrum therapy",

    # Gram-positive reserve (linezolid, vancomycin — without carbapenem)
    str_detect(x_lo, "linezolid|vancomycin|vanco|daptomycin|staphylex|flucloxacillin") &
      !str_detect(x_lo, "meropenem|menopenem|imipenem")
      ~ "Gram-positive reserve therapy",

    # Clindamycin / toxin-suppression
    str_detect(x_lo, "clinda") &
      !str_detect(x_lo, "meropenem|pip|taz|ceftri|ceftaz|levo|cipro|ampi|amp|vanco|linezolid")
      ~ "Clindamycin / toxin-suppression coverage",

    # Broad-spectrum beta-lactam (pip/taz, without carbapenem)
    str_detect(x_lo, "pip.*taz|piperacillin.*tazobactam|pip/taz|piptaz") &
      !str_detect(x_lo, "meropenem|menopenem|imipenem")
      ~ "Broad-spectrum beta-lactam therapy",

    # Cephalosporin + anaerobic coverage
    (str_detect(x_lo, "ceftriax|ceftr|cetfr|cetriax") &
       str_detect(x_lo, "metro|metronidazol|metronidaxol|metromidazle|clont|met\\b")) |
    str_detect(x_lo, "cef/clont|cef/met|cetfr.*clont")
      ~ "Cephalosporin / anaerobic coverage",

    # Quinolone / macrolide / atypical
    str_detect(x_lo, "levo|cipro|ciproflox|roxi|roxithromycin|clarithromycin|erythromycin|doxycyclin|cotrimoxazol|cotrim|azithromycin") &
      !str_detect(x_lo, "meropenem|pip.*taz|piperacillin") &
      !str_detect(x_lo, "ceftriax|ceftr|amp.*sulb|ampi.*sulb|amoclav|cefazolin|cephazolin|cefuroxim")
      ~ "Quinolone / macrolide / atypical coverage",

    # Standard beta-lactam
    str_detect(x_lo, "amp.*sulb|ampi.*sulb|amoclav|ampicillin.*sulbactam|amp/sulb|ampi/sulb") |
    str_detect(x_lo, "^ceftriax|^ceftr|^cetriax|^cetfr") & !str_detect(x_lo, "metro|clont|met\\b|mero") |
    str_detect(x_lo, "cefazolin|cephazolin|cefuroxim|cefurox|cef\\b") & !str_detect(x_lo, "metro|clont") |
    str_detect(x_lo, "^penicillin|^pen g|^pen\\b") & !str_detect(x_lo, "clinda|mero") |
    str_detect(x_lo, "^flucloxacillin$|^staphylex$")
      ~ "Standard beta-lactam therapy",

    TRUE ~ "Other / unclear antimicrobial therapy"
  )
}

df <- df %>%
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

# Show any "Other" rows for review
other_rows <- df %>%
  filter(abx_at_sampling == "Other / unclear antimicrobial therapy") %>%
  select(job_id, antiinfektiva, abx_at_sampling)
if (nrow(other_rows) > 0) {
  cat("\nRows classified as 'Other' — review:\n")
  print(other_rows)
}

# ------------------------------------------------
# Pathogen family grouping
# ------------------------------------------------

df <- df %>%
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
# Build publication dataframe
# ------------------------------------------------

pub_df <- df %>%
  transmute(
    fallnummer = as.character(fallnummer),          # unique patient ID — not in table
    therapy_group = factor(therapy_group,
                           levels = c("Adequate", "Inadequate / NGS-guided")),
    age      = to_num(age),
    sex      = sex,
    icu_days = to_num(icu_days),
    hosp_days = to_num(hosp_days),

    # ── 28-day ICU mortality ──────────────────────────────────
    # Definition: patient died (Died == Yes) AND ICU days == Hosp. days
    # (never left ICU — death occurred within the ICU stay)
    icu_mortality_28d = factor(
      case_when(
        # Died AND spent entire hospitalisation in ICU → ICU death
        mortality_to_binary(died) == 1L &
          !is.na(to_num(icu_days)) & !is.na(to_num(hosp_days)) &
          to_num(icu_days) == to_num(hosp_days)          ~ 1L,
        # Survived → no ICU death
        mortality_to_binary(died) == 0L                  ~ 0L,
        # Died but ICU days < Hosp. days → transferred out of ICU first
        mortality_to_binary(died) == 1L                  ~ 0L,
        TRUE                                              ~ NA_integer_
      ),
      levels = c(0L, 1L),
      labels = c("Survived", "Died")
    ),

    # ── Hospital mortality ────────────────────────────────────
    hospital_mortality = factor(
      mortality_to_binary(died),
      levels = c(0L, 1L),
      labels = c("Survived", "Died")
    ),

    infect_type           = factor(infect_type),
    abx_at_sampling       = abx_at_sampling,
    pathogen_family_group = pathogen_family_group
  ) %>%
  filter(!is.na(therapy_group))

cat("\nFinal N by group:\n")
print(table(pub_df$therapy_group, useNA = "ifany"))
cat("Total N:", nrow(pub_df), "\n")

# ── Fallnummer uniqueness check ──────────────────────────────
n_unique <- length(unique(pub_df$fallnummer))
dups     <- pub_df$fallnummer[duplicated(pub_df$fallnummer)]
if (length(dups) == 0) {
  cat(sprintf("\nUniqueness check: ALL %d Fallnummern unique — no duplicates. Analysis is on unique patients.\n",
              n_unique))
} else {
  cat(sprintf("\nWARNING: %d duplicate Fallnummer(n) detected — review before publication:\n",
              length(dups)))
  print(pub_df %>% filter(fallnummer %in% dups) %>%
          select(fallnummer, therapy_group, icu_days, hosp_days, hospital_mortality))
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
  select(-fallnummer) %>%                          # ID column — not a table variable
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
      infect_type           ~ "Infection type",
      abx_at_sampling       ~ "Antibiotic therapy at sampling",
      pathogen_family_group ~ "Detected pathogen group"
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
# Export → Word
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

print(doc, target = "Table1_Publication_Summary.docx")
cat("\nDone. File written: Table1_Publication_Summary.docx\n")
