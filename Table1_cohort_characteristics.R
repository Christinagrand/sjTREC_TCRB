# ============================================================
# TABLE 1: FINAL LOCALIZED RCC COHORT
# ============================================================

library(tidyverse)
library(gtsummary)
library(flextable)


# ============================================================
# 1. PATHS
# ============================================================

data_dir <- "data"
out_dir <- file.path(
  "results",
  "tables"
)

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# 2. LOAD FINAL DATASET
# ============================================================

dat <- readRDS(
  file.path(
    data_dir,
    "RCC_sjTREC_TRB.rds"
  )
)

# ============================================================
# 3. PREPARE TABLE VARIABLES
# ============================================================

histology_labels <- c(
  ccRCC = "Clear cell",
  pRCC = "Papillary",
  chRCC = "Chromophobe"
)


patients <- dat %>%
  transmute(
    age = as.numeric(age),
    sex = factor(
      gender,
      levels = c(
        "Man",
        "Woman"
      ),
      labels = c(
        "Men",
        "Women"
      )
    ),
    histology = coalesce(
      histology_labels[subtype_PatoMcode],
      "Not recorded"
    ),
    stage = case_when(
      as.character(stadium_con) == "1" ~ "I",
      as.character(stadium_con) == "2" ~ "II",
      as.character(stadium_con) == "3" ~ "III",
      as.character(stadium_con) == "local" ~ "Localized, stage not specified"
    ),
    leibovich = case_when(
      as.character(Leibscore_gr) == "low" ~ "Low",
      as.character(Leibscore_gr) == "intermediate" ~ "Intermediate",
      as.character(Leibscore_gr) == "high" ~ "High"
    ),
    follow_up = as.numeric(time_year),
    recurrence = as.numeric(rec_at_any_time_sum) == 1,
    death =
      str_to_lower(
        str_trim(
          survival_status
        )
      ) == "dead",
    has_trec = has_sjTREC,
    has_trb = has_TRB,
    histology = factor(
      histology,
      levels = c(
        "Clear cell",
        "Papillary",
        "Chromophobe",
        "Not recorded"
      )
    ),
    stage = factor(
      stage,
      levels = c(
        "I",
        "II",
        "III",
        "Localized, stage not specified"
      )
    ),
    leibovich = factor(
      leibovich,
      levels = c(
        "Low",
        "Intermediate",
        "High"
      )
    )
  )


# ============================================================
# 4. TABLE 1
# ============================================================

table1 <- patients %>%
  select(
    age,
    sex,
    histology,
    stage,
    leibovich,
    follow_up,
    recurrence,
    death,
    has_trec,
    has_trb
  ) %>%
  tbl_summary(
    label = list(
      age ~ "Age, yr",
      sex ~ "Sex",
      histology ~ "Histology",
      stage ~ "Tumour stage",
      leibovich ~ "Leibovich risk group",
      follow_up ~ "Follow-up, yr",
      recurrence ~ "Recurrence",
      death ~ "Death",
      has_trec ~ "sjTREC measurement available",
      has_trb ~ "TCRβ repertoire available"
    ),
    type = list(
      sex ~ "categorical",
      histology ~ "categorical",
      stage ~ "categorical",
      leibovich ~ "categorical",
      recurrence ~ "dichotomous",
      death ~ "dichotomous",
      has_trec ~ "dichotomous",
      has_trb ~ "dichotomous"
    ),
    value = list(
      recurrence ~ TRUE,
      death ~ TRUE,
      has_trec ~ TRUE,
      has_trb ~ TRUE
    ),
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = list(
      age ~ 1,
      follow_up ~ 1,
      all_categorical() ~ c(0, 1)
    ),
    missing = "ifany",
    missing_text = "Not recorded"
  ) %>%
  modify_header(
    label = "**Characteristic**",
    stat_0 = "**N = {N}**"
  ) %>%
  bold_labels()


# Show table in RStudio
table1


# ============================================================
# 5. CONVERT TO FLEXTABLE
# ============================================================

table1_ft <- table1 %>%
  as_flex_table() %>%
  flextable::font(
    fontname = "Times New Roman",
    part = "all"
  ) %>%
  fontsize(
    size = 10,
    part = "all"
  ) %>%
  padding(
    padding.top = 2.5,
    padding.bottom = 2.5,
    part = "all"
  ) %>%
  line_spacing(
    space = 1,
    part = "all"
  ) %>%
  autofit()


# Show formatted flextable
table1_ft



# ============================================================
# 6. SAVE
# ============================================================

save_as_docx(
  table1_ft,
  path = file.path(
    out_dir,
    "Table1_cohort_characteristics.docx"
  )
)
