# ============================================================
# COMPETING-RISK SENSITIVITY ANALYSIS
#
# Event of interest:
#   Recurrence
#
# Competing event:
#   Death without prior recurrence
#
# Fine-Gray regression using tidycmprsk
# ============================================================


library(tidyverse)
library(tidycmprsk)
library(broom)
library(officer)
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
  recursive = TRUE,
  showWarnings = FALSE
)

dat <- readRDS(
  file.path(
    data_dir,
    "RCC_sjTREC_TRB.rds"
  )
)



# ============================================================
# 3. CREATE COMPETING-RISK OUTCOME
#
# cr_status:
#   0 = censored
#   1 = recurrence
#   2 = death without prior recurrence
#
# cr_time:
#   time to recurrence, competing death, or censoring
# ============================================================

dat_cr <- dat %>%
  mutate(
    cr_status = case_when(
      # Missing event indicators remain missing
      is.na(rec) | is.na(death) ~ NA_character_,

      # Any observed recurrence is the event of interest
      rec == 1 ~ "1",

      # Death without prior recurrence is the competing event
      rec == 0 & death == 1 ~ "2",

      # Alive and recurrence-free
      rec == 0 & death == 0 ~ "0"
    ),
    cr_time = case_when(
      cr_status == "1" ~ t_rec,
      cr_status == "2" ~ t_os,
      cr_status == "0" ~ t_rec,
      TRUE ~ NA_real_
    ),
    cr_status = factor(
      cr_status,
      levels = c("0", "1", "2")
    )
  )


# ============================================================
# 6. sjTREC FINE-GRAY MODEL
# ============================================================

sjtrec_cr <- dat_cr %>%
  filter(
    !is.na(sjTREC_group),
    !is.na(cr_time),
    !is.na(cr_status)
  )

fg_sjtrec <- crr(
  Surv(
    cr_time,
    cr_status
  ) ~
    sjTREC_group +
    Age +
    Sex +
    Leibovich,
  data = sjtrec_cr
)




# ============================================================
# 7. TRB FINE-GRAY MODEL
# ============================================================

trb_cr <- dat_cr %>%
  filter(
    !is.na(TRB_z),
    !is.na(cr_time),
    !is.na(cr_status)
  )


fg_trb <- crr(
  Surv(
    cr_time,
    cr_status
  ) ~
    TRB_z +
    Age +
    Sex +
    Leibovich,
  data = trb_cr
)

# ============================================================
# SUPPLEMENTARY TABLE
# ============================================================

# ============================================================
# 1. EXTRACT sjTREC MODEL
# ============================================================

sjtrec_tab <- tidy(
  fg_sjtrec,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(
    Variable = case_when(
      term == "sjTREC_groupHigh sjTREC" ~
        "High vs low sjTREC",
      term == "Age" ~
        "Age, per year",
      term == "SexWoman" ~
        "Sex: woman vs man",
      term == "Leibovichintermediate" ~
        "Leibovich: intermediate vs low",
      term == "Leibovichhigh" ~
        "Leibovich: high vs low",
      TRUE ~ term
    ),
    SHR = sprintf(
      "%.2f",
      estimate
    ),
    CI = sprintf(
      "%.2f–%.2f",
      conf.low,
      conf.high
    ),
    P = case_when(
      p.value < 0.001 ~
        "<0.001",
      TRUE ~
        sprintf(
          "%.3f",
          p.value
        )
    )
  ) %>%
  select(
    Variable,
    SHR,
    CI,
    P
  )


# ============================================================
# 2. EXTRACT TRB MODEL
# ============================================================

trb_tab <- tidy(
  fg_trb,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  mutate(
    Variable = case_when(
      term == "TRB_z" ~
        "TCRβ diversity per 1 SD",
      term == "Age" ~
        "Age, per year",
      term == "SexWoman" ~
        "Sex: woman vs man",
      term == "Leibovichintermediate" ~
        "Leibovich: intermediate vs low",
      term == "Leibovichhigh" ~
        "Leibovich: high vs low",
      TRUE ~ term
    ),
    SHR = sprintf(
      "%.2f",
      estimate
    ),
    CI = sprintf(
      "%.2f–%.2f",
      conf.low,
      conf.high
    ),
    P = case_when(
      p.value < 0.001 ~
        "<0.001",
      TRUE ~
        sprintf(
          "%.3f",
          p.value
        )
    )
  ) %>%
  select(
    Variable,
    SHR,
    CI,
    P
  )


# ============================================================
# 3. EVENT COUNTS
# ============================================================

sjtrec_N <- nrow(
  sjtrec_cr
)

sjtrec_rec <- sum(
  sjtrec_cr$cr_status == "1"
)

sjtrec_death <- sum(
  sjtrec_cr$cr_status == "2"
)


trb_N <- nrow(
  trb_cr
)

trb_rec <- sum(
  trb_cr$cr_status == "1"
)

trb_death <- sum(
  trb_cr$cr_status == "2"
)


# ============================================================
# 4. ADD GROUP HEADERS
# ============================================================

sjtrec_label <- paste0(
  "sjTREC model (n = ",
  sjtrec_N,
  "; ",
  sjtrec_rec,
  " recurrences; ",
  sjtrec_death,
  " competing deaths)"
)


trb_label <- paste0(
  "TCRβ model (n = ",
  trb_N,
  "; ",
  trb_rec,
  " recurrences; ",
  trb_death,
  " competing deaths)"
)


tab_out <- bind_rows(
  tibble(
    Variable = sjtrec_label,
    SHR = "",
    CI = "",
    P = "",
    is_group = TRUE
  ),
  sjtrec_tab %>%
    mutate(
      is_group = FALSE
    ),
  tibble(
    Variable = trb_label,
    SHR = "",
    CI = "",
    P = "",
    is_group = TRUE
  ),
  trb_tab %>%
    mutate(
      is_group = FALSE
    )
)


group_rows <- which(
  tab_out$is_group
)


# ============================================================
# 5. FLEXTABLE
# ============================================================

line <- fp_border(
  color = "black",
  width = 1
)


ft <- flextable(
  tab_out,
  col_keys = c(
    "Variable",
    "SHR",
    "CI",
    "P"
  )
) %>%
  set_header_labels(
    Variable = "Variable",
    SHR = "SHR",
    CI = "95% CI",
    P = "p value"
  ) %>%
  flextable::font(
    fontname = "Arial",
    part = "all"
  ) %>%
  fontsize(
    size = 10,
    part = "all"
  ) %>%
  bold(
    part = "header"
  ) %>%
  border_remove() %>%
  hline_top(
    part = "header",
    border = line
  ) %>%
  hline_bottom(
    part = "header",
    border = line
  ) %>%
  hline_bottom(
    part = "body",
    border = line
  ) %>%
  merge_h_range(
    i = group_rows,
    j1 = "Variable",
    j2 = "P"
  ) %>%
  bold(
    i = group_rows,
    part = "body"
  ) %>%
  bg(
    i = group_rows,
    bg = "#EDEDED",
    part = "body"
  ) %>%
  align(
    j = 1,
    align = "left",
    part = "all"
  ) %>%
  align(
    j = 2:4,
    align = "center",
    part = "all"
  ) %>%
  width(
    j = 1,
    width = 3.3
  ) %>%
  width(
    j = 2,
    width = 0.75
  ) %>%
  width(
    j = 3,
    width = 1.1
  ) %>%
  width(
    j = 4,
    width = 0.75
  ) %>%
  set_table_properties(
    layout = "fixed"
  )

# ============================================================
# 6. FOOTNOTES
# ============================================================

ft <- ft %>%
  add_footer_lines(
    values = c(
      paste0(
        "Fine-Gray regression was used with recurrence as the event of interest ",
        "and death without prior recurrence as the competing event."
      ),
      paste0(
        "Both models were adjusted for age, sex, and Leibovich risk group."
      ),
      paste0(
        "sjTREC was dichotomised at the cohort median. ",
        "TCRβ diversity was modelled continuously per 1 SD increase."
      ),
      paste0(
        "Low sjTREC, male sex, and low Leibovich risk were reference categories."
      ),
      "SHR = subdistribution hazard ratio; CI = confidence interval."
    )
  ) %>%
  flextable::font(
    fontname = "Arial",
    part = "footer"
  ) %>%
  fontsize(
    size = 8,
    part = "footer"
  )


# ============================================================
# 7. TITLE
# ============================================================

title <- fpar(
  ftext(
    paste0(
      "Fine-Gray competing-risk sensitivity analyses ",
      "for recurrence"
    ),
    fp_text(
      bold = TRUE,
      font.size = 11,
      font.family = "Arial"
    )
  )
)


# ============================================================
# 8. SAVE WORD FILE
# ============================================================

doc <- read_docx() %>%
  body_add_fpar(
    title
  ) %>%
  body_add_flextable(
    ft
  )

print(
  doc,
  target = file.path(
    out_dir,
    "Supplementary_competing_risk_full_models.docx"
  )
)
