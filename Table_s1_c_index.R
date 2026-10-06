# ============================================================
# sjTREC C-INDEX ANALYSIS
#
# Recurrence-free survival
#
# Models:
#   1. Leibovich only
#   2. Leibovich + sjTREC high/low
#
# Uses all eligible patients with sjTREC.
# TRB is NOT required.
#
# C-index is optimism-corrected using bootstrap resampling.
# ============================================================

library(tidyverse)
library(survival)
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


dat <- dat %>%
  filter(
    !is.na(sjTREC),
    !is.na(t_rec),
    !is.na(rec)
  )


# ============================================================
# 4. MODEL FITTING FUNCTION
# ============================================================

fit_models <- function(d) {
  list(
    "Leibovich only" =
      coxph(
        Surv(t_rec, rec) ~ Leibovich,
        data = d,
        x = TRUE
      ),
    "Leibovich + sjTREC" =
      coxph(
        Surv(t_rec, rec) ~ Leibovich + sjTREC_group,
        data = d,
        x = TRUE
      )
  )
}


# ============================================================
# 5. C-INDEX FUNCTION
#
# Predict risk from fitted model in dataset d and calculate
# Harrell's concordance index.
# ============================================================

get_cindex <- function(fit, d) {
  risk <- predict(
    fit,
    newdata = d,
    type = "lp"
  )

  concordance(
    Surv(t_rec, rec) ~ risk,
    data = d,
    reverse = TRUE
  )$concordance
}


# ============================================================
# 6. BOOTSTRAP OPTIMISM CORRECTION
#
# Corrected C =
# apparent C - mean(bootstrap C - test C)
# ============================================================

get_corrected_cindex <- function(
  d,
  B = 2000,
  seed = 1
) {
  set.seed(seed)


  # ----------------------------------------------------------
  # Apparent performance
  # Fit on original data and evaluate on original data
  # ----------------------------------------------------------

  models <- fit_models(d)

  apparent <- sapply(
    models,
    function(fit) {
      get_cindex(
        fit,
        d
      )
    }
  )


  # ----------------------------------------------------------
  # Bootstrap optimism
  # ----------------------------------------------------------

  optimism <- matrix(
    NA_real_,
    nrow = B,
    ncol = length(models),
    dimnames = list(
      NULL,
      names(models)
    )
  )


  for (b in seq_len(B)) {
    boot_index <- sample(
      seq_len(nrow(d)),
      size = nrow(d),
      replace = TRUE
    )

    boot <- d[
      boot_index, ,
      drop = FALSE
    ]


    boot_models <- fit_models(
      boot
    )


    optimism[b, ] <- sapply(
      boot_models,
      function(fit) {
        c_boot <- get_cindex(
          fit,
          boot
        )

        c_test <- get_cindex(
          fit,
          d
        )

        c_boot - c_test
      }
    )
  }


  # ----------------------------------------------------------
  # Mean optimism and corrected C
  # ----------------------------------------------------------

  mean_optimism <- colMeans(
    optimism,
    na.rm = TRUE
  )


  tibble(
    Model = names(models),
    Apparent_C = apparent,
    Optimism = mean_optimism,
    Corrected_C =
      apparent - mean_optimism
  )
}


# ============================================================
# 7. RUN ANALYSIS FUNCTION
# ============================================================

run_analysis <- function(
  d,
  analysis_name,
  B = 2000
) {
  models <- fit_models(
    d
  )


  cindex_results <- get_corrected_cindex(
    d,
    B = B
  )


  # ----------------------------------------------------------
  # Likelihood ratio test
  # Does sjTREC improve model fit beyond Leibovich?
  # ----------------------------------------------------------

  lrt_p <- anova(
    models[[1]],
    models[[2]],
    test = "LRT"
  )[2, "Pr(>|Chi|)"]


  # ----------------------------------------------------------
  # Baseline corrected C
  # ----------------------------------------------------------

  baseline_c <- cindex_results$Corrected_C[
    cindex_results$Model ==
      "Leibovich only"
  ]


  # ----------------------------------------------------------
  # Results
  # ----------------------------------------------------------

  cindex_results %>%
    mutate(
      Analysis = analysis_name,
      N = nrow(d),
      Recurrences = sum(d$rec),
      Delta_corrected_C =
        Corrected_C - baseline_c,
      LRT_p = lrt_p
    ) %>%
    select(
      Analysis,
      Model,
      N,
      Recurrences,
      Apparent_C,
      Optimism,
      Corrected_C,
      Delta_corrected_C,
      LRT_p
    )
}


# ============================================================
# 8. ALL LOCALISED PATIENTS
# ============================================================

results_all <- run_analysis(
  dat,
  analysis_name = "All patients",
  B = 2000
)


# ============================================================
# 9. INTERMEDIATE + HIGH LEIBOVICH
#
# Uses same sjTREC median cut-off as full cohort.
# ============================================================

dat_ih <- dat %>%
  filter(
    Leibovich %in%
      c(
        "intermediate",
        "high"
      )
  ) %>%
  droplevels()

results_ih <- run_analysis(
  dat_ih,
  analysis_name = "Intermediate/high",
  B = 2000
)


# ============================================================
# 10. COMBINE RESULTS
# ============================================================

results <- bind_rows(
  results_all,
  results_ih
)


# ============================================================
# 12. FORMAT WORD TABLE
# ============================================================

tab <- results %>%
  mutate(
    Apparent_C =
      sprintf(
        "%.4f",
        Apparent_C
      ),
    Optimism =
      sprintf(
        "%.4f",
        Optimism
      ),
    Corrected_C =
      sprintf(
        "%.4f",
        Corrected_C
      ),
    Delta =
      case_when(
        Model == "Leibovich only" ~
          "Reference",
        Delta_corrected_C > 0 ~
          sprintf(
            "+%.4f",
            Delta_corrected_C
          ),
        Delta_corrected_C < 0 ~
          sprintf(
            "−%.4f",
            abs(
              Delta_corrected_C
            )
          ),
        TRUE ~
          "0.0000"
      ),
    LRT_p =
      case_when(
        Model == "Leibovich only" ~
          "—",
        LRT_p < 0.001 ~
          "<0.001",
        TRUE ~
          sprintf(
            "%.3f",
            LRT_p
          )
      )
  )


# ============================================================
# 13. ADD POPULATION LABEL ROWS
# ============================================================

tab_list <- list()


for (a in unique(tab$Analysis)) {
  rows_a <- tab %>%
    filter(
      Analysis == a
    )


  label <- paste0(
    if_else(
      a == "Intermediate/high",
      "Intermediate/high Leibovich",
      a
    ),
    " (n = ",
    rows_a$N[1],
    "; ",
    rows_a$Recurrences[1],
    " recurrences)"
  )


  tab_list[[a]] <- bind_rows(
    tibble(
      Model = label,
      is_group = TRUE
    ),
    rows_a %>%
      mutate(
        is_group = FALSE
      )
  )
}


tab_out <- bind_rows(
  tab_list
) %>%
  mutate(
    across(
      c(
        Apparent_C,
        Optimism,
        Corrected_C,
        Delta,
        LRT_p
      ),
      ~ replace_na(
        .x,
        ""
      )
    )
  )


group_rows <- which(
  tab_out$is_group
)


# ============================================================
# 14. BUILD FLEXTABLE
# ============================================================

line <- fp_border(
  color = "black",
  width = 1
)


ft <- flextable(
  tab_out,
  col_keys = c(
    "Model",
    "Apparent_C",
    "Optimism",
    "Corrected_C",
    "Delta",
    "LRT_p"
  )
) %>%
  set_header_labels(
    Model =
      "Model",
    Apparent_C =
      "Apparent C",
    Optimism =
      "Optimism",
    Corrected_C =
      "Corrected C",
    Delta =
      "\u0394C vs Leibovich",
    LRT_p =
      "LRT p"
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
    j1 = "Model",
    j2 = "LRT_p"
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
    j = 2:6,
    align = "center",
    part = "all"
  ) %>%
  align(
    j = 1,
    align = "left",
    part = "all"
  ) %>%
  width(
    j = 1,
    width = 2.1
  ) %>%
  width(
    j = 2:6,
    width = 0.83
  ) %>%
  set_table_properties(
    layout = "fixed"
  ) %>%
  padding(
    padding.top = 3,
    padding.bottom = 3,
    part = "all"
  )


# ============================================================
# 15. FOOTNOTES
# ============================================================

ft <- ft %>%
  add_footer_lines(
    values = c(
      paste0(
        "Corrected C = apparent C minus mean optimism ",
        "from 2,000 bootstrap resamples."
      ),
      paste0(
        "\u0394C = corrected C of the model minus corrected C ",
        "of the Leibovich-only model in the same population."
      ),
      paste0(
        "LRT p: likelihood ratio test comparing the model ",
        "including sjTREC with the Leibovich-only model."
      ),
      paste0(
        "sjTREC dichotomised at the median of all patients ",
        "in the sjTREC analysis cohort; the same cut-off ",
        "is used in both populations."
      )
    )
  ) %>%
  flextable::font(
    fontname = "Arial",
    part = "footer"
  ) %>%
  fontsize(
    size = 8,
    part = "footer"
  ) %>%
  align(
    align = "left",
    part = "footer"
  )


# ============================================================
# 16. SAVE WORD TABLE
# ============================================================

title <- fpar(
  ftext(
    paste0(
      "C-index for recurrence: ",
      "Leibovich score alone and with sjTREC"
    ),
    fp_text(
      bold = TRUE,
      font.size = 11,
      font.family = "Arial"
    )
  )
)


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
    "C_index_sjTREC_only.docx"
  )
)
