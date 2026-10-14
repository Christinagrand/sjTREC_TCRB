# ============================================================
# CONTINUOUS sjTREC ANALYSIS
# OS + RFS
# Adjusted for age, sex, and Leibovich score
# ============================================================

library(tidyverse)
library(survival)
library(broom)
library(ggplot2)
library(patchwork)



# ============================================================
# 1. LOAD DATA
# ============================================================

data_dir <- "data"
out_dir <- file.path(
  "results",
  "figures"
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
# 2. OS DATA
# ============================================================

os_data <- dat %>%
  filter(
    !is.na(log2_sjTREC),
    !is.na(t_os),
    !is.na(death)
  )

# ============================================================
# 3. CONTINUOUS COX - OS
# ============================================================

cox_os <- coxph(
  Surv(
    t_os,
    death
  ) ~
    log2_sjTREC +
    Age +
    Sex +
    Leibovich,
  data = os_data
)



# sjTREC result only
os_result <- tidy(
  cox_os,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(
    term == "log2_sjTREC"
  ) %>%
  select(
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )


# ============================================================
# 4. RFS DATA
# ============================================================

rfs_data <- dat %>%
  filter(
    !is.na(log2_sjTREC),
    !is.na(t_rec),
    !is.na(rec)
  )


# ============================================================
# 5. CONTINUOUS COX - RFS
# ============================================================

cox_rfs <- coxph(
  Surv(
    t_rec,
    rec
  ) ~
    log2_sjTREC +
    Age +
    Sex +
    Leibovich,
  data = rfs_data
)



# sjTREC result only
rfs_result <- tidy(
  cox_rfs,
  exponentiate = TRUE,
  conf.int = TRUE
) %>%
  filter(
    term == "log2_sjTREC"
  ) %>%
  select(
    term,
    estimate,
    conf.low,
    conf.high,
    p.value
  )


# ============================================================
# FOREST PLOTS - FULL MULTIVARIABLE MODELS
# Forest plot + separate HR/CI and p-value columns
# ============================================================


# ============================================================
# 6. PREPARE COX RESULTS
# ============================================================

prepare_forest_data <- function(model) {
  tidy(
    model,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    mutate(
      Variable = case_when(
        term == "log2_sjTREC" ~ "sjTREC, log2",
        term == "Age" ~ "Age",
        term == "SexWoman" ~ "Sex: Woman vs man",
        term == "SexFemale" ~ "Sex: Woman vs man",
        term == "Leibovichintermediate" ~
          "Leibovich: Intermediate vs low",
        term == "Leibovichhigh" ~
          "Leibovich: High vs low",
        TRUE ~ term
      ),
      HR_CI = sprintf(
        "%.2f (%.2f–%.2f)",
        estimate,
        conf.low,
        conf.high
      ),
      P_value = case_when(
        p.value < 0.001 ~ "<0.001",
        TRUE ~ sprintf(
          "%.3f",
          p.value
        )
      )
    )
}


os_forest_data <- prepare_forest_data(
  cox_os
)

rfs_forest_data <- prepare_forest_data(
  cox_rfs
)


# ============================================================
# 7. FUNCTION TO MAKE ONE FOREST + TABLE
# ============================================================

make_forest <- function(
  df,
  title
) {
  # Keep model order, with first variable at top
  df <- df %>%
    mutate(
      Variable = factor(
        Variable,
        levels = rev(
          Variable
        )
      )
    )


  # ----------------------------------------------------------
  # FOREST PANEL
  # ----------------------------------------------------------

  p_forest <- ggplot(
    df,
    aes(
      x = estimate,
      y = Variable
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      linewidth = 0.5
    ) +
    geom_errorbarh(
      aes(
        xmin = conf.low,
        xmax = conf.high
      ),
      height = 0.14,
      linewidth = 0.7
    ) +
    geom_point(
      size = 3
    ) +
    scale_x_log10() +
    labs(
      title = title,
      x = "Hazard ratio",
      y = NULL
    ) +
    theme_classic(
      base_size = 12
    ) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 12,
        hjust = 0
      ),
      axis.text.y = element_text(
        size = 10
      ),
      axis.title.x = element_text(
        size = 10
      ),
      plot.margin = margin(
        t = 10,
        r = 10,
        b = 10,
        l = 5
      )
    )


  # ----------------------------------------------------------
  # RESULTS TABLE
  # ----------------------------------------------------------

  p_table <- ggplot(
    df,
    aes(
      y = Variable
    )
  ) +
    geom_text(
      aes(
        x = 1,
        label = HR_CI
      ),
      hjust = 0.5,
      size = 3.4
    ) +
    geom_text(
      aes(
        x = 2,
        label = P_value
      ),
      hjust = 0.5,
      size = 3.4
    ) +
    annotate(
      "text",
      x = 1,
      y = length(
        unique(
          df$Variable
        )
      ) + 0.65,
      label = "HR (95% CI)",
      fontface = "bold",
      size = 3.4
    ) +
    annotate(
      "text",
      x = 2,
      y = length(
        unique(
          df$Variable
        )
      ) + 0.65,
      label = "p value",
      fontface = "bold",
      size = 3.4
    ) +
    scale_x_continuous(
      limits = c(
        0.5,
        2.5
      )
    ) +
    scale_y_discrete(
      limits = levels(
        df$Variable
      )
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_void() +
    theme(
      plot.margin = margin(
        t = 10,
        r = 5,
        b = 10,
        l = 0
      )
    )


  # ----------------------------------------------------------
  # COMBINE FOREST + TABLE
  # ----------------------------------------------------------

  p_forest +
    p_table +
    plot_layout(
      widths = c(
        1.7,
        1
      )
    )
}


# ============================================================
# 8. MAKE RFS AND OS PLOTS
# ============================================================

plot_rfs <- make_forest(
  rfs_forest_data,
  "Recurrence-free survival"
)

plot_os <- make_forest(
  os_forest_data,
  "Overall survival"
)



# ============================================================
# 9. COMBINE AS SUPPLEMENTARY FIGURE
#
# A = Recurrence-free survival
# B = Overall survival
# ============================================================

combined_forest <-
  wrap_elements(
    full = plot_rfs
  ) /
  wrap_elements(
    full = plot_os
  ) +
  plot_annotation(
    tag_levels = "A"
  )


ggsave(
  filename = file.path(
    out_dir,
    "Continuous_sjTREC_full_model_forest_combined.pdf"
  ),
  plot = combined_forest,
  width = 8,
  height = 9
)
