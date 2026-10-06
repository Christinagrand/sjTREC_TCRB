library(tidyverse)
library(survival)
library(broom)


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


dat_ih <- dat %>%
  filter(
    tolower(Leibscore_gr) %in% c(
      "intermediate",
      "high"
    )
  )


# ============================================================
# COMPLETE-CASE DATA
# ============================================================

os_full <- dat %>%
  filter(
    complete.cases(
      TRB,
      normalized_shannon,
      t_os,
      death
    )
  )

rfs_full <- dat %>%
  filter(
    complete.cases(
      TRB,
      normalized_shannon,
      t_rec,
      rec
    )
  )

os_ih <- dat_ih %>%
  filter(
    complete.cases(
      TRB,
      normalized_shannon,
      t_os,
      death
    )
  )

rfs_ih <- dat_ih %>%
  filter(
    complete.cases(
      TRB,
      normalized_shannon,
      t_rec,
      rec
    )
  )


# ============================================================
# FUNCTION
# ============================================================

fit_marker <- function(
  data,
  time,
  event,
  marker,
  metric,
  outcome,
  cohort
) {
  fit <- coxph(
    as.formula(
      paste0(
        "Surv(", time, ",", event, ") ~ ",
        marker, " + Age + Sex + Leibovich"
      )
    ),
    data = data
  )

  tidy(
    fit,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term == marker
    ) %>%
    transmute(
      Cohort = cohort,
      Outcome = outcome,
      Metric = metric,
      N = fit$n,
      Events = fit$nevent,
      HR = estimate,
      CI_low = conf.low,
      CI_high = conf.high,
      p = p.value
    )
}


# ============================================================
# FULL ANALYSES
# ============================================================

results <- bind_rows(
  # Full cohort: recurrence-free survival
  fit_marker(
    rfs_full,
    "t_rec",
    "rec",
    "TRB_z",
    "Simpson",
    "RFS",
    "Full cohort"
  ),
  fit_marker(
    rfs_full,
    "t_rec",
    "rec",
    "shannon_z",
    "Normalized Shannon",
    "RFS",
    "Full cohort"
  ),

  # Full cohort: overall survival
  fit_marker(
    os_full,
    "t_os",
    "death",
    "TRB_z",
    "Simpson",
    "OS",
    "Full cohort"
  ),
  fit_marker(
    os_full,
    "t_os",
    "death",
    "shannon_z",
    "Normalized Shannon",
    "OS",
    "Full cohort"
  ),

  # Intermediate/high: recurrence-free survival
  fit_marker(
    rfs_ih,
    "t_rec",
    "rec",
    "TRB_z",
    "Simpson",
    "RFS",
    "Intermediate/high"
  ),
  fit_marker(
    rfs_ih,
    "t_rec",
    "rec",
    "shannon_z",
    "Normalized Shannon",
    "RFS",
    "Intermediate/high"
  ),

  # Intermediate/high: overall survival
  fit_marker(
    os_ih,
    "t_os",
    "death",
    "TRB_z",
    "Simpson",
    "OS",
    "Intermediate/high"
  ),
  fit_marker(
    os_ih,
    "t_os",
    "death",
    "shannon_z",
    "Normalized Shannon",
    "OS",
    "Intermediate/high"
  )
) %>%
  mutate(
    Metric = factor(
      Metric,
      levels = c(
        "Normalized Shannon",
        "Simpson"
      )
    ),
    Panel = case_when(
      Cohort == "Full cohort" &
        Outcome == "RFS" ~
        "A. Full cohort - recurrence-free survival",
      Cohort == "Full cohort" &
        Outcome == "OS" ~
        "B. Full cohort - overall survival",
      Cohort == "Intermediate/high" &
        Outcome == "RFS" ~
        "C. Intermediate/high - recurrence-free survival",
      Cohort == "Intermediate/high" &
        Outcome == "OS" ~
        "D. Intermediate/high - overall survival"
    ),
    Panel = factor(
      Panel,
      levels = c(
        "A. Full cohort - recurrence-free survival",
        "B. Full cohort - overall survival",
        "C. Intermediate/high - recurrence-free survival",
        "D. Intermediate/high - overall survival"
      )
    ),
    hr_label = sprintf(
      "%.2f (%.2f–%.2f)",
      HR,
      CI_low,
      CI_high
    ),
    p_label = case_when(
      p < 0.001 ~ "<0.001",
      TRUE ~ sprintf(
        "%.3f",
        p
      )
    )
  )


# ============================================================
# POSITIONS FOR TEXT COLUMNS
# ============================================================

hr_x <- max(
  results$CI_high,
  na.rm = TRUE
) * 1.30

p_x <- hr_x * 1.95

x_upper <- p_x * 1.45


# ============================================================
# HEADER DATA
# ============================================================

header_df <- tibble(
  Panel = factor(
    levels(results$Panel),
    levels = levels(results$Panel)
  ),
  Metric = factor(
    "Normalized Shannon",
    levels = levels(results$Metric)
  )
)


# ============================================================
# FOREST PLOT
# ============================================================

fig_trb <- ggplot(
  results,
  aes(
    x = HR,
    y = Metric
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    colour = "grey60"
  ) +
  geom_errorbar(
    aes(
      xmin = CI_low,
      xmax = CI_high
    ),
    orientation = "y",
    width = 0.08
  ) +
  geom_point(
    size = 2.8
  ) +

  # ----------------------------------------------------------
  # HR (95% CI) column
  # ----------------------------------------------------------

  geom_text(
    aes(
      x = hr_x,
      label = hr_label
    ),
    hjust = 0,
    size = 2.7
  ) +

  # ----------------------------------------------------------
  # p value column
  # ----------------------------------------------------------

  geom_text(
    aes(
      x = p_x,
      label = p_label
    ),
    hjust = 0,
    size = 2.7
  ) +

  # ----------------------------------------------------------
  # Column headers
  # ----------------------------------------------------------

  geom_text(
    data = header_df,
    aes(
      x = hr_x,
      y = Metric,
      label = "HR (95% CI)"
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -17.5,
    fontface = "bold",
    size = 2.8
  ) +
  geom_text(
    data = header_df,
    aes(
      x = p_x,
      y = Metric,
      label = "p value"
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -17.5,
    fontface = "bold",
    size = 2.8
  ) +
  facet_wrap(
    ~Panel,
    ncol = 2
  ) +
  scale_x_log10(
    limits = c(
      min(
        results$CI_low,
        na.rm = TRUE
      ) * 0.8,
      x_upper
    )
  ) +
  labs(
    x = "Hazard ratio per 1 SD higher TCRβ diversity",
    y = NULL
  ) +
  theme_classic(
    base_size = 10
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 10
    ),
    axis.text.y = element_text(
      size = 9
    ),
    plot.margin = margin(
      t = 8,
      r = 10,
      b = 5,
      l = 5
    )
  )


fig_trb


# ============================================================
# SAVE FOREST PLOT
# ============================================================

ggsave(
  filename = file.path(
    out_dir,
    "Figure_TCRbeta_combined_sensitivity.pdf"
  ),
  plot = fig_trb,
  width = 14,
  height = 6.5
)
