# ============================================================
# ccRCC SENSITIVITY ANALYSIS
#
# A: sjTREC high vs low - recurrence-free survival
# B: sjTREC high vs low - overall survival
# C: adjusted sjTREC Cox - RFS + OS
# D: adjusted TCRβ Cox - RFS + OS
#
# Clear-cell RCC only
# TCRβ diversity modelled per 1 SD
# Cox models adjusted for age, sex, and Leibovich risk group
# ============================================================

library(tidyverse)
library(survival)
library(survminer)
library(broom)
library(patchwork)
library(grid)


# ============================================================
# 1. PATHS
# ============================================================

data_dir <- "data"
out_dir <- file.path(
  "results",
  "figures"
)

dir.create(
  out_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# sjTREC
col_high <- "#2F6F73"
col_low <- "#B55F72"

# TRB
col_trb <- "#6E8FB3"


# ============================================================
# 2. LOAD + PREPARE
# ============================================================

dat <- readRDS(
  file.path(
    data_dir,
    "RCC_sjTREC_TRB.rds"
  )
)

# Restrict to clear-cell RCC
dat <- dat %>%
  filter(
    subtype_PatoMcode == "ccRCC"
  )
# ============================================================
# 3. P VALUE FUNCTION
# ============================================================

format_p <- function(p) {
  ifelse(
    p < 0.001,
    "p < 0.001",
    sprintf(
      "p = %.3f",
      p
    )
  )
}


# ============================================================
# 4. sjTREC OS DATA
# ============================================================


os_data <- dat %>%
  filter(
    !is.na(sjTREC),
    !is.na(t_os),
    !is.na(death)
  )


# ============================================================
# 5. sjTREC RECURRENCE-FREE SURVIVAL DATA
# ============================================================

rfs_data <- dat %>%
  filter(
    !is.na(sjTREC),
    !is.na(t_rec),
    !is.na(rec)
  )

# ============================================================
# 6. KM PANEL FUNCTION
#    Log-rank p + unadjusted Cox HR
# ============================================================

make_km_panel <- function(
  data,
  time,
  event,
  group,
  panel_letter,
  title,
  ylab,
  palette,
  legend_labs
) {
  d <- data %>%
    mutate(
      .time = .data[[time]],
      .event = .data[[event]],
      .group = .data[[group]]
    )


  # ----------------------------------------------------------
  # Kaplan-Meier fit
  # ----------------------------------------------------------

  fit <- survfit(
    Surv(
      .time,
      .event
    ) ~ .group,
    data = d
  )


  # ----------------------------------------------------------
  # Log-rank test
  # ----------------------------------------------------------

  p <- survdiff(
    Surv(
      .time,
      .event
    ) ~ .group,
    data = d
  )$pvalue


  # ----------------------------------------------------------
  # Unadjusted Cox regression
  #
  # Low sjTREC = reference
  # HR therefore = High vs Low
  # ----------------------------------------------------------

  cox_unadjusted <- coxph(
    Surv(
      .time,
      .event
    ) ~ .group,
    data = d
  )


  cox_result <- tidy(
    cox_unadjusted,
    exponentiate = TRUE,
    conf.int = TRUE
  )


  hr_text <- sprintf(
    "HR = %.2f (95%% CI %.2f–%.2f)",
    cox_result$estimate[1],
    cox_result$conf.low[1],
    cox_result$conf.high[1]
  )


  # ----------------------------------------------------------
  # KM plot
  # ----------------------------------------------------------

  g <- ggsurvplot(
    fit,
    data = d,
    fun = "pct",
    conf.int = TRUE,
    conf.int.alpha = 0.06,
    risk.table = TRUE,
    risk.table.height = 0.18,
    risk.table.title = NULL,
    risk.table.fontsize = 2.5,
    palette = palette,
    legend.title = "",
    legend.labs = legend_labs,
    xlab = "Years",
    ylab = ylab,
    break.x.by = 2.5,
    censor.size = 1.5,
    size = 0.8,
    ggtheme = theme_classic(
      base_size = 9
    )
  )


  # ----------------------------------------------------------
  # Survival plot
  # ----------------------------------------------------------

  g$plot <- g$plot +
    labs(
      title = paste0(
        panel_letter,
        "   ",
        title
      ),
      subtitle = paste0(
        "Log-rank ", format_p(p),
        "\n",
        hr_text
      )
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 10
      ),
      plot.subtitle = element_text(
        size = 8,
        colour = "grey30",
        lineheight = 1.1
      ),
      axis.title = element_text(
        size = 8.5
      ),
      axis.text = element_text(
        size = 7.5,
        colour = "black"
      ),
      legend.position = c(
        0.76,
        0.94
      ),
      legend.direction = "vertical",
      legend.text = element_text(
        size = 7
      ),
      legend.key.width = unit(
        0.8,
        "cm"
      ),
      legend.key.height = unit(
        0.30,
        "cm"
      ),
      legend.background = element_rect(
        fill = scales::alpha(
          "white",
          0.75
        ),
        colour = NA
      ),
      plot.margin = margin(
        3,
        3,
        0,
        3
      )
    )


  # ----------------------------------------------------------
  # Risk table
  # ----------------------------------------------------------

  g$table <- g$table +
    theme_classic(
      base_size = 7
    ) +
    theme(
      axis.text.x = element_text(
        size = 6.7,
        colour = "black"
      ),
      axis.text.y = element_text(
        size = 6.5,
        colour = "black"
      ),
      axis.title.x = element_text(
        size = 7.5
      ),
      axis.title.y = element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      plot.margin = margin(
        0,
        3,
        2,
        3
      )
    )


  # ----------------------------------------------------------
  # Combine curve + risk table
  # ----------------------------------------------------------

  wrap_elements(
    full = wrap_plots(
      g$plot,
      g$table,
      ncol = 1,
      heights = c(
        4.4,
        1
      )
    )
  )
}


# ============================================================
# 7. PANELS A + B
# ============================================================

panel_B <- make_km_panel(
  os_data,
  "t_os",
  "death",
  "sjTREC_group",
  "B",
  "Overall survival",
  "Overall survival (%)",

  # Factor order is Low, High
  c(
    col_low,
    col_high
  ),
  c(
    "Low sjTREC",
    "High sjTREC"
  )
)


panel_A <- make_km_panel(
  rfs_data,
  "t_rec",
  "rec",
  "sjTREC_group",
  "A",
  "Recurrence-free survival",
  "Recurrence-free survival (%)",

  # Factor order is Low, High
  c(
    col_low,
    col_high
  ),
  c(
    "Low sjTREC",
    "High sjTREC"
  )
)


# ============================================================
# 8. ADJUSTED sjTREC COX MODELS
# ============================================================

cox_os <- coxph(
  Surv(
    t_os,
    death
  ) ~
    sjTREC_group +
    Age +
    Sex +
    Leibovich,
  data = os_data
)


cox_rfs <- coxph(
  Surv(
    t_rec,
    rec
  ) ~
    sjTREC_group +
    Age +
    Sex +
    Leibovich,
  data = rfs_data
)


# ============================================================
# 9. PANEL C: sjTREC COX
# ============================================================

forest_dat <- bind_rows(
  tidy(
    cox_rfs,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term == "sjTREC_groupHigh sjTREC"
    ) %>%
    mutate(
      Outcome = "Recurrence-free survival"
    ),
  tidy(
    cox_os,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term == "sjTREC_groupHigh sjTREC"
    ) %>%
    mutate(
      Outcome = "Overall survival"
    )
) %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = c(
        "Overall survival",
        "Recurrence-free survival"
      )
    ),
    hr_text = sprintf(
      "%.2f (%.2f–%.2f)",
      estimate,
      conf.low,
      conf.high
    ),
    p_text = if_else(
      p.value < 0.001,
      "p < 0.001",
      sprintf(
        "p = %.3f",
        p.value
      )
    ),

    # Put HR/CI and p-value on separate lines
    stat_text = paste0(
      hr_text,
      "\n",
      p_text
    )
  )

# ------------------------------------------------------------
# Panel C: forest part
# ------------------------------------------------------------

forest_plot <- ggplot(
  forest_dat,
  aes(
    x = estimate,
    y = Outcome
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey60"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.08,
    linewidth = 0.7,
    colour = col_high
  ) +
  geom_point(
    size = 2.6,
    colour = col_high
  ) +
  scale_x_log10() +
  labs(
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 7.5,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 8.5
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(
      2,
      3,
      2,
      2
    )
  )


# ------------------------------------------------------------
# Panel C: statistics text
# ------------------------------------------------------------

forest_text <- ggplot(
  forest_dat,
  aes(
    x = 0,
    y = Outcome
  )
) +
  geom_text(
    aes(
      label = stat_text
    ),
    hjust = 0,
    vjust = 0.5,
    lineheight = 1.15,
    size = 2.7
  ) +
  annotate(
    "text",
    x = 0,
    y = 2.45,
    label = "HR (95% CI)",
    fontface = "bold",
    hjust = 0,
    size = 2.8
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    expand = c(
      0,
      0
    )
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      2,
      3,
      2,
      2
    )
  )


# ------------------------------------------------------------
# Combine panel C
# ------------------------------------------------------------

panel_C_inside <- wrap_plots(
  forest_plot,
  forest_text,
  nrow = 1,
  widths = c(
    1.55,
    0.95
  )
) +

  plot_annotation(
    title = "C   sjTREC Cox regression",
    subtitle = paste0(
      "High vs low sjTREC\n",
      "Adjusted for age, sex, and Leibovich"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 10
      ),
      plot.subtitle = element_text(
        size = 8,
        colour = "grey35"
      )
    )
  )


panel_C <- wrap_elements(
  full = panel_C_inside
)


# ============================================================
# 10. TRB DATA
# ============================================================

trb_data <- dat %>%
  filter(!is.na(TRB))


# ============================================================
# 11. TRB OS + RECURRENCE-FREE SURVIVAL ANALYSIS DATA
# ============================================================

trb_os <- trb_data %>%
  filter(
    !is.na(t_os),
    !is.na(death)
  )


trb_rfs <- trb_data %>%
  filter(
    !is.na(t_rec),
    !is.na(rec)
  )

# ============================================================
# 12. CONTINUOUS TRB COX MODELS
# ============================================================

cox_trb_os <- coxph(
  Surv(
    t_os,
    death
  ) ~
    TRB_z +
    Age +
    Sex +
    Leibovich,
  data = trb_os
)


cox_trb_rfs <- coxph(
  Surv(
    t_rec,
    rec
  ) ~
    TRB_z +
    Age +
    Sex +
    Leibovich,
  data = trb_rfs
)

# ============================================================
# 13. PANEL D: CONTINUOUS TRB COX
# ============================================================

trb_forest_dat <- bind_rows(
  tidy(
    cox_trb_os,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term == "TRB_z"
    ) %>%
    mutate(
      Outcome = "Overall survival",
      N = nrow(trb_os),
      Events = sum(trb_os$death)
    ),
  tidy(
    cox_trb_rfs,
    exponentiate = TRUE,
    conf.int = TRUE
  ) %>%
    filter(
      term == "TRB_z"
    ) %>%
    mutate(
      Outcome = "Recurrence-free survival",
      N = nrow(trb_rfs),
      Events = sum(trb_rfs$rec)
    )
) %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = c(
        "Overall survival",
        "Recurrence-free survival"
      )
    ),
    hr_text = sprintf(
      "%.2f (%.2f–%.2f)",
      estimate,
      conf.low,
      conf.high
    ),
    p_text = if_else(
      p.value < 0.001,
      "p < 0.001",
      sprintf(
        "p = %.3f",
        p.value
      )
    ),

    # Put HR/CI and p-value on separate lines
    stat_text = paste0(
      hr_text,
      "\n",
      p_text
    )
  )

# ------------------------------------------------------------
# Panel D: forest part
# ------------------------------------------------------------

trb_forest_plot <- ggplot(
  trb_forest_dat,
  aes(
    x = estimate,
    y = Outcome
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey60"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    height = 0.08,
    linewidth = 0.7,
    colour = col_trb
  ) +
  geom_point(
    size = 2.6,
    colour = col_trb
  ) +
  scale_x_log10() +
  labs(
    x = "Hazard ratio",
    y = NULL
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text.y = element_text(
      size = 8,
      colour = "black"
    ),
    axis.text.x = element_text(
      size = 7.5,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 8.5
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(
      2,
      3,
      2,
      2
    )
  )


# ------------------------------------------------------------
# Panel D: statistics text
# ------------------------------------------------------------

trb_forest_text <- ggplot(
  trb_forest_dat,
  aes(
    x = 0,
    y = Outcome
  )
) +
  geom_text(
    aes(
      label = stat_text
    ),
    hjust = 0,
    vjust = 0.5,
    lineheight = 1.15,
    size = 2.7
  ) +
  annotate(
    "text",
    x = 0,
    y = 2.45,
    label = "HR (95% CI)",
    fontface = "bold",
    hjust = 0,
    size = 2.8
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    expand = c(
      0,
      0
    )
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      2,
      3,
      2,
      2
    )
  )


# ------------------------------------------------------------
# Combine panel D
# ------------------------------------------------------------

panel_D_inside <- wrap_plots(
  trb_forest_plot,
  trb_forest_text,
  nrow = 1,
  widths = c(
    1.55,
    0.95
  )
) +

  plot_annotation(
    title = "D   TCRβ Cox regression",
    subtitle = paste0(
      "TCRβ per 1 SD\n",
      "Adjusted for age, sex, and Leibovich"
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 10
      ),
      plot.subtitle = element_text(
        size = 8,
        colour = "grey35"
      )
    )
  )


panel_D <- wrap_elements(
  full = panel_D_inside
)


# ============================================================
# 14. FINAL ccRCC SENSITIVITY FIGURE
# ============================================================

top_row <- wrap_plots(
  panel_A,
  panel_B,
  nrow = 1,
  widths = c(
    1,
    1
  )
)


bottom_row <- wrap_plots(
  panel_C,
  panel_D,
  nrow = 1,
  widths = c(
    1,
    1
  )
)


figure_main <- wrap_plots(
  top_row,
  bottom_row,
  ncol = 1,
  heights = c(
    2.3,
    1.25
  )
)


figure_main


# ============================================================
# 15. SAVE
# ============================================================

ggsave(
  file.path(
    out_dir,
    "FigureS_ccRCC_sensitivity.pdf"
  ),
  figure_main,
  ...
)
