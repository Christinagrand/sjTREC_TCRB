# ============================================================
# SUPPLEMENTARY FIGURE
# FULL ADJUSTED COX MODELS
#
# A: sjTREC - recurrence-free survival
# B: sjTREC - overall survival
# C: TCRβ - recurrence-free survival
# D: TCRβ - overall survival
#
# All model covariates shown
# ============================================================

library(tidyverse)
library(survival)
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
# 3. sjTREC DATA
# ============================================================

trec_data <- dat %>%
  filter(
    !is.na(sjTREC)
  )

trec_os <- trec_data %>%
  filter(
    !is.na(t_os),
    !is.na(death)
  )

trec_rfs <- trec_data %>%
  filter(
    !is.na(t_rec),
    !is.na(rec)
  )


# ============================================================
# 4. TCRβ DATA
# ============================================================

trb_data <- dat %>%
  filter(
    !is.na(TRB_z)
  ) %>%
  droplevels()

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
# 5. COX MODELS
# ============================================================

cox_trec_os <- coxph(
  Surv(t_os, death) ~
    sjTREC_group + Age + Sex + Leibovich,
  data = trec_os
)

cox_trec_rfs <- coxph(
  Surv(t_rec, rec) ~
    sjTREC_group + Age + Sex + Leibovich,
  data = trec_rfs
)

cox_trb_os <- coxph(
  Surv(t_os, death) ~
    TRB_z + Age + Sex + Leibovich,
  data = trb_os
)

cox_trb_rfs <- coxph(
  Surv(t_rec, rec) ~
    TRB_z + Age + Sex + Leibovich,
  data = trb_rfs
)


# ============================================================
# 6. EXTRACT COX RESULTS
# ============================================================

extract_cox <- function(fit) {
  s <- summary(fit)

  coef_df <- as.data.frame(
    s$coefficients
  )

  ci_df <- as.data.frame(
    s$conf.int
  )

  tibble(
    term = rownames(coef_df),
    HR = coef_df[["exp(coef)"]],
    low = ci_df[["lower .95"]],
    high = ci_df[["upper .95"]],
    p = coef_df[["Pr(>|z|)"]]
  )
}


# ============================================================
# 7. FOREST-PLOT FUNCTION
# ============================================================

make_forest <- function(
  fit,
  title,
  subtitle,
  term_labels,
  colour
) {
  d <- extract_cox(fit) %>%
    mutate(
      label = recode(
        term,
        !!!term_labels
      ),
      result = if_else(
        p < 0.001,
        sprintf(
          "%.2f (%.2f–%.2f), p < 0.001",
          HR,
          low,
          high
        ),
        sprintf(
          "%.2f (%.2f–%.2f), p = %.3f",
          HR,
          low,
          high,
          p
        )
      ),
      y = rev(
        seq_len(
          n()
        )
      )
    )

  forest <- ggplot(
    d,
    aes(
      x = HR,
      y = y
    )
  ) +
    geom_vline(
      xintercept = 1,
      linetype = "dashed",
      colour = "grey60",
      linewidth = 0.4
    ) +
    geom_segment(
      aes(
        x = low,
        xend = high,
        yend = y
      ),
      linewidth = 0.65,
      colour = colour
    ) +
    geom_point(
      size = 2.4,
      colour = colour
    ) +
    scale_x_log10() +
    scale_y_continuous(
      breaks = d$y,
      labels = d$label,
      limits = c(
        0.5,
        max(d$y) + 0.5
      )
    ) +
    labs(
      x = "Hazard ratio",
      y = NULL
    ) +
    theme_classic(
      base_size = 9
    ) +
    theme(
      axis.text.y = element_text(
        size = 7.5,
        colour = "black"
      ),
      axis.text.x = element_text(
        size = 7,
        colour = "black"
      ),
      axis.title.x = element_text(
        size = 8
      ),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(
        2,
        2,
        2,
        2
      )
    )

  text <- ggplot(
    d,
    aes(
      x = 0,
      y = y
    )
  ) +
    geom_text(
      aes(
        label = result
      ),
      hjust = 0,
      size = 2.5
    ) +
    annotate(
      "text",
      x = 0,
      y = max(d$y) + 0.35,
      label = "HR (95% CI), p value",
      hjust = 0,
      fontface = "bold",
      size = 2.6
    ) +
    scale_x_continuous(
      limits = c(
        0,
        1
      )
    ) +
    scale_y_continuous(
      limits = c(
        0.5,
        max(d$y) + 0.5
      )
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_void()

  inside <- wrap_plots(
    forest,
    text,
    nrow = 1,
    widths = c(
      1.35,
      1
    )
  ) +
    plot_annotation(
      title = title,
      subtitle = subtitle,
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 10
        ),
        plot.subtitle = element_text(
          size = 7.5,
          colour = "grey35"
        )
      )
    )

  wrap_elements(
    full = inside
  )
}


# ============================================================
# 8. TERM LABELS
# ============================================================

trec_labels <- c(
  "sjTREC_groupHigh sjTREC" = "High vs low sjTREC",
  "Age" = "Age (per year)",
  "SexWoman" = "Woman vs man",
  "Leibovichintermediate" = "Leibovich: intermediate vs low",
  "Leibovichhigh" = "Leibovich: high vs low"
)

trb_labels <- c(
  "TRB_z" = "TCRβ diversity (per 1 SD)",
  "Age" = "Age (per year)",
  "SexWoman" = "Woman vs man",
  "Leibovichintermediate" = "Leibovich: intermediate vs low",
  "Leibovichhigh" = "Leibovich: high vs low"
)


# ============================================================
# 9. PANELS
# ============================================================

panel_A <- make_forest(
  cox_trec_rfs,
  "A   sjTREC: recurrence-free survival",
  paste0(
    "N = ", nrow(trec_rfs),
    ", events = ", sum(trec_rfs$rec)
  ),
  trec_labels,
  "#2F6F73"
)

panel_B <- make_forest(
  cox_trec_os,
  "B   sjTREC: overall survival",
  paste0(
    "N = ", nrow(trec_os),
    ", events = ", sum(trec_os$death)
  ),
  trec_labels,
  "#2F6F73"
)

panel_C <- make_forest(
  cox_trb_rfs,
  "C   TCRβ: recurrence-free survival",
  paste0(
    "N = ",
    nrow(trb_rfs),
    ", events = ",
    sum(trb_rfs$rec)
  ),
  trb_labels,
  "#6E8FB3"
)

panel_D <- make_forest(
  cox_trb_os,
  "D   TCRβ: overall survival",
  paste0(
    "N = ",
    nrow(trb_os),
    ", events = ",
    sum(trb_os$death)
  ),
  trb_labels,
  "#6E8FB3"
)


# ============================================================
# 10. COMBINE 2 x 2
# ============================================================

figure_supp_cox <- wrap_plots(
  panel_A,
  panel_B,
  panel_C,
  panel_D,
  ncol = 2
)

figure_supp_cox


# ============================================================
# 11. SAVE
# ============================================================

ggsave(
  file.path(
    out_dir,
    "Supplementary_Figure_full_Cox_models.pdf"
  ),
  figure_supp_cox,
  width = 240,
  height = 160,
  units = "mm",
  device = cairo_pdf
)


# ============================================================
# 12. MULTIVARIABLE COX WITH BOTH sjTREC + TCRβ
#
# A: recurrence-free survival
# B: overall survival
#
# sjTREC = median split
# TCRβ = continuous, per 1 SD
# Adjusted for age, sex and Leibovich
# ============================================================

both_data <- dat %>%
  filter(
    !is.na(sjTREC),
    !is.na(TRB_z)
  )

both_os <- both_data %>%
  filter(
    !is.na(t_os),
    !is.na(death)
  )

both_rfs <- both_data %>%
  filter(
    !is.na(t_rec),
    !is.na(rec)
  )


# ============================================================
# 12.1 COX MODELS
# ============================================================

cox_both_os <- coxph(
  Surv(t_os, death) ~
    sjTREC_group +
    TRB_z +
    Age +
    Sex +
    Leibovich,
  data = both_os
)

cox_both_rfs <- coxph(
  Surv(t_rec, rec) ~
    sjTREC_group +
    TRB_z +
    Age +
    Sex +
    Leibovich,
  data = both_rfs
)


# ============================================================
# 12.2 TERM LABELS
# ============================================================

both_labels <- c(
  "sjTREC_groupHigh sjTREC" = "High vs low sjTREC",
  "TRB_z" = "TCRβ diversity (per 1 SD)",
  "Age" = "Age (per year)",
  "SexWoman" = "Woman vs man",
  "Leibovichintermediate" = "Leibovich: intermediate vs low",
  "Leibovichhigh" = "Leibovich: high vs low"
)


# ============================================================
# 12.3 FOREST PLOTS
# ============================================================

panel_both_A <- make_forest(
  cox_both_rfs,
  "A   sjTREC + TCRβ: recurrence-free survival",
  paste0(
    "N = ", nrow(both_rfs),
    ", events = ", sum(both_rfs$rec)
  ),
  both_labels,
  "#6E8FB3"
)

panel_both_B <- make_forest(
  cox_both_os,
  "B   sjTREC + TCRβ: overall survival",
  paste0(
    "N = ", nrow(both_os),
    ", events = ", sum(both_os$death)
  ),
  both_labels,
  "#6E8FB3"
)


# ============================================================
# 12.4 COMBINE
# ============================================================

figure_both_cox <- wrap_plots(
  panel_both_A,
  panel_both_B,
  ncol = 2
)

figure_both_cox


# ============================================================
# 12.5 SAVE
# ============================================================

ggsave(
  file.path(
    out_dir,
    "Supplementary_Figure_sjTREC_TRB_multivariable_Cox.pdf"
  ),
  figure_both_cox,
  width = 260,
  height = 85,
  units = "mm",
  device = cairo_pdf
)
