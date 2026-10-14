# ============================================================
# TCRβ SEQUENCING KIT COMPARISON
#
# A: Kit 1 vs Kit 2 Gini-Simpson diversity
# B: Kit 1 vs Kit 2 normalized Shannon diversity
# C: Kit 1 vs Kit 2 Gini-Simpson ECDF
#
# Kit differences assessed using Wilcoxon rank-sum tests
# ============================================================

library(tidyverse)
library(patchwork)

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


# ============================================================
# TRB AXIS BREAKS
#
# Coordinates remain transformed,
# but labels show values on the original Simpson scale.
# ============================================================

trb_breaks_original <- c(
  0.65,
  0.90,
  0.99,
  0.999,
  1
)

trb_breaks_transformed <- -log10(
  1 - trb_breaks_original
)

trb_labels <- c(
  "0.65",
  "0.90",
  "0.99",
  "0.999",
  "1"
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
# 3. PREPARE VARIABLES
# ============================================================

dat <- dat %>%
  mutate(
    cohort = factor(
      cohort,
      levels = c(
        "Kit 1",
        "Kit 2"
      )
    ),

    # Transform TRB Simpson for visualisation
    TRB_plot = -log10(
      1 - TRB
    )
  )


# ============================================================
# 4. TRB DISTRIBUTION DATA
# ============================================================

trb_dist <- dat %>%
  filter(
    !is.na(TRB),
    !is.na(normalized_shannon),
    !is.na(cohort)
  )


# ============================================================
# 5. DISTRIBUTION TESTS
# ============================================================

wilcox_simpson <- wilcox.test(
  TRB ~ cohort,
  data = trb_dist,
  exact = FALSE
)


wilcox_shannon <- wilcox.test(
  normalized_shannon ~ cohort,
  data = trb_dist,
  exact = FALSE
)
# ============================================================
# 6. HELPER FOR P-VALUE LABELS
# ============================================================

format_p <- function(p) {
  if (p < 0.001) {
    "p < 0.001"
  } else {
    paste0(
      "p = ",
      formatC(
        p,
        format = "f",
        digits = 3
      )
    )
  }
}

# ============================================================
# FIGURE S1
# SEQUENCING KIT COMPARISON
# ============================================================


# ============================================================
# S1A. SIMPSON
# ============================================================

p_s1a <- ggplot(
  trb_dist,
  aes(
    x = cohort,
    y = TRB_plot,
    fill = cohort
  )
) +
  geom_violin(
    width = 0.8,
    trim = FALSE,
    alpha = 0.35,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    alpha = 0.45,
    size = 1.3
  ) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = format_p(
      wilcox_simpson$p.value
    ),
    vjust = 1.5,
    size = 3.5
  ) +
  scale_y_continuous(
    breaks = trb_breaks_transformed,
    labels = trb_labels
  ) +
  labs(
    x = NULL,
    y = "TCRβ Simpson diversity"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none"
  )


# ============================================================
# S1B. NORMALIZED SHANNON
# ============================================================

p_s1b <- ggplot(
  trb_dist,
  aes(
    x = cohort,
    y = normalized_shannon,
    fill = cohort
  )
) +
  geom_violin(
    width = 0.8,
    trim = FALSE,
    alpha = 0.35,
    colour = NA
  ) +
  geom_boxplot(
    width = 0.22,
    outlier.shape = NA,
    alpha = 0.65
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    alpha = 0.45,
    size = 1.3
  ) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = format_p(
      wilcox_shannon$p.value
    ),
    vjust = 1.5,
    size = 3.5
  ) +
  labs(
    x = NULL,
    y = "TCRβ normalized Shannon diversity"
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "none"
  )


# ============================================================
# S1C. SIMPSON ECDF
# ============================================================

p_s1c <- ggplot(
  trb_dist,
  aes(
    x = TRB_plot,
    colour = cohort
  )
) +
  stat_ecdf(
    linewidth = 0.9
  ) +
  scale_x_continuous(
    breaks = trb_breaks_transformed,
    labels = trb_labels
  ) +
  labs(
    x = "TCRβ Simpson diversity",
    y = "Cumulative proportion",
    colour = NULL
  ) +
  theme_classic(
    base_size = 11
  ) +
  theme(
    legend.position = "top"
  )


# ============================================================
# COMBINE FIGURE S1
# ============================================================

fig_s1 <- (
  p_s1a |
    p_s1b |
    p_s1c
) +
  plot_annotation(
    tag_levels = "A"
  )


ggsave(
  filename = file.path(
    out_dir,
    "FigureS_TCRbeta_kit_comparison.pdf"
  ),
  plot = fig_s1,
  width = 10.5,
  height = 3.8
)
