# ============================================================
# SUPPLEMENTARY OVERVIEW FIGURE
#
# A: sjTREC distribution
# B: TCRβ diversity distribution
# C: TCRβ diversity by sjTREC group
# D: sjTREC vs TCRβ diversity
# E: sjTREC vs age
# F: TCRβ diversity vs age
#
# Correlations assessed using Spearman's rank correlation
# ============================================================

library(tidyverse)
library(patchwork)


# ============================================================
# PATHS
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
# COLOURS
# ============================================================

col_dark <- "#7A5490"
col_light <- "#C9AED6"
col_grey <- "grey80"
col_point <- "grey35"


# ============================================================
# LOAD DATA
# ============================================================

dat <- readRDS(
  file.path(
    data_dir,
    "RCC_sjTREC_TRB.rds"
  )
)

# ============================================================
# TRANSFORM TRB SIMPSON FOR VISUALISATION
# ============================================================

trec_offset <- min(
  dat$sjTREC[dat$sjTREC > 0],
  na.rm = TRUE
) / 2

dat <- dat %>%
  mutate(
    TRB_plot = -log10(
      1 - TRB
    )
  )


# ============================================================
# AXIS BREAKS
#
# Coordinates remain transformed,
# but labels show values on the original scale.
# ============================================================

# sjTREC values we want displayed
sjtrec_breaks_original <- c(
  0,
  0.00001,
  0.0001,
  0.0005,
  0.0015
)

sjtrec_breaks_transformed <- log2(
  sjtrec_breaks_original + trec_offset
)

sjtrec_labels <- c(
  "0",
  "0.00001",
  "0.0001",
  "0.0005",
  "0.0015"
)


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
# SPEARMAN CORRELATION
# ============================================================

cor_dat <- dat %>%
  filter(
    !is.na(TRB),
    !is.na(log2_sjTREC),
    !is.na(TRB_plot)
  )

cor_res <- cor.test(
  cor_dat$TRB,
  cor_dat$sjTREC,
  method = "spearman"
)

rho_txt <- sprintf(
  "%.2f",
  unname(cor_res$estimate)
)

p_cor_txt <- ifelse(
  cor_res$p.value < 0.001,
  "p < 0.001",
  paste0(
    "p = ",
    sprintf(
      "%.3f",
      cor_res$p.value
    )
  )
)

cor_label <- paste0(
  "Spearman \u03C1 = ",
  rho_txt,
  "\n",
  p_cor_txt
)


# ============================================================
# WILCOXON: TRB BY sjTREC GROUP
# ============================================================

wilcox_res <- wilcox.test(
  TRB ~ sjTREC_group,
  data = dat,
  exact = FALSE
)

p_group_txt <- ifelse(
  wilcox_res$p.value < 0.001,
  "p < 0.001",
  paste0(
    "p = ",
    sprintf(
      "%.3f",
      wilcox_res$p.value
    )
  )
)


# ============================================================
# COMMON THEME
# ============================================================

common_theme <- theme_classic(
  base_size = 11
) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    axis.title = element_text(
      size = 11
    ),
    axis.text = element_text(
      size = 10,
      colour = "black"
    )
  )


# ============================================================
# PANEL A
# sjTREC DISTRIBUTION
# ============================================================


panel_A <- ggplot(
  dat %>%
    filter(
      !is.na(log2_sjTREC)
    ),
  aes(
    x = log2_sjTREC
  )
) +
  geom_histogram(
    aes(
      y = after_stat(density)
    ),
    bins = 30,
    fill = col_grey,
    colour = NA
  ) +
  geom_density(
    linewidth = 0.9,
    colour = "black"
  ) +
  scale_x_continuous(
    breaks = sjtrec_breaks_transformed,
    labels = sjtrec_labels
  ) +
  labs(
    title = "A   sjTREC levels",
    x = "sjTREC per cell",
    y = "Density"
  ) +
  common_theme


# ============================================================
# PANEL B
# TRB SIMPSON DISTRIBUTION
# ============================================================

panel_B <- ggplot(
  dat %>%
    filter(
      !is.na(TRB)
    ),
  aes(
    x = TRB
  )
) +
  geom_histogram(
    bins = 30,
    fill = col_grey,
    colour = NA
  ) +
  scale_x_continuous(
    breaks = c(
      0.65,
      0.75,
      0.85,
      0.90,
      0.95,
      0.99
    )
  ) +
  labs(
    title = "B   TCRβ diversity",
    x = "TCRβ Simpson diversity",
    y = "Count"
  ) +
  common_theme


# ============================================================
# PANEL C
# TRB DIVERSITY BY sjTREC GROUP
# ============================================================

panel_C <- ggplot(
  dat %>%
    filter(
      !is.na(TRB_plot),
      !is.na(sjTREC_group)
    ),
  aes(
    x = sjTREC_group,
    y = TRB_plot,
    fill = sjTREC_group
  )
) +
  geom_boxplot(
    width = 0.52,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  geom_jitter(
    width = 0.11,
    alpha = 0.40,
    size = 1.4,
    colour = "grey30"
  ) +
  scale_fill_manual(
    values = c(
      "Low sjTREC" = col_light,
      "High sjTREC" = col_dark
    )
  ) +
  scale_y_continuous(
    breaks = trb_breaks_transformed,
    labels = trb_labels
  ) +
  annotate(
    "text",
    x = 1.5,
    y = Inf,
    label = p_group_txt,
    vjust = 1.5,
    size = 3.5
  ) +
  labs(
    title = "C   TCRβ diversity by sjTREC",
    x = NULL,
    y = "TCRβ Simpson diversity"
  ) +
  common_theme +
  theme(
    legend.position = "none"
  )


# ============================================================
# PANEL D
# sjTREC vs TRB DIVERSITY
# ============================================================

panel_D <- ggplot(
  cor_dat,
  aes(
    x = TRB_plot,
    y = log2_sjTREC
  )
) +
  geom_point(
    size = 1.6,
    alpha = 0.50,
    colour = col_point
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.9,
    colour = col_dark,
    fill = "grey85"
  ) +
  scale_x_continuous(
    breaks = trb_breaks_transformed,
    labels = trb_labels
  ) +
  scale_y_continuous(
    breaks = sjtrec_breaks_transformed,
    labels = sjtrec_labels
  ) +
  annotate(
    "text",
    x = min(
      cor_dat$TRB_plot,
      na.rm = TRUE
    ),
    y = max(
      cor_dat$log2_sjTREC,
      na.rm = TRUE
    ),
    label = cor_label,
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  # Panel D
  labs(
    title = "D  sjTREC vs TCRβ diversity",
    x = "TCRβ Simpson diversity",
    y = "sjTREC per cell"
  ) +
  common_theme


# ============================================================
# AGE DATA
# ============================================================

age_trec_dat <- dat %>%
  filter(
    !is.na(Age),
    !is.na(log2_sjTREC)
  )

age_trb_dat <- dat %>%
  filter(
    !is.na(Age),
    !is.na(TRB),
    !is.na(TRB_plot)
  )


# ============================================================
# sjTREC vs AGE
# ============================================================

cor_age_trec <- cor.test(
  age_trec_dat$Age,
  age_trec_dat$sjTREC,
  method = "spearman"
)

rho_age_trec <- sprintf(
  "%.2f",
  unname(cor_age_trec$estimate)
)

p_age_trec <- ifelse(
  cor_age_trec$p.value < 0.001,
  "p < 0.001",
  paste0(
    "p = ",
    sprintf(
      "%.3f",
      cor_age_trec$p.value
    )
  )
)

label_age_trec <- paste0(
  "Spearman \u03C1 = ",
  rho_age_trec,
  "\n",
  p_age_trec
)


# ============================================================
# TRB vs AGE
# ============================================================

cor_age_trb <- cor.test(
  age_trb_dat$Age,
  age_trb_dat$TRB,
  method = "spearman"
)

rho_age_trb <- sprintf(
  "%.2f",
  unname(cor_age_trb$estimate)
)

p_age_trb <- ifelse(
  cor_age_trb$p.value < 0.001,
  "p < 0.001",
  paste0(
    "p = ",
    sprintf(
      "%.3f",
      cor_age_trb$p.value
    )
  )
)

label_age_trb <- paste0(
  "Spearman \u03C1 = ",
  rho_age_trb,
  "\n",
  p_age_trb
)


# ============================================================
# PANEL E
# AGE vs sjTREC
# ============================================================


panel_E <- ggplot(
  age_trec_dat,
  aes(
    x = Age,
    y = log2_sjTREC
  )
) +
  geom_point(
    size = 1.6,
    alpha = 0.50,
    colour = col_point
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.9,
    colour = col_dark,
    fill = "grey85"
  ) +
  scale_y_continuous(
    breaks = sjtrec_breaks_transformed,
    labels = sjtrec_labels
  ) +
  annotate(
    "text",
    x = min(
      age_trec_dat$Age,
      na.rm = TRUE
    ),
    y = max(
      age_trec_dat$log2_sjTREC,
      na.rm = TRUE
    ),
    label = label_age_trec,
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  labs(
    title = "E   sjTREC and age",
    x = "Age (years)",
    y = "sjTREC per cell"
  ) +
  common_theme


# ============================================================
# PANEL F
# AGE vs TRB DIVERSITY
# ============================================================

panel_F <- ggplot(
  age_trb_dat,
  aes(
    x = Age,
    y = TRB_plot
  )
) +
  geom_point(
    size = 1.6,
    alpha = 0.50,
    colour = col_point
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 0.9,
    colour = col_dark,
    fill = "grey85"
  ) +
  scale_y_continuous(
    breaks = trb_breaks_transformed,
    labels = trb_labels
  ) +
  annotate(
    "text",
    x = min(
      age_trb_dat$Age,
      na.rm = TRUE
    ),
    y = max(
      age_trb_dat$TRB_plot,
      na.rm = TRUE
    ),
    label = label_age_trb,
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  labs(
    title = "F   TCRβ diversity and age",
    x = "Age (years)",
    y = "TCRβ Simpson diversity"
  ) +
  common_theme


# ============================================================
# COMBINE FIGURE
# ============================================================

overview_plot <-
  (panel_A | panel_B) /
    (panel_C | panel_D) /
    (panel_E | panel_F)

overview_plot


# ============================================================
# SAVE
# ============================================================

ggsave(
  filename = file.path(
    out_dir,
    "FigureS_overview_sjTREC_TCRbeta.pdf"
  ),
  plot = overview_plot,
  width = 180,
  height = 200,
  units = "mm",
  device = cairo_pdf
)
