# ============================================================
# AALEN-JOHANSEN CURVE
#
# Cumulative incidence of recurrence
# Low vs high sjTREC
#
# Competing event: death without recurrence
# 95% CI + Gray's test p-value
# ============================================================

library(tidyverse)
library(survival)
library(tidycmprsk)
library(ggsurvfit)

# ============================================================
# 1. PATHS + COLOURS
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

col_low <- "#B55F72"
col_high <- "#2F6F73"

# ============================================================
# 2. LOAD + PREPARE DATA
# ============================================================

dat <- readRDS(
  file.path(
    data_dir,
    "RCC_sjTREC_TRB.rds"
  )
)

d <- dat %>%
  filter(
    !is.na(sjTREC),
    !is.na(t_rec),
    !is.na(rec),
    !is.na(death),
    !is.na(sjTREC_group)
  ) %>%
  mutate(
    sjTREC_group = factor(
      sjTREC_group,
      levels = c(
        "Low sjTREC",
        "High sjTREC"
      )
    ),
    status_code = case_when(
      rec == 1 ~ 1L,
      rec == 0 & death == 1 ~ 2L,
      TRUE ~ 0L
    ),
    status = factor(
      status_code,
      levels = c(0, 1, 2),
      labels = c(
        "Censored",
        "Recurrence",
        "Death without recurrence"
      )
    )
  )

# ============================================================
# 3. VALIDATE EVENT CODING
# ============================================================

stopifnot(
  all(d$rec %in% c(0, 1)),
  all(d$death %in% c(0, 1)),
  all(is.finite(d$t_rec)),
  all(d$t_rec >= 0),
  !anyNA(d$sjTREC_group)
)



# For deaths without recurrence, t_rec must be
# the time to death, not an earlier censoring time.

competing_deaths <- d %>%
  filter(status_code == 2)

if (nrow(competing_deaths) > 0) {
  stopifnot(
    !anyNA(competing_deaths$t_os),
    all(
      abs(
        competing_deaths$t_rec -
          competing_deaths$t_os
      ) < 0.01
    )
  )
}

cat("\nEVENT COUNTS\n")

print(
  d %>%
    count(sjTREC_group, status)
)


# ============================================================
# 4. AALEN-JOHANSEN CURVE WITH 95% CI + GRAY'S TEST
# ============================================================

# Calculate cumulative incidence
fit_aj <- tidycmprsk::cuminc(
  Surv(t_rec, status) ~ sjTREC_group,
  data = d
)

# Print Gray's test results
print(fit_aj)
p_gray <- tidycmprsk::glance(fit_aj)$p.value_1
p_label <- if (p_gray < 0.001) "p < 0.001" else sprintf("p = %.3f", p_gray)

# Create figure
p_aj <- fit_aj %>%
  ggcuminc(
    outcome = "Recurrence",
    linewidth = 0.8
  ) +
  add_confidence_interval(
    type = "ribbon",
    alpha = 0.10
  ) +
  add_risktable(
    times = seq(0, 12.5, 2.5),
    risktable_stats = "n.risk",
    stats_label = list(
      n.risk = "Number at risk"
    ),
    size = 3.3
  ) +
  scale_colour_manual(
    values = c(
      "Low sjTREC" = col_low,
      "High sjTREC" = col_high
    )
  ) +
  scale_fill_manual(
    values = c(
      "Low sjTREC" = col_low,
      "High sjTREC" = col_high
    )
  ) +
  scale_x_continuous(
    breaks = seq(0, 12.5, 2.5),
    limits = c(0, 12.5)
  ) +
  scale_y_continuous(
    labels = scales::label_percent(),
    breaks = seq(0, 0.6, 0.1),
    limits = c(0, 0.6)
  ) +
  labs(
    title = "A   Cumulative incidence of recurrence",
    subtitle = paste("Gray's test", p_label),
    x = "Years",
    y = "Cumulative incidence of recurrence (%)",
    colour = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    plot.title = element_text(
      size = 10,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 8,
      colour = "grey40"
    ),
    axis.title = element_text(size = 8.5),
    axis.text = element_text(size = 7.5),
    legend.text = element_text(size = 7)
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.76, 0.90)
  )

print(p_aj)

# ============================================================
# 5. SAVE
# ============================================================

ggsave(
  file.path(
    out_dir,
    "Aalen_Johansen_sjTREC.pdf"
  ),
  p_aj,
  width = 140,
  height = 140,
  units = "mm"
)

ggsave(
  file.path(
    out_dir,
    "Aalen_Johansen_sjTREC.png"
  ),
  p_aj,
  width = 140,
  height = 140,
  units = "mm",
  dpi = 600
)
