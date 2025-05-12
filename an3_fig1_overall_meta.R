## Create plot of overall treatment by stage, plus meta-analysis
library(tidyverse)
library(gt)
library(broom)
library(metafor)

#
# Norway - "#e41a1c"
# Australia - "#377eb8"
# Canada - "#4daf4a"
# UK - "#984ea3"
colours <- c("#984ea3", "#e41a1c", "#4daf4a", "#377eb8")

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

# Load in 'raw' data -------------------------------------------------------

# Jdn y-pos
plot_ordering <- function(data) {
  data |>
    mutate(ordering    = as.numeric(ctry_order)) |>
    mutate(ordering    = ifelse(ordering >=  6, ordering+1, ordering)) |>
    mutate(ordering    = ifelse(ordering >= 16, ordering+1, ordering)) 
}

meta_dat <- readRDS("Rdata/dat_ind.RDS") |>
  filter(
      (
        stage %in% c("All stages", "Stages 1-3", "Stage 4")
        &
        variable == "all"
      ) |
      (
        stage == "All stages" 
        &
        variable == "stage"
      ) 
    ) |>
  mutate(stage = 
           case_when(
             stage_x %in% c("1") ~ "1",
             stage_x %in% c("2") ~ "2",
             stage_x %in% c("3") ~ "3",
             stage_x %in% c("4") ~ "4",
             stage_x %in% c("L") ~ "L",
             stage_x %in% c("R") ~ "R",
             stage_x %in% c("D") ~ "D",
             stage_x %in% c("X", "X SEER") ~ "Missing",
             TRUE ~ stage
           )) |>
  mutate(stage = factor(stage, levels = c("All stages", "Stages 1-3", "Stage 4", "Missing", "1", "2", "3", "4", "L", "R", "D"))) |> 
  mutate(stage = factor(stage, labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant", "No recorded stage", "Stage 1", "Stage 2", "Stage 3", "Stage 4", "Localised", "Regional", "Distant"))) |>
  mutate(trt = case_when(
    trt == "Chemo" ~ "Chemotherapy",
    trt == "Radio" ~ "Radiotherapy"
  )) |> 
  plot_ordering()

load("Rdata/useful_vectors.Rdata")

# Create the meta-analysis table --------------------------  ----------------
table <- meta_dat |>
  filter(stage != "No recorded stage") |>
  filter(stage != "Stage 4") |>
  filter(!(stage %in% c("Localised", "Regional", "Distant"))) |>
  group_by(trt, stage) |>
  mutate(minimum = min(prop)) |>
  mutate(maximum = max(prop)) |>
  group_by(trt, stage, minimum, maximum) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop = estimate, lower = conf.low, upper = conf.high, se = std.error))
  ) |>
  mutate(metan_summ = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      glance() |>
      select(i2 = i.squared, tau2 = tau.squared, nobs = nobs)
  )) |>
  select(trt, stage, minimum, maximum, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  mutate(
    pred_lb = expit(prop-qt(.975, df = nobs)*sqrt(tau2+(se^2))),
    pred_ub = expit(prop+qt(.975, df = nobs)*sqrt(tau2+(se^2)))
  ) |>
  mutate(
    prop  = expit(prop),
    lower = expit(lower),
    upper = expit(upper),
    tau = sqrt(tau2)
  ) |> 
  select(trt, stage, prop, lower, upper, pred_lb, pred_ub, tau, i2, minimum, maximum) |>
  group_by() |>
  arrange(trt, stage)
table

overall_analysis1_table <- table |>
  gt() |>
  gt_metan_common() |>
  cols_align(align = "left", columns = c(trt, stage)) |>
  fmt_percent(
    columns = c(prop, lower, upper, pred_lb, pred_ub, minimum, maximum),
    decimals = 1
  ) |>
  cols_label (
    trt = "Treatment",
    stage = "Stage",
    prop = "Pooled estimate %",
    lower = "(95% confidence interval)",
    pred_lb = "95% prediction interval *",
    minimum = "Observed jurisdictional range",
    tau = md("*Tau*, log-odds scale **"),
    i2 = md("*I<sup>2</sup>* ***")
  ) 

overall_analysis1_table |> saveRDS("results/table3_overall_meta.RDS")

# Meta-analysis data processing for plot ------------------------------------------
poly_metan <- 
  meta_dat |>
  filter(stage != "No recorded stage") |>
  group_by(trt, stage) |>
  nest() |>
  mutate(metan = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      tidy(conf.int = TRUE) |>
      select(prop0 = conf.low, prop1 = estimate, prop2 = conf.high, prop3 = estimate, se = std.error))
  ) |>
  mutate(metan_summ = map(
    data, ~ rma(xi = n_trt, ni = n, data = .x, measure = "PLO", method = "REML") |>
      glance() |>
      select(tau2 = tau.squared, nobs = nobs)
  )) |>
  select(trt, stage, metan, metan_summ) |>
  unnest(cols = c(metan, metan_summ)) |>
  group_by(trt) |>
  mutate(
    pred_lb = expit(prop1-qt(.975, df = nobs)*sqrt(tau2+(se^2))),
    pred_ub = expit(prop1+qt(.975, df = nobs)*sqrt(tau2+(se^2)))
  ) |>
  select(trt, stage, pred_lb, pred_ub, prop0, prop1, prop2, prop3) |>
  common_meta_format() |>
  mutate(x = expit(x))

# Draw the forest plot ------------------------------------------

meta_dat_all <- meta_dat
poly_metan_all <- poly_metan

meta_dat   <- meta_dat_all   |> 
  filter(stage %in% c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant", "No recorded stage")) |>
  mutate(stage = factor(
    stage, 
    levels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant", "No recorded stage"),
    labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or D", "No recorded stage")
    ))
poly_metan <- poly_metan_all |> 
  filter(stage %in% c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant", "No recorded stage")) |>
  mutate(stage = factor(
    stage, 
    levels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant", "No recorded stage"),
    labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or D", "No recorded stage")
  ))

# x-axis for proportions
prop_xscale <- function() {
  scale_x_continuous(
    limits = c(0, 0.75), 
    breaks = c(0, 0.25, 0.5), 
    labels = c("0%", "25%", "50%"), 
    minor_breaks = c(0.125, 0.375, 0.625), 
    expand = c(0,0), 
    name = "Receiving treatment (%)"
  )
}

p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(yvar = ordering, data = meta_dat |> mutate(include = ctry_class))
p <- p + incl_point(yvar = ordering, xvar = prop, data = meta_dat |> mutate(include = ctry_class))
# Meta-analysis
p <- p + meta_error(y_shift = 20) 
p <- p + meta_diamond(y_shift = 20, groupvar = stage)
# y axis
p <- p + jdn_yscale()
# x axis
p <- p + prop_xscale()
# legend
p <- p + guides(
  alpha = "none",
  colour = "none"
)
# faceting
p <- p + facet_grid(trt ~ stage)
# themes
p <- p + theme_bw() +
  theme_icbp() +
  scale_color_manual(values = colours) +
  scale_fill_manual(values = colours)
p
p |> saveRDS("results/figure1_lung_meta.RDS")

# Appendix version of plot ------------------------------------------------
meta_dat   <- meta_dat_all   |> 
  filter(stage %in% c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "Localised", "Regional", "Distant", "No recorded stage")) |>
  mutate(
    include = case_when(
      stage %in% c("Localised", "Regional", "Distant") ~ factor(0),
      TRUE ~ ctry_class
    )
  ) |>
  mutate(
    include = factor(
      include,
      levels = c("0", "United Kingdom", "Norway", "Canada", "Australia")
    )) |>
  mutate(
    stage = case_when(
      stage == "Localised" ~ "Stage 1",
      stage == "Regional" ~ "Stage 3",
      stage == "Distant" ~ "Stage 4",
      stage == "No recorded stage" ~ "Not recorded",
      TRUE ~ stage
    )
  ) |>
  mutate(stage = factor(
    stage, 
    levels = c(
      "Stage 1", 
      "Stage 2", 
      "Stage 3", 
      "Stage 4",
      "Not recorded"
    )
  )
  )

poly_metan <- poly_metan_all |> 
  filter(stage %in% c("Stage 1", "Stage 2", "Stage 3", "Stage 4", "No recorded stage")) |>
  mutate(
    stage = case_when(
      stage == "No recorded stage" ~ "Not recorded",
      TRUE ~ stage
    )
  ) |>
  mutate(stage = factor(
    stage, 
    levels = c(
      "Stage 1", 
      "Stage 2", 
      "Stage 3", 
      "Stage 4",
      "Not recorded"
    )
  )
  )

# x-axis for proportions
prop_xscale <- function() {
  scale_x_continuous(
    limits = c(0, 0.75), 
    breaks = c(0, 0.25, 0.5), 
    labels = c("0%", "25%", "50%"), 
    minor_breaks = c(0.125, 0.375, 0.625), 
    expand = c(0,0), 
    name = "Receiving treatment (%)"
  )
}

p <- ggplot()
# Included jurisdictions and estimates
p <- p + incl_error(
    yvar = ordering, 
    data = meta_dat
  )
p <- p + incl_point(yvar = ordering, xvar = prop)
# Meta-analysis
p <- p + meta_error(y_shift = 20) 
p <- p + meta_diamond(y_shift = 20, groupvar = stage)
# y axis
p <- p + jdn_yscale()
# x axis
p <- p + prop_xscale()
# legend
p <- p + guides(
  alpha = "none",
  colour = "none"
)
# faceting
p <- p + facet_grid(trt ~ stage)
# themes
p <- p + theme_bw() +
  theme_icbp() +
  scale_color_manual(values = colours) + 
  scale_fill_manual(values = c("white", colours))

p |> saveRDS("results/figure1_appendix_lung_meta.RDS")

meta_dat |> select(include) |> unique() |> mutate(i_n = as.numeric(include)) |> print()

# Data table for appendix -------------------------------------------------
appdx_table <- meta_dat_all |>
  group_by() |>
  filter(
    stage %in% c(
      "Stage 1", 
      "Stage 2", 
      "Stage 3", 
      "Stage 4", 
      "Localised",
      "Regional",
      "Distant",
      "No recorded stage"
    )
  ) |>
  mutate(
    stage = factor(
      stage,
      levels =  c(
        "Stage 1", 
        "Stage 2", 
        "Stage 3", 
        "Stage 4", 
        "Localised",
        "Regional",
        "Distant",
        "No recorded stage"
      )
    )
  ) |>
  select(trt, ctry_class, ctry_order, stage, n, n_trt, prop, lower, upper) |>
  arrange(trt, ctry_class, ctry_order, stage)

appdx_table |> 
  group_by( trt) |>
  gt() |>
  cols_align(align = "left", columns = c( trt, ctry_class, ctry_order, stage)) |>
  fmt_percent(
    columns = c(prop, lower, upper),
    decimals = 1
  ) |>
  fmt_number(
    columns = c(n, n_trt),
    decimals = 0
  ) |>
  tab_spanner(
    columns = c(n_trt, prop, lower, upper),
    label = "Received treatment"
  ) |>
  cols_merge(
    columns = c(lower, upper),
    pattern = "({1}, {2})"
  ) |>
  cols_label (
    trt = "Treatment",
    stage = "Stage at diagnosis",
    ctry_class = "Country",
    ctry_order = "Jurisdiction",
    prop = "%",
    n = "Patients",
    n_trt = "N",
    lower = "(95% CI)"
  ) |>
  saveRDS("results/z_appdxtable_figure1_lung_meta.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())

