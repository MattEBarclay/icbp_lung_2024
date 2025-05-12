## Create plot of treatment by age, stage
library(tidyverse)
library(gt)
library(patchwork)

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

# Load in 'raw' data -------------------------------------------------------
dat <- readRDS("Rdata/dat.RDS") |>
  mutate(ctry_order = factor(ctry_order, labels = c(
    "England",              
    "NI",     
    "Scotland",             
    "Wales",                
    "Norway",               
    "Alberta",              
    "BC",     
    "Ontario",              
    "SK-MB",
    "Atl. Canada",      
    "NSW",      
    "Victoria"    
  )))

# Overall age treatment trends --------------------------------------------
chemo0 <- trend_plot(
  df = dat,
  variable_touse = "age",
  trt_touse = "Chemo",
  panel_title = "Chemotherapy, all stages",
  stages = "All stages",
  x_var = age,
  x_name = "Age at diagnosis",
  x_break = c(1, 2, 3, 4),
  x_limit = c(0.5, 4.5),
  x_label = c("15-64", "65-74", "75-84", "85-99"),
  x_angle = 30,
  y_name = "Receiving treatment (%)",
  labfill = "grey"
)

chemo1 <- trend_plot(
  df = dat,
  variable_touse = "age",
  trt_touse = "Chemo",
  panel_title = "Chemotherapy, stages 1-3 or L-R",
  stages = "Stages 1-3",
  x_var = age,
  x_name = "Age at diagnosis",
  x_break = c(1, 2, 3, 4),
  x_limit = c(0.5, 4.5),
  x_label = c("15-64", "65-74", "75-84", "85-99"),
  x_angle = 30,
  y_name = "Receiving treatment (%)",
  labfill = "grey"
)

chemo2 <- trend_plot(
  df = dat,
  variable_touse = "age",
  trt_touse = "Chemo",
  panel_title = "Chemotherapy, stage 4 or distant",
  stages = "Stage 4",
  x_var = age,
  x_name = "Age at diagnosis",
  x_break = c(1, 2, 3, 4),
  x_limit = c(0.5, 4.5),
  x_label = c("15-64", "65-74", "75-84", "85-99"),
  x_angle = 30,
  y_name = "Receiving treatment (%)",
  labfill = "grey"
)

radio0 <- trend_plot(
  df = dat,
  variable_touse = "age",
  trt_touse = "Radio",
  panel_title = "Radiotherapy, all stages",
  stages = "All stages",
  x_var = age,
  x_name = "Age at diagnosis",
  x_break = c(1, 2, 3, 4),
  x_limit = c(0.5, 4.5),
  x_label = c("15-64", "65-74", "75-84", "85-99"),
  x_angle = 30,
  y_name = "Receiving treatment (%)",
  labfill = "grey"
)

radio1 <- trend_plot(
  df = dat,
  variable_touse = "age",
  trt_touse = "Radio",
  panel_title = "Radiotherapy, stages 1-3 or L-R",
  stages = "Stages 1-3",
  x_var = age,
  x_name = "Age at diagnosis",
  x_break = c(1, 2, 3, 4),
  x_limit = c(0.5, 4.5),
  x_label = c("15-64", "65-74", "75-84", "85-99"),
  x_angle = 30,
  y_name = "Receiving treatment (%)",
  labfill = "grey"
)

radio2 <- trend_plot(
  df = dat,
  variable_touse = "age",
  trt_touse = "Radio",
  panel_title = "Radiotherapy, stage 4 or distant",
  stages = "Stage 4",
  x_var = age,
  x_name = "Age at diagnosis",
  x_break = c(1, 2, 3, 4),
  x_limit = c(0.5, 4.5),
  x_label = c("15-64", "65-74", "75-84", "85-99"),
  x_angle = 30,
  y_name = "Receiving treatment (%)",
  labfill = "grey"
)

six_panel = chemo0 + chemo1 + chemo2 + radio0 + radio1 + radio2 + plot_layout(axis_titles = "collect")
six_panel

six_panel |> saveRDS("results/figure2_age_stage.RDS")

# Data table for appendix -------------------------------------------------
appdx_table <- dat |>
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |>
  group_by() |>
  select(trt, ctry_class, ctry_order, stage, age, n, n_trt, prop, lower, upper) |>
  filter(!is.na(age)) |>
  mutate(stage = factor(
    stage,
    levels = c("All stages", "Stages 1-3", "Stage 4"),
    labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant")
  )) |>
  arrange(trt, ctry_class, ctry_order, stage, age) 

appdx_table |> 
  group_by( trt, stage) |>
  gt() |>
  cols_align(align = "left", columns = c( trt, ctry_class, ctry_order, stage, age)) |>
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
    age = "Age at diagnosis",
    ctry_class = "Country",
    ctry_order = "Jurisdiction",
    prop = "%",
    n = "Patients",
    n_trt = "N",
    lower = "(95% CI)"
  ) |>
  saveRDS("results/z_appdxtable_figure2_age_stage.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())
