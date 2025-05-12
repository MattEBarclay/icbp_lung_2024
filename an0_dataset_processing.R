## Turn raw datasets into analysis files
library(tidyverse)

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

# Load in 'raw' data -------------------------------------------------------
trt_unadj_n <- readRDS("Rdata/trt_unadj_n.RDS")
trt_q_perc  <- readRDS("Rdata/trt_q_perc.RDS")
load("Rdata/useful_vectors.Rdata")

trt_unadj_n |> 
  filter(trt == "Chemo") |>
  filter(jurisdiction == "Prince Edward Island") |>
  print(n = 100)


# Exclusions for lung cancer ----------------------------------------------
site_select <- function(data) {
  
  data |>
    filter(cancer == "Lung")
  
}

site_select_time <- function(data) {
  
  data |>
    filter(variable == "tumour_topography_group") |>
    filter(variable_value == "Lung") 
  
}

site_exclusions <- function(data, jurisdiction = jurisdiction, trt = trt) {
  
  data |>
    filter(
      # Data completeness problem, remove
      !(
        ({{jurisdiction}} == "Newfoundland & Labrador")  & 
        ({{trt}} == "Chemo")
      )
    )
      
}

# Create dataset for analysis --------------------------------------------------
dat <- trt_unadj_n |>
  site_select() |>
  site_exclusions() |>
  mutate(ctry_order = factor(
    jurisdiction, 
    levels = c(
      "England",
      "Northern Ireland",
      "Scotland",
      "Wales",
      "Norway",
      "Alberta",
      "British Columbia",
      "Ontario",
      "Saskatchewan",
      "Manitoba",
      "Prince Edward Island",
      "New Brunswick",
      "Newfoundland & Labrador",
      "Nova Scotia",
      "New South Wales",
      "Victoria"
    )
  )) |>
  select(trt, ctry_order, cancer, stage, variable, variable_value, n, n_trt, prop, lower, upper) |> 
  exclusion()


dat |> 
  filter(trt == "Chemo") |>
  filter(ctry_order == "Prince Edward Island") |>
  print(n = 100)

# Rationalise stage
dat <- dat |> rationalise_stage()

# add jurisdiction class
dat <- dat |> assign_ctry_class()

# group canada
dat <- dat |> group_canada()

# Summarise by grouped jurisdictions
dat <- dat |>
  select(  trt, ctry_class, ctry_order, cancer, stage, variable, variable_value, n, n_trt) |>
  group_by(trt, ctry_class, ctry_order, cancer, stage, variable, variable_value) |>
  summarise_all(sum) |>
  mutate(prop  = n_trt / n , 
         lower = (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2))),
         upper = (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2)))
  ) |>
  ungroup() 

# Add variables for plotting
dat <- dat |> plotting_variables()

dat |> saveRDS(file = "Rdata/dat.RDS")

# Create dataset for analysis - ungrouped jurisdictions ------------------------
dat_ind <- trt_unadj_n |>
  site_select() |>
  site_exclusions() |>
  mutate(ctry_order = factor(
    jurisdiction, 
    levels = c(
      "England",
      "Northern Ireland",
      "Scotland",
      "Wales",
      "Norway",
      "Alberta",
      "British Columbia",
      "Ontario",
      "Saskatchewan",
      "Manitoba",
      "Prince Edward Island",
      "New Brunswick",
      "Newfoundland & Labrador",
      "Nova Scotia",
      "New South Wales",
      "Victoria"
    )
  )) |>
  select(trt, ctry_order, cancer, stage, variable, variable_value, n, n_trt, prop, lower, upper) |> 
  exclusion()

# Rationalise stage
dat_ind <- dat_ind |> rationalise_stage()

# add jurisdiction class
dat_ind <- dat_ind |> assign_ctry_class()

# Summarise by grouped jurisdictions
dat_ind <- dat_ind |>
  select(  trt, ctry_class, ctry_order, cancer, stage, variable, variable_value, n, n_trt) |>
  group_by(trt, ctry_class, ctry_order, cancer, stage, variable, variable_value) |>
  summarise_all(sum) |>
  mutate(prop  = n_trt / n , 
         lower = (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) - (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2))),
         upper = (1/(1+(qnorm(0.975)^2)/n))*((n_trt/n)+(qnorm(0.975)^2)/(2*n)) + (qnorm(0.975)/(1+((qnorm(0.975)^2)/n)))*sqrt((n_trt/n)*(1-n_trt/n)/n + (qnorm(0.975)^2)/(4*(n^2)))
  ) |>
  ungroup() 

# Add variables for plotting
dat_ind <- dat_ind |> plotting_variables()

dat_ind |> saveRDS(file = "Rdata/dat_ind.RDS")

# Time-to-treat dataset ---------------------------------------------------
dat_time <- trt_q_perc |>
  site_select_time() |>
  site_exclusions() |>
  mutate(ctry_order = factor(
    jurisdiction, 
    levels = c(
      "England",
      "Northern Ireland",
      "Scotland",
      "Wales",
      "Norway",
      "Alberta",
      "British Columbia",
      "Ontario",
      "Saskatchewan",
      "Manitoba",
      "Prince Edward Island",
      "New Brunswick",
      "Newfoundland & Labrador",
      "Nova Scotia",
      "New South Wales",
      "Victoria"
    )
  )) |>
  select(trt, ctry_order, stage, variable, variable_value, percentile, days) |> 
  exclusion() |>
  assign_ctry_class() |>
  filter(
    percentile == 5  |
      percentile == 25 |
      percentile == 50 |
      percentile == 75 |
      percentile == 95
  ) |>
  select(trt, ctry_class, ctry_order, stage, percentile, days)

dat_time <- dat_time |> 
  pivot_wider(
    id_cols = c(trt, ctry_class, ctry_order, stage), 
    values_from = days, 
    names_from = percentile, 
    names_prefix = "pct_"
  )

dat_time <- dat_time |>
  mutate(stage = factor(stage, levels = c("All stages", "Stages 1-2", "Stages 1-3", "Stages 3-4", "Stage 4")))

dat_time |> saveRDS(file = "Rdata/dat_time.RDS")


# Cumulative treatment dataset --------------------------------------------

### Percent treated
dat_ind2 <- trt_unadj_n |>
  site_select() |>
  site_exclusions() |>
  mutate(ctry_order = factor(
    jurisdiction, 
    levels = c(
      "England",
      "Northern Ireland",
      "Scotland",
      "Wales",
      "Norway",
      "Alberta",
      "British Columbia",
      "Ontario",
      "Saskatchewan",
      "Manitoba",
      "Prince Edward Island",
      "New Brunswick",
      "Newfoundland & Labrador",
      "Nova Scotia",
      "New South Wales",
      "Victoria"
    )
  )) |>
  select(trt, ctry_order, cancer, stage, variable, variable_value, n, n_trt, prop, lower, upper) |> 
  exclusion()

# Rationalise stage
dat_ind2 <- dat_ind2 |> rationalise_stage()

# add jurisdiction class
dat_ind2 <- dat_ind2 |> assign_ctry_class()

#### Time to treatment
dat_time2 <- trt_q_perc |>
  site_select_time() |>
  site_exclusions() |>
  mutate(ctry_order = factor(
    jurisdiction, 
    levels = c(
      "England",
      "Northern Ireland",
      "Scotland",
      "Wales",
      "Norway",
      "Alberta",
      "British Columbia",
      "Ontario",
      "Saskatchewan",
      "Manitoba",
      "Prince Edward Island",
      "New Brunswick",
      "Newfoundland & Labrador",
      "Nova Scotia",
      "New South Wales",
      "Victoria"
    )
  )) |>
  select(trt, ctry_order, stage, variable, variable_value, percentile, days) |> 
  exclusion() |>
  assign_ctry_class() |>
  rename(cancer = variable_value) |>
  select(cancer, trt, ctry_class, ctry_order, stage, percentile, days)

# Pull separate bits together
a1 <- dat_ind2 |> 
  filter(variable == "all") |> 
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |> 
  group_by() |>
  select(cancer, trt, stage, ctry_class, ctry_order, prop, n, n_trt)
a2 <- dat_time2 |> 
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |> 
  group_by() |> 
  select(cancer, trt, stage, ctry_order, percentile, days)

x <- a1 |> inner_join(a2) |>
  mutate(n_trt_bit = case_when(
    percentile == 0 ~ 0,
    TRUE ~ n_trt*0.05
  )) |> 
  group_by(ctry_order, cancer, trt, stage) |> 
  arrange(ctry_order, cancer, trt, stage, days) |>
  mutate(n_trt_cum = cumsum(n_trt_bit)) |>
  rename(n_trt_total = n_trt)

dat_cumulative <- x |> 
  mutate(n_trt = n_trt_cum) |>
  mutate(prop = n_trt/n) |>
  mutate(percentile = round(100*n_trt/n_trt_total,0)) |>
  select(cancer, trt, stage, ctry_order, ctry_class, percentile, days, n, n_trt, prop) |> 
  arrange(cancer, trt, stage, ctry_order, ctry_class, percentile, days) 

rm(a1, a2, dat_ind2, dat_time2, x)

dat_cumulative |> saveRDS(file = "Rdata/dat_cumulative.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())

