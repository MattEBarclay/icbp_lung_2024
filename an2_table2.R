## Produce a table 2
## Stage distributions by jurisdiction

library(tidyverse)
library(gt)

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

# Load in 'processed' data -------------------------------------------------------
dat_ind <- readRDS("Rdata/dat_ind.RDS")
load("Rdata/useful_vectors.Rdata")


# Table 2 - jurisdictional stage distributions ---------------------------------

# 2A - TNM jurisdictions
lung_stage_dist_table1 <- dat_ind |> 
  filter(ctry_order != "Norway" & ctry_order != "New South Wales" & ctry_order != "Victoria") |>
  filter(stage == "All stages") |> 
  filter(variable == "stage") |>
  filter(trt == "Radio" | ctry_order == "Northern Ireland") |>
  select(trt, ctry_order, variable_value, n) |> 
  pivot_wider(id_cols = c(trt, ctry_order), values_from = n, names_from = variable_value, values_fill = 0, names_prefix = "tnm") |>
  arrange(ctry_order) |>
  mutate(tnm1_pct = tnm1/(tnm1+tnm2+tnm3+tnm4+tnmX)) |>
  mutate(tnm2_pct = tnm2/(tnm1+tnm2+tnm3+tnm4+tnmX)) |>
  mutate(tnm3_pct = tnm3/(tnm1+tnm2+tnm3+tnm4+tnmX)) |>
  mutate(tnm4_pct = tnm4/(tnm1+tnm2+tnm3+tnm4+tnmX)) |>
  mutate(tnmX_pct = tnmX/(tnm1+tnm2+tnm3+tnm4+tnmX)) |>
  mutate(tumours = tnm1+tnm2+tnm3+tnm4+tnmX) |>
  mutate(ctry_name = as.character(ctry_order)) |>
  mutate(ctry_name = case_when(
    ctry_order == "Northern Ireland" ~ paste0(ctry_name, " (", trt, ")"), 
    TRUE ~ ctry_name)
  ) |>
  arrange(ctry_order, trt) |>
  select(ctry_name, tumours, tnm1, tnm1_pct, tnm2, tnm2_pct, tnm3, tnm3_pct, tnm4, tnm4_pct, tnmX, tnmX_pct) |>
  gt() |>
  cols_label (
    ctry_name = "Jurisdiction",
    tumours = "Lung cancers",
    tnm1 = "N",
    tnm2 = "N",
    tnm3 = "N",
    tnm4 = "N",
    tnmX = "N",
    tnm1_pct = "(%)",
    tnm2_pct = "(%)",
    tnm3_pct = "(%)",
    tnm4_pct = "(%)",
    tnmX_pct = "(%)",
  ) |>
  tab_spanner(
    label = "TNM 1",
    columns = c(tnm1, tnm1_pct)
  ) |>
  tab_spanner(
    label = "TNM 2",
    columns = c(tnm2, tnm2_pct)
  ) |>
  tab_spanner(
    label = "TNM 3",
    columns = c(tnm3, tnm3_pct)
  ) |>
  tab_spanner(
    label = "TNM 4",
    columns = c(tnm4, tnm4_pct)
  ) |>
  tab_spanner(
    label = "Unknown",
    columns = c(tnmX, tnmX_pct)
  ) |>
  fmt_percent (
    columns = c(tnm1_pct, tnm2_pct, tnm3_pct, tnm4_pct, tnmX_pct),
    decimals = 1,
    pattern = "({x})"
  ) |>
  fmt_number(
    columns = c(tumours, tnm1, tnm2, tnm3, tnm4, tnmX),
    decimals = 0
  ) |>
  cols_align(align = "center") |>
  cols_align(align = "left", columns = c(ctry_name)) |>
  tab_style(
    style = cell_text(
      font = "Arial",
      size = 8
    ),
    locations = list(
      cells_body(), 
      cells_column_labels()
    )
  )


# 2A - SEER staging jurisdictions
lung_stage_dist_table2 <- dat_ind |> 
  filter(ctry_order == "Norway" | ctry_order == "New South Wales" | ctry_order == "Victoria") |>
  filter(stage == "All stages") |> 
  filter(variable == "stage") |>
  filter(trt == "Radio" | ctry_order == "Northern Ireland") |>
  select(trt, ctry_order, variable_value, n) |> 
  pivot_wider(id_cols = c(trt, ctry_order), values_from = n, names_from = variable_value, values_fill = 0) |>
  arrange(ctry_order) |>
  mutate(L_pct = L/(L+R+D+X)) |>
  mutate(R_pct = R/(L+R+D+X)) |>
  mutate(D_pct = D/(L+R+D+X)) |>
  mutate(X_pct = X/(L+R+D+X)) |>
  mutate(tumours = L+R+D+X) |>
  mutate(ctry_name = as.character(ctry_order)) |>
  mutate(ctry_name = case_when(
    ctry_order == "Northern Ireland" ~ paste0(ctry_name, " (", trt, ")"), 
    TRUE ~ ctry_name)
  ) |>
  arrange(ctry_order, trt) |>
  select(ctry_name, tumours, L, L_pct, R, R_pct, D, D_pct, X, X_pct) |>
  gt() |>
  cols_label (
    ctry_name = "Jurisdiction",
    tumours = "Lung cancers",
    L = "N",
    R = "N",
    D = "N",
    X = "N",
    L_pct = "(%)",
    R_pct = "(%)",
    D_pct = "(%)",
    X_pct = "(%)",
  ) |>
  tab_spanner(
    label = "Localised",
    columns = c(L, L_pct)
  ) |>
  tab_spanner(
    label = "Regional",
    columns = c(R, R_pct)
  ) |>
  tab_spanner(
    label = "Distant",
    columns = c(D, D_pct)
  ) |>
  tab_spanner(
    label = "Unknown",
    columns = c(X, X_pct)
  ) |>
  fmt_percent (
    columns = c(L_pct, R_pct, D_pct, X_pct),
    decimals = 1,
    pattern = "({x})"
  ) |>
  fmt_number(
    columns = c(tumours, L, R, D, X),
    decimals = 0
  ) |>
  cols_align(align = "center") |>
  cols_align(align = "left", columns = c(ctry_name)) |>
  tab_style(
    style = cell_text(
      font = "Arial",
      size = 8
    ),
    locations = list(
      cells_body(), 
      cells_column_labels()
    )
  )


lung_stage_dist_table1 |> saveRDS("results/table2A_TNM.RDS")
lung_stage_dist_table2 |> saveRDS("results/table2B_SEER.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())


