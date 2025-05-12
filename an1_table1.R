## Produce a table 1

library(tidyverse)
library(gt)

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

# Load in 'processed' data -------------------------------------------------------
dat_ind <- readRDS("Rdata/dat_ind.RDS")
load("Rdata/useful_vectors.Rdata")


# A table 1 ---------------------------------------------------------------

# Chemo patients | Radio patients
# Overall
# Jurisdiction (individual)
# Sex
# Age
# Stage
# Diagnosis year

overall <- dat_ind |> 
  mutate(variable_value = ifelse(variable == "sex", str_to_upper(variable_value), variable_value)) |>
  select(trt, stage, variable, variable_value, n, n_trt) |> 
  group_by(trt, stage, variable, variable_value) |> 
  summarise_all(sum)

jdns <- dat_ind |> 
  filter(variable == "all") |> 
  mutate(variable = "Jurisdiction") |>
  select(trt, stage, variable, ctry_order, n, n_trt) |> 
  group_by(trt, stage, variable, ctry_order) |> 
  summarise_all(sum) |>
  mutate(variable_value = as.character(ctry_order)) |>
  select(trt, stage, variable, variable_value, n, n_trt)

lung_full_table1 <- rbind(overall, jdns)
rm(overall, jdns)
lung_full_table1 |> ungroup() |> select(variable) |> unique() |> print(n = 100)
lung_full_table1 |> ungroup() |> select(variable_value) |> unique() |> print(n = 100)

lung_full_table1 <- lung_full_table1 |> 
  mutate(
    variable = factor(
      variable,
      levels = c(
        "all",
        "Jurisdiction",
        "sex",
        "age",
        "stage",
        "diagnosis_year"
      ),
      labels = c(
        "All",
        "Jurisdiction",
        "Sex",
        "Age group",
        "Stage",
        "Diagnosis year"
      )
    )
  )

lung_full_table1 <- lung_full_table1 |> 
  mutate(
    variable_value = factor(
      variable_value,
      levels = c(
        "all",           
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
        "New Brunswick",
        "Newfoundland & Labrador",        
        "Prince Edward Island",
        "Nova Scotia",         
        "New South Wales",     
        "Victoria",
        "F",
        "M",
        "64",
        "74",
        "84",
        "99",
        "1",
        "2",
        "3",
        "4",
        "L",
        "R",
        "D",
        "X",
        "2012",                
        "2013",              
        "2014",                
        "2015",               
        "2016",              
        "2017"
      ),
      labels = c(
        " ",
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
        "New Brunswick",
        "Newfoundland & Labrador",
        "Prince Edward Island",
        "Nova Scotia",         
        "New South Wales",     
        "Victoria",
        "Women",
        "Men",
        "15-64",
        "65-74",
        "75-84",
        "85-99",
        "TNM 1",
        "TNM 2",
        "TNM 3",
        "TNM 4",
        "Localised",
        "Regional",
        "Distant",
        "Not recorded",
        "2012",                
        "2013",              
        "2014",                
        "2015",               
        "2016",              
        "2017"
      )
    )
  )

lung_full_table1 <- lung_full_table1 |> arrange(variable, variable_value) |> ungroup()

# Look at stage cuts for all stages only
lung_full_table1 <- lung_full_table1 |>
  filter(stage == "All stages" | variable != "Stage")

lung_full_table1 <- lung_full_table1 |> 
  mutate(
    trt = case_when(
      trt == "Chemo" ~ "0",
      trt == "Radio" ~ "1"
    )
  ) |>
  pivot_wider(names_from = trt, values_from = c(n, n_trt)) |>
  select(stage, variable, variable_value, n_0, n_trt_0, n_1, n_trt_1)

lung_full_table1 <- lung_full_table1 |> 
  ungroup() |> 
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |>
  mutate(stage = factor(stage, levels = c("All stages", "Stages 1-3", "Stage 4"))) |>
  group_by(stage) |>
  arrange(stage, variable, variable_value) 
  
lung_full_table1 <- lung_full_table1 |> 
  mutate(prop_0 = n_trt_0/n_0) |>
  mutate(prop_1 = n_trt_1/n_1)

table1_processing <- function(data) {
  data |> 
    gt(rowname_col = "variable") |>
    cols_label(
      variable_value = "",
      n_0 = "Patients",
      n_trt_0 = "Treated",
      prop_0 = "(%)",
      n_1 = "Patients",
      n_trt_1 = "Treated",
      prop_1 = "(%)"
    ) |>
    tab_spanner(
      label = "Chemotherapy sample",
      columns = c(n_0, n_trt_0, prop_0)
    ) |>
    tab_spanner(
      label = "Radiotherapy sample",
      columns = c(n_1, n_trt_1, prop_1)
    ) |>
    fmt_number(
      columns = c(n_0, n_trt_0, n_1, n_trt_1),
      decimals = 0
    )  |>
    fmt_percent(
      columns = c(prop_0, prop_1),
      decimals = 1,
      pattern = "({x})"
    ) 
}

table1_all <- lung_full_table1 |> 
  table1_processing()

table1 <- lung_full_table1 |> filter(stage == "All stages") |> ungroup() |> select(-stage) |>
  table1_processing()


table1_all |> saveRDS("results/table1_all.RDS")
table1     |> saveRDS("results/table1.RDS")


# Clean up ----------------------------------------------------------------
rm(list = ls())

