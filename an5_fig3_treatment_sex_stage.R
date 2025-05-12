## Create plot of treatment by age, stage
library(tidyverse)

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

# Load in 'raw' data -------------------------------------------------------

# Jdn y-pos
plot_ordering <- function(data) {
  data |>
    mutate(ordering    = as.numeric(ctry_order)) |>
    mutate(ordering    = ifelse(ordering >=  6, ordering+1, ordering)) |>
    mutate(ordering    = ifelse(ordering >= 12, ordering+1, ordering)) 
}

# Jdn labels for plot
country_order_grp <- c(
    "England",              
    "Northern Ireland",     
    "Scotland",             
    "Wales",                
    "Norway",               
    "Alberta",              
    "British Columbia",     
    "Ontario",              
    "SK-MB",
    "Atlantic Canada",      
    "New South Wales",      
    "Victoria"
  )

# data processing
dat <- readRDS("Rdata/dat.RDS") |>
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |>
  mutate(stage = factor(stage, levels = c("All stages", "Stages 1-3", "Stage 4"))) |>
  mutate(trt = case_when(
    trt == "Chemo" ~ "Chemotherapy",
    trt == "Radio" ~ "Radiotherapy"
  )) |> 
  plot_ordering() |> 
  filter(variable == "sex") |>
  mutate(offset = case_when(
    sex == "F" ~  0.1,
    sex == "M" ~ -0.1,
    TRUE ~ 0
  )) |>
  mutate(stage = factor(stage, labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant")))
  
# Draw plot --------------------------------------------
p <- ggplot(data = dat)
# error bars
p <- p + geom_errorbarh(
  linetype = 1,
  size = 0.3,
  height = 0,
  aes(
    xmin = lower,
    xmax = upper,
    y    = ordering+offset,
    colour = sex
  )
)
# observed value
p <- p + geom_point(
  aes(
    x = prop,
    y = ordering+offset,
    colour = sex,
    shape = sex
  )
)
# X scale - proportions
p <- p + scale_x_continuous(
  name = "Receiving treatment (%)",
  limits = c(0, 0.625),
  breaks = c(0, 0.25, 0.5),
  minor_breaks = c(0.125, 0.375, 0.625),
  labels = c("0%", "25%", "50%"),
  expand = c(0,0)
)
# Y scale - countries
p <- p + scale_y_continuous(
  breaks = c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 14),
  minor_breaks = NULL,
  labels = country_order_grp,
  name = "",
  transform = "reverse"
  )
# Theme
p <- p +
  theme_bw() +
  theme_dotplots() +
  guides(
    alpha = "none",
    size = "none",
    fill = "none"
  ) +
  scale_color_manual(
    name = "",
    breaks = c("F", "M"),
    labels = c("Women", "Men"),
    values = c("#BC1C80", "#18ACBD")
  ) +
  scale_shape_manual(
    name = "",
    breaks = c("F", "M"),
    labels = c("Women", "Men"),
    values = c(15,16)
  ) +
  theme(legend.position = "top", legend.justification = c(1,0)) 
# One column for each stage
p <- p + facet_grid(trt ~ stage)
p 

p |> saveRDS("results/figure3_sex_stage.RDS")

# Data table for appendix -------------------------------------------------
appdx_table <- dat |>
  group_by() |>
  select(trt, ctry_class, ctry_order, stage, sex, n, n_trt, prop, lower, upper) |>
  filter(!is.na(sex)) |>
  arrange(trt, ctry_class, ctry_order, stage, sex) 

appdx_table |> 
  group_by( trt, stage) |>
  gt() |>
  cols_align(align = "left", columns = c( trt, ctry_class, ctry_order, stage, sex)) |>
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
    sex = "Sex",
    ctry_class = "Country",
    ctry_order = "Jurisdiction",
    prop = "%",
    n = "Patients",
    n_trt = "N",
    lower = "(95% CI)"
  ) |>
  saveRDS("results/z_appdxtable_figure3_sex_stage.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())

