## Create plot of chemo vs radio, incl by stage
library(tidyverse)

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

dat <- readRDS("Rdata/dat_ind.RDS") |>
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |>
  filter(variable == "all") |>
  select(
    trt,
    ctry_class,
    ctry_order,
    stage,
    n,
    n_trt,
    prop,
    lower,
    upper
  ) |>
  mutate(
    stage = factor(
      stage, 
      levels = c("All stages", "Stages 1-3", "Stage 4"),
      labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant"),
      ))

dat <- dat |>
  pivot_wider(
    id_cols = c(ctry_class, ctry_order, stage), 
    names_from = trt, 
    values_from = c(n, n_trt, prop, lower, upper)
  )

load("Rdata/useful_vectors.Rdata")

# Calculate Pearson correlation -------------------------------------------
calc_cor <- function(df, x, y) {
  
  x1 <- df |> select({{x}}) |> pull()
  y1 <- df |> select({{y}}) |> pull()
  
  cor.test(x1, y1) |> 
    broom::tidy() |>
    select(estimate, conf.low, conf.high, p.value)
}

dat |> select(stage) |> unique()

# All
c1 <- dat |> 
  filter(stage == "All stages") |>
  calc_cor(x = prop_Chemo, y = prop_Radio)

c1 <- cbind(tibble(compare = "Chemotherapy vs radiotherapy", stage = "All stages", year = " "), c1)

# Non-advanced
c2 <- dat |> 
  filter(stage == "Stages 1-3 or L-R") |>
  calc_cor(x = prop_Chemo, y = prop_Radio)

c2 <- cbind(tibble(compare = "Chemotherapy vs radiotherapy", stage = "Stages 1-3 or L-R", year = " "), c2)

# Advanced
c3 <- dat |> 
  filter(stage == "Stage 4 or distant") |>
  calc_cor(x = prop_Chemo, y = prop_Radio)

c3 <- cbind(tibble(compare = "Chemotherapy vs radiotherapy", stage = "Stage 4 or distant", year = " "), c3)

rbind(c1, c2, c3) |> saveRDS("results/corr_trt.RDS")
c1 |> rm()
c2 |> rm()
c3 |> rm()

# Draw plots --------------------------------------------------------------
p <- ggplot(data = dat, aes(x = prop_Chemo, y = prop_Radio, label = ctry_order, colour = ctry_class))
p <- p + geom_point()
#p <- p + geom_text(nudge_y = 0.005, size = 3)
p <- p + scale_x_continuous(
  name = "Receiving chemotherapy (%)",
  limits = c(0.2, 0.6),
  breaks = c(0.2, 0.3, 0.4, 0.5, 0.6),
  labels = c("20%", "30%", "40%", "50%", "60%"),
  minor_breaks = NULL
)
p <- p + scale_y_continuous(
  name = "Receiving radiotherapy (%)",
  limits = c(0.2, 0.6),
  breaks = c(0.2, 0.3, 0.4, 0.5, 0.6),
  labels = c("20%", "30%", "40%", "50%", "60%"),
  minor_breaks = NULL
)
p <- p + facet_wrap(~stage)
p <- p +
  theme_bw() +
  theme_dotplots() + 
  theme(legend.position = "bottom") + 
  guides(
    alpha = "none",
    size = "none",
    fill = "none",
    shape = "none"
  ) +
  scale_colour_manual('Country', values = c("#984ea3", "#e41a1c", "#4daf4a", "#377eb8")) +
  scale_fill_manual(values = c("#984ea3", "#e41a1c", "#4daf4a", "#377eb8"))

p |> saveRDS("results/figure6_vs_chemo_radio.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())

