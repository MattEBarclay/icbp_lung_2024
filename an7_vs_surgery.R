## Create plot of non-surgical treat vs surgery
library(tidyverse)

# Load in functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

dat <- readRDS("Rdata/dat_ind.RDS") |>
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |>
  filter(variable == "all") |>
  filter(stage == "All stages") |>
  mutate(trt = case_when(
    trt == "Chemo" ~ "Surgery vs chemotherapy",
    trt == "Radio" ~ "Surgery vs radiotherapy"
  )) |>
  filter(ctry_order %in% c("England", "New South Wales", "Ontario", "Norway", "Victoria", "Scotland")) |>
  mutate(surgery = case_when(
    ctry_order == "England"         ~ 0.149507,
    ctry_order == "New South Wales" ~ 0.193,
    ctry_order == "Ontario"         ~ 0.194896743,
    ctry_order == "Norway"          ~ 0.200755602,
    ctry_order == "Victoria"        ~ 0.259259259,
    ctry_order == "Scotland"        ~ 0.123311463
  ))

load("Rdata/useful_vectors.Rdata")

# England surgery = 2013-2016, NSCLC or SCLC, all cases

# Calculate Pearson correlation -------------------------------------------
calc_cor <- function(df, x, y) {
  
  x1 <- df |> select({{x}}) |> pull()
  y1 <- df |> select({{y}}) |> pull()
  
  cor.test(x1, y1) |> 
    broom::tidy() |>
    select(estimate, conf.low, conf.high, p.value)
}

c1 <- dat |> 
  filter(trt == "Surgery vs chemotherapy") |>
  calc_cor(x = prop, y = surgery)

c1 <- cbind(tibble(compare = "Surgery vs chemotherapy", stage = "All stages", year = " "), c1)

c2 <- dat |> 
  filter(trt == "Surgery vs radiotherapy") |>
  calc_cor(x = prop, y = surgery)

c2 <- cbind(tibble(compare = "Surgery vs radiotherapy", stage = "All stages", year = " "), c2)

rbind(c1, c2) |> saveRDS("results/corr_surg.RDS")
c1 |> rm()
c2 |> rm()

# Draw plots --------------------------------------------------------------
p <- ggplot(
  data = dat, 
  aes(
      x = prop, 
      y = surgery, 
      label = ctry_order, 
      colour = ctry_class,
      hjust = if_else(ctry_order %in% c("New South Wales"),  0.04,  0.9),
      vjust = if_else(ctry_order %in% c("New South Wales"),  2.5 ,  0)
    )
  )
p <- p + geom_point()
p <- p + geom_text(nudge_y = 0.005, size = 2.5)
p <- p + scale_x_continuous(
  name = "Receiving non-surgical treatment (%)",
  limits = c(0.2, 0.6),
  breaks = c(0.2, 0.3, 0.4, 0.5, 0.6),
  labels = c("20%", "30%", "40%", "50%", "60%"),
  minor_breaks = NULL
)
p <- p + scale_y_continuous(
  name = "Receiving surgery (%)",
  limits = c(0.1, 0.3),
  breaks = c(0.1, 0.15, 0.2, 0.25, 0.3),
  labels = c("10%", "15%", "20%", "25%", "30%"),
  minor_breaks = NULL
)
p <- p + facet_wrap(~trt)
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

p |> saveRDS("results/figure5_vs_surgery.RDS")


# Clean up ----------------------------------------------------------------
rm(list = ls())

