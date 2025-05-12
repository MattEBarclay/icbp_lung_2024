# Draw cumulative treatment plots...

# Idea here is to draw a separate plot for each country, with
# the jdns for other countries in grey

library(tidyverse)
library(gt)
#library(gganimate)
#library(transformr)

## colour scheme
#colours <- c("#BC1C80","#18ACBD", "#BC1C80", "#BC1C80")
#colours <- c("#54278f", "#756bb1", "#9e9ac8", "#cbc9e2",
#             "#fd8d3c", 
#             "#08519c", "#3182bd", "#6baed6", "#31a354", "#a1d99b",
#             "#BC1C80", "#f09ad0")
#colours <- c("#003566", "#005fb6", "#0084ff", "#49a7ff",
#             "#524A6F", 
#             "#0d5e67", "#128390", "#18a8b9", "#1dcee2", "#46d7e7",
#             "#BC1C80", "#f09ad0")
colours <- c("#08519c", "#08519c", "#9ecae1", "#9ecae1",
             "#df65b0", 
             "#006d2c", "#006d2c", "#006d2c", "#a1d99b", "#a1d99b",
             "#fb6a4a", "#fb6a4a")
linetypes <- c(1, 2, 1, 2,
             1, 
             1, 2, 4, 1, 2,
             1, 2)

# Helpful functions for processing time-to-treat data ---------------------
source("icbp_cancerspec_lung_functions_script.R")

# Load data, final processing ---------------------------------------------
dat_cumulative <- readRDS("Rdata/dat_cumulative.RDS") |>
  filter(stage %in% c("All stages", "Stages 1-3", "Stage 4")) |>
  mutate(stage = factor(
    stage, 
    levels = c("All stages", "Stages 1-3", "Stage 4"),
    labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant")
  )) |>
  mutate(trt = case_when(
    trt == "Chemo" ~ "Chemotherapy",
    trt == "Radio" ~ "Radiotherapy"
  )) |> ungroup()

dat_cumulative <- dat_cumulative |> 
  mutate(ctry_order2 = ctry_order) |>
  group_canada() |>
  select(-jdn_class, -chk) |>
  rename(groups = ctry_order) |>
  rename(ctry_order = ctry_order2)

# Load data, final processing ---------------------------------------------

# Draw the static cumulative plot -----------------------------------------

# Data to plot
plot_data <- dat_cumulative |> 
  calculate_grouped_centiles() |>
  improve_centiles() |>
  filter_cumulative(
    centile = 100,
    plot_these = dat_cumulative |> select(ctry_class) |> unique() |> arrange() |> pull(),
    to_label = dat_cumulative |> select(ctry_class) |> unique() |> arrange() |> pull(),
    nudge = NULL
  ) |>
  rename(ctry_order = groups)

plot_function <- function(ctry_plot = "Canada") {
  
  # Fiddle with the data
  plot_data_new <- plot_data |> 
    mutate(ctry_order2 = case_when(
      ctry_class == ctry_plot ~ ctry_order,
      TRUE ~ "Other jurisdictions"
    ))
  factor_labels <- plot_data_new |> select(ctry_order2) |> unique() |> filter(ctry_order2 != "Other jurisdictions") |> pull()
  plot_data_new <- plot_data_new |>
    mutate(ctry_order2 = factor(ctry_order2, levels = c( factor_labels, "Other jurisdictions"))) |>
    arrange(ctry_order2, days)
  
  # Colour scheme
  c1 <- tibble(ctry = "Canada"        , colours = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "lightgrey") )
  c2 <- tibble(ctry = "Norway"        , colours = c("#e41a1c", "lightgrey"))
  c3 <- tibble(ctry = "Australia"     , colours = c("#e41a1c", "#377eb8", "lightgrey") )
  c4 <- tibble(ctry = "United Kingdom", colours = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "lightgrey") )
  colour_schemes <- rbind(c1, c2, c3, c4)
  rm(c1, c2, c3, c4)
  colours <- colour_schemes |> filter(ctry == ctry_plot) |> select(colours) |> pull()

  p <- ggplot(
    data = plot_data_new,
    aes(
      x = days, 
      y = prop, 
      group = ctry_order, 
      col = ctry_order2,
      linewidth = 1,
      label = label
    ),
  )
  # draw in separate passes to try to sort out layering...
  p <- p + geom_line(data = plot_data_new |> filter(ctry_order2 == "Other jurisdictions"), col = "lightgrey")
  p <- p + geom_line(data = plot_data_new |> filter(ctry_order2 != "Other jurisdictions"), aes(col = ctry_order2))
  p <- p + scale_x_continuous(
    name = "Days from diagnosis",
    limits = c(-31, 365),
    breaks = c(0, 90, 180, 270, 365),
    minor_breaks = c(0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330)
  )
  p <- p + scale_y_continuous(
    limits = c(0, 0.625),
    breaks = c(0, 0.25, 0.5),
    labels = c("0%", "25%", "50%"),
    minor_breaks = c(0.125, 0.375, 0.625),
    expand = c(0,0), 
    name = "Receiving treatment (%)"
  )
  p <- p +
    theme_bw() +
    # standard theme elements
    theme(
      strip.text.x = element_text(size = 8, color = "black", hjust= 0 ),
      strip.background = element_rect(fill = "grey"),
      axis.text.x = element_text(size = 7, color = "black", vjust = 0.5, angle = 0),
      axis.text.y = element_text(size = 7, color = "black"),
      axis.ticks.length = unit(0, "cm"),
      panel.spacing.x = unit(0.5, "lines"),
      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
      panel.spacing.y = unit(0.3,"line")
    ) +
    guides(
      alpha = "none",
      size = "none",
      fill = "none",
      shape = "none",
      linewidth = "none"
    )  +
    scale_linewidth( range = c(0.5, 1.3) )
  p <- p + facet_grid(trt ~ stage)
  p <- p + scale_color_manual(values = colours)
  p <- p + theme(legend.title = element_blank())
  p <- p + geom_vline(xintercept = 0)
  p
}

plot_function(ctry = "Canada") |> saveRDS("results/figure4b_cumulative_jdn_canada.RDS")
plot_function(ctry = "Australia") |> saveRDS("results/figure4b_cumulative_jdn_australia.RDS")
plot_function(ctry = "Norway") |> saveRDS("results/figure4b_cumulative_jdn_norway.RDS")
plot_function(ctry = "United Kingdom") |> saveRDS("results/figure4b_cumulative_jdn_uk.RDS")

# Appendix table showing key percentages ----------------------------------
table <- dat_cumulative |> 
  calculate_grouped_centiles() |>
  improve_centiles() |>
  filter(percentile %in% c(25,50,75,90)) |>
  rename(ctry_order = groups)

table |> 
  arrange(trt, stage, ctry_order, percentile) |> 
  filter(ctry_order %in% c("Saskatchewan-Manitoba", "Atlantic Canada")) |> 
  print(n=100)

table <- table |> 
  select(cancer, trt, stage, ctry_class, ctry_order, percentile, days, prop) |>
  group_by(
  ) |>
  mutate(stage_n = as.numeric(stage)) |> select(-stage) |>
  pivot_wider(
    id_cols = c(cancer, trt, ctry_class, ctry_order, percentile),
    names_from = stage_n,
    values_from = c(days, prop)
    ) |>
  select(cancer, trt, ctry_class, ctry_order, percentile, days_1, prop_1, days_2, prop_2, days_3, prop_3)

table <- table |> 
  group_by(cancer, trt) |>
  gt() |>
  cols_align(align = "left", columns = c(cancer, trt, ctry_class, ctry_order)) |>
  fmt_percent(
    columns = c(prop_1, prop_2, prop_3),
    decimals = 1
  ) |>
  fmt_number(
    columns = c(percentile, days_1, days_2, days_3),
    decimals = 0
  ) |>
  tab_spanner(
    columns = c(days_1, prop_1),
    label = "All lung cancers"
  ) |>
  tab_spanner(
    columns = c(days_2, prop_2),
    label = "Stages 1-3 or L-R"
  ) |>
  tab_spanner(
    columns = c(days_3, prop_3),
    label = "Stage 4 or distant"
  ) |>
  cols_label (
    cancer = "Cancer",
    trt = "Treatment",
    ctry_class = "Country",
    ctry_order = "Jurisdiction",
    percentile = "Percentile of those treated",
    days_1 = "Days",
    prop_1 = "% treated",
    days_2 = "Days",
    prop_2 = "% treated",
    days_3 = "Days",
    prop_3 = "% treated"
  ) 

table |> saveRDS("results/z_appdxtable_figure4b_cumulative_jdn.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())
