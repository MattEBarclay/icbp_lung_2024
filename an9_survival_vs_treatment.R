# Plots survival vs treatment rates, by stage if possible
library(tidyverse)

# Repeated data functions -------------------------------------------------------
source("icbp_cancerspec_lung_functions_script.R")

data_clean <- function(df) {
  # Make values nicer
  df <- df |> mutate(country = str_to_title(country))
  df <- df |> mutate(jurisdiction = str_to_title(jurisdiction))
  df <- df |> mutate(cancer = str_to_title(cancer))
  
  # Correct some formatting / typos
  df <- df |>
    mutate(
      jurisdiction = case_when(
        jurisdiction == "Newfoundland And Labrador" ~ "Newfoundland & Labrador",
        TRUE ~ jurisdiction
      ))
  
  df <- df |>
    mutate(
      country = case_when(
        country == "Uk" ~ "UK",
        country == "Austrlia" ~ "Australia",
        TRUE ~ country
      ))
  
  df <- df |>
    mutate(
      country = case_when(
        country == "Uk" ~ "UK",
        TRUE ~ country
      ))
  
}

# Reshape dataset
survival_reshape <- function(df) {
  df |>
    pivot_longer(
      cols = starts_with("netsurv"),
      names_to = c(".value", "est", "year"),
      names_sep = "_"
    ) |>
    mutate(netsurv = as.numeric(netsurv)) |>
    mutate(year = as.numeric(year)) |>
    pivot_wider(
      values_from = netsurv,
      names_from = est,
      names_prefix = "surv_"
    ) |>
    rename(surv = surv_pt)
}

factor_ctry <- function(df) {
# Make factors for countries
  df |>
    mutate(
      country = factor(
        country, 
        levels = c(
          "UK", 
          "Norway", 
          "Canada", 
          "Australia"
        )
      )
    )
}

factor_jdns <- function(df) {
  df |> 
    mutate(
      jurisdiction = factor(
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
      )
    )
}

# Load in survival data ------------------------------------------------------------

# Load all-lung survival data and fix data frame
surv_all <- readxl::read_excel("Rdata/survival_icbp_estimates-lung-ovarian.xlsx", sheet = "all_stage_1y_and_5y", .name_repair = "universal")
surv_all <- surv_all |> 
  rename(netsurv_pt_1 = ..1y_net_surv   ) |> 
  rename(netsurv_lb_1 = ..1y_net_surv_lb) |>
  rename(netsurv_ub_1 = ..1y_net_surv_ub) |>
  rename(netsurv_pt_5 = ..5y_net_surv   ) |>
  rename(netsurv_lb_5 = ..5y_net_surv_lb) |>
  rename(netsurv_ub_5 = ..5y_net_surv_ub) 

surv_all <- surv_all |> 
  data_clean() |> 
  factor_jdns() |>
  factor_ctry() |>
  survival_reshape()

# Load NSCLC stage-spec survival data and fix data frame
surv_stage <- readxl::read_excel("Rdata/survival_icbp_estimates-lung-ovarian.xlsx", sheet = "NSCLC_stage_specific_1y_and_3y", .name_repair = "universal")
surv_stage <- surv_stage |> 
  rename(netsurv_pt_1 = ..1y_net_surv   ) |> 
  rename(netsurv_lb_1 = ..1y_net_surv_lb) |>
  rename(netsurv_ub_1 = ..1y_net_surv_ub) |>
  rename(netsurv_pt_3 = ..3y_net_surv   ) |>
  rename(netsurv_lb_3 = ..3y_net_surv_lb) |>
  rename(netsurv_ub_3 = ..3y_net_surv_ub) 

surv_stage <- surv_stage |> 
  data_clean() |> 
  factor_jdns() |>
  factor_ctry() |>
  survival_reshape()

surv_all   <- surv_all   |> filter(cancer == "Lung") |> select(-cancer) 
surv_stage <- surv_stage |> filter(cancer == "Lung") |> select(-cancer, -source)

surv <- surv_all |> mutate(stage = "all") |> rbind(surv_stage)
surv <- surv |>
  mutate(stage = factor(stage, levels = c("all", "1", "2", "3", "4", "L", "R", "D"))) |> 
  select(stage, country, jurisdiction, year, surv, surv_lb, surv_ub) |>
  arrange(stage, country, jurisdiction, year)

surv_all |> rm()
surv_stage |> rm()

# Load in treatment data --------------------------------------------------
trt <- readRDS("Rdata/dat_ind.RDS") |>
  filter(variable %in% c("stage", "all")) |>
  filter(stage == "All stages") |>
  filter(cancer == "Lung") |>
  filter(variable_value != "X") |> ## don't care about missing stage
  select(-ctry_class, -variable, -cancer, -age, -year, -stage, -stage_x, -sex) |>
  mutate(stage = factor(variable_value, levels = c("all", "1", "2", "3", "4", "L", "R", "D"))) |> 
  rename(jurisdiction = ctry_order) |>
  rename(trt_prop = prop) |>
  rename(trt_prop_lb = lower) |>
  rename(trt_prop_ub = upper) |>
  select(trt, stage, jurisdiction, n, n_trt, trt_prop, trt_prop_lb, trt_prop_ub) |>
  arrange(stage, jurisdiction) |>
  mutate(trt = case_when(
    trt == "Chemo" ~ "Chemotherapy",
    trt == "Radio" ~ "Radiotherapy",
    TRUE ~ trt
  ))

surg <- trt |> 
  filter(trt == "Chemotherapy") |>
  filter(stage == "all") |> 
  filter(jurisdiction %in% c("England", "New South Wales", "Ontario", "Norway", "Victoria", "Scotland")) |>
  mutate(trt = "Surgery") |>
  mutate(trt_prop = case_when(
    jurisdiction == "England"         ~ 0.149507,
    jurisdiction == "New South Wales" ~ 0.193,
    jurisdiction == "Ontario"         ~ 0.194896743,
    jurisdiction == "Norway"          ~ 0.200755602,
    jurisdiction == "Victoria"        ~ 0.259259259,
    jurisdiction == "Scotland"        ~ 0.123311463
  )) |>
  mutate(n = NA) |>
  mutate(n_trt = NA) |>
  mutate(trt_prop_lb = NA) |>
  mutate(trt_prop_ub = NA)

trt <- trt |> rbind(surg)
surg |> rm()

# merge treatment and survival data ---------------------------------------
plot_data <- trt |> merge(surv)

# Calculate Pearson correlation -------------------------------------------
calc_cor <- function(df, x, y) {
  cor.test(x, y) |> 
    broom::tidy() |>
    select(estimate, conf.low, conf.high, p.value)
}

corr_data <- plot_data |> select(-n, -n_trt, -trt_prop_lb, -trt_prop_ub, -surv_lb, -surv_ub) 

corr_data |> 
  mutate(
    stage = case_when(
      stage == "all" ~ "All",
      stage == "L" ~ "1",
      stage == "R" ~ "3",
      stage == "D" ~ "4",
      TRUE ~ stage
    )
  ) |>
  mutate(stage = factor(
    as.character(stage), 
    levels = c("All", "1", "2", "3", "4")
  )
  ) |>
  mutate(stage = factor(
    stage, 
    labels = c("All", "1 or L", "2", "3 or R", "4 or D")
  )
  ) |>
  mutate(trt = str_to_lower(trt)) |>
  mutate(trt = paste0("Net survival vs ",trt)) |>
  filter(year != 1) |>
  group_by(trt, stage, year) |>
  summarise(
      estimate  = calc_cor(x = trt_prop, y = surv)$estimate,
      conf.low  = calc_cor(x = trt_prop, y = surv)$conf.low,
      conf.high = calc_cor(x = trt_prop, y = surv)$ conf.high,
      p.value   = calc_cor(x = trt_prop, y = surv)$p.value
    ) |>
  rename(compare = trt) |>
  saveRDS("results/corr_surv.RDS")

# Plot - all --------------------------------------------------------------------
p <- ggplot(
  data = plot_data |> filter(stage == "all") |> filter(year == 5), 
  aes(
      x = trt_prop,
      xmin = trt_prop_lb,
      xmax = trt_prop_ub,
      y = surv,
      ymin = surv_lb,
      ymax = surv_ub,
      colour = country
    )
  )
p <- p + geom_point() ##+ geom_errorbarh(height = 0) + geom_errorbar(width = 0)
p <- p + facet_wrap(~trt)
p <- p + scale_x_continuous(
  name = "Receiving treatment (%)",
  limits = c(0, 0.7),
  breaks = c(0, 0.2, 0.4, 0.6),
  labels = c("0%", "20%", "40%", "60%"),
  expand = c(0,0),
  minor_breaks = NULL
)
p <- p + scale_y_continuous(
  name = "Five-year net survival (%)",
  limits = c(9, 28),
  breaks = c(10, 15, 20, 25),
  labels = c("10%", "15%", "20%", "25%"),
  minor_breaks = NULL
)
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
p
p |> saveRDS("results/figure7a_trt_vs_surv_all.RDS")


# Plot - all - labelled --------------------------------------------------------------------
p <- ggplot(
  data = plot_data |> filter(stage == "all") |> filter(year == 5), 
  aes(
    x = trt_prop,
    xmin = trt_prop_lb,
    xmax = trt_prop_ub,
    y = surv,
    ymin = surv_lb,
    ymax = surv_ub,
    colour = jurisdiction, 
    shape = jurisdiction
  )
)
p <- p + geom_point() ##+ geom_errorbarh(height = 0) + geom_errorbar(width = 0)
p <- p + facet_wrap(~trt)
p <- p + scale_x_continuous(
  name = "Receiving treatment (%)",
  limits = c(0, 0.7),
  breaks = c(0, 0.2, 0.4, 0.6),
  labels = c("0%", "20%", "40%", "60%"),
  expand = c(0,0),
  minor_breaks = NULL
)
p <- p + scale_y_continuous(
  name = "Five-year net survival (%)",
  limits = c(9, 28),
  breaks = c(10, 15, 20, 25),
  labels = c("10%", "15%", "20%", "25%"),
  minor_breaks = NULL
)
p <- p +
  theme_bw() +
  theme_dotplots() + 
  theme(legend.position = "bottom") + 
  guides(
    alpha = "none",
    size = "none",
    fill = "none",
    #shape = "none"
  ) +
  scale_colour_manual(
    'Jurisdiction', 
    values = c(
      "#984ea3", "#984ea3", "#984ea3", "#984ea3", 
      "#e41a1c", 
      "#4daf4a", "#4daf4a", "#4daf4a", "#4daf4a", "#4daf4a", "#4daf4a", "#c5e6c4", "#c5e6c4", "#c5e6c4",
      "#377eb8", "#377eb8"
      )
    ) +
  scale_fill_manual(
    'Jurisdiction', 
    values = c(
      "#984ea3", "#984ea3", "#984ea3", "#984ea3", 
      "#e41a1c", 
      "#4daf4a", "#4daf4a", "#4daf4a", "#4daf4a", "#4daf4a", "#4daf4a", "#c5e6c4", "#c5e6c4", "#c5e6c4",
      "#377eb8", "#377eb8"
    )
  ) +
  scale_shape_manual(
    'Jurisdiction', 
    values = c(
      1, 2, 3, 4,
      1,
      1, 2, 3, 4, 5, 6, 1, 2, 3,
      1, 2))
p
p |> saveRDS("results/figure7a_trt_vs_surv_all_alt.RDS")
ggsave(
       "results/figure7a_trt_vs_surv_all_alt.png",
       plot = p,
       width = 15,
       height = 13,
       units = "cm")

# Plot - stage specific ---------------------------------------------------
p <- ggplot(
  data = plot_data |> 
    filter(stage != "all") |> 
    filter(year == 3) |>
    mutate(
      stage = case_when(
        stage == "L" ~ "1",
        stage == "R" ~ "3",
        stage == "D" ~ "4",
        TRUE ~ stage
      )
    ) |>
    mutate(stage = factor(
        as.character(stage), 
        levels = c("1", "2", "3", "4")
      )
    ) |>
    mutate(stage = factor(
        stage, 
        labels = c("1 or L", "2", "3 or R", "4 or D")
      )
    ) |>
    mutate(trt_prop    = 100*trt_prop) |>
    mutate(trt_prop_lb = 100*trt_prop_lb) |>
    mutate(trt_prop_ub = 100*trt_prop_ub) 
    , 
  aes(
    x = trt_prop,
    xmin = trt_prop_lb,
    xmax = trt_prop_ub,
    y = surv,
    ymin = surv_lb,
    ymax = surv_ub,
    colour = country
  )
)
p <- p + geom_point() ##+ geom_errorbarh(height = 0) + geom_errorbar(width = 0)
p <- p + facet_grid(stage~trt, scales = "free_y")
p <- p + scale_x_continuous(
  name = "Receiving treatment (%)",
  limits = c(0, 100),
  breaks = c(0, 25, 50, 75),
  expand = c(0, 0),
  minor_breaks = NULL
)
p <- p + scale_y_continuous(
  name = "Three-year net survival (%)",
  #limits = c(0, 100),
  #breaks = c(0, 25, 50, 75),
  #expand = c(0, 0),
  minor_breaks = NULL
)
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
p
p |> saveRDS("results/figure7b_trt_vs_surv_stagespec.RDS")

# Clean up ----------------------------------------------------------------
rm(list = ls())


