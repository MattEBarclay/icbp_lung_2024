

# expit function ----------------------------------------------------------
# inverse logistic function for use with proportions data
expit <- function(x) {
  
  x = exp(x)/(1+exp(x))
  
}



# Helpful functions for processing time-to-treat data ---------------------

# calculate group centiles
calculate_grouped_centiles <- function(data) {
  # prepare data
  dat0 <- data |> 
    group_by(cancer, trt, stage, ctry_order) |>
    arrange(cancer, trt, stage, ctry_order, days) |>
    # pull out 'overall' patient and treatment counts
    # at 'max' days, so should be all treated patients
    mutate(n1 = case_when(row_number() == n() ~ n, TRUE ~ 0)) |>
    mutate(n_trt1 = case_when(row_number() == n() ~ n_trt, TRUE ~ 0)) |>
    mutate(n_trt_bit = case_when(
      percentile == 0 ~ 0,
      TRUE ~ sum(n_trt1)/20
    )
    )|>
    group_by() 
  
  # group data
  dat1 <- dat0 |>
    group_by(cancer, trt, stage, groups) |>
    mutate(n_total = sum(n1)) |>
    mutate(n_trt_total = sum(n_trt1)) |>
    arrange(cancer, trt, stage, groups, days) |>
    select(cancer, trt, stage, ctry_order, ctry_class, groups, percentile, 
           days, n_total, n_trt_total, n_trt_bit)
  
  # calculate empirical centiles
  dat2 <- dat1 |>
    # calculate cumulative treated within grouped jurisdictions
    mutate(n_trt_new = cumsum(n_trt_bit)) |>
    mutate(pct_new = 100*n_trt_new/n_trt_total) |>
    # get rid of all but first '0' row
    filter(pct_new != 0 | row_number() == 1)
  
  # restructure to match original
  dat3 <- dat2 |>
    mutate(n = n_total) |>
    mutate(n_trt = n_trt_new) |>
    mutate(prop = n_trt/n) |>
    select(-percentile) |>
    rename(percentile = pct_new) |>
    group_by(
      cancer,
      trt,
      stage, 
      ctry_class,
      groups
    ) |>
    select(
      cancer,
      trt,
      stage, 
      ctry_class,
      groups,
      percentile,
      days,
      n,
      n_trt,
      prop
    )
  
  # return
  dat3
  
}

# improve group centiles
improve_centiles <- function(data) {
  data |>
    mutate(pct_new = 5*round(percentile/5,0)) |> 
    mutate(diff = abs(pct_new - percentile)) |>
    group_by(trt, stage, groups, pct_new) |>
    arrange(diff) |>
    slice(1) |>
    ungroup() |>
    mutate(percentile = pct_new) |>
    select(-diff, -pct_new)
}

# stuff to plot
filter_cumulative <- function(centile, data, site, plot_these, to_label, nudge) {
  data |> 
    filter(ctry_class %in% plot_these) |>
    filter(percentile <= centile) |>
    mutate(label = case_when(percentile == centile & ctry_class %in% to_label ~ ctry_class, TRUE ~ NA)) |>
    mutate(label_x = days+30) |>
    mutate(ynudge = case_when(ctry_class %in% nudge ~ 1, TRUE ~ 0))
}

# GT meta-analysis common -------------------------------------------------
# Common options for meta-analysis tables
gt_metan_common <- function(tabdata) {
  tabdata |>
    cols_merge(
      columns = c(lower, upper),
      pattern = "({1}, {2})"
    ) |>
    cols_merge(
      columns = c(pred_lb, pred_ub),
      pattern = md("{1}-{2}")
    ) |>
    cols_merge(
      columns = c(minimum, maximum),
      pattern = md("{1}-{2}")
    ) |>
    fmt_number(
      columns = c(i2),
      decimals = 1
    ) |>
    fmt_number(
      columns = c(tau),
      decimals = 3
    ) 
}


# Meta-analysis data formatting -------------------------------------------
# formatting of the poly_metan dataset that is consistent across analyses
common_meta_format <- function(data) {
  
  data |>
    pivot_longer(cols = starts_with("prop"), names_to = "prop", values_to = "x") |> 
    mutate(y = ifelse(prop == "prop0" | prop == "prop2", 0, ifelse(prop == "prop1", 0.3, ifelse(prop == "prop3", -0.3, NA)))) |>
    mutate(trt_order   = factor(trt, levels = c("Chemo", "Radio"))) 
  
}

# Functions for working with countries and stages -------------------------

assign_ctry_class <- function(data) {
  
  data |> mutate(
    ctry_class = case_when(
      ctry_order == "England"          ~ 1,
      ctry_order == "Northern Ireland" ~ 1,
      ctry_order == "Scotland"         ~ 1,
      ctry_order == "Wales"            ~ 1,
      ctry_order == "Norway"           ~ 2,
      ctry_order == "Alberta"          ~ 3,
      ctry_order == "British Columbia" ~ 3,
      ctry_order == "Ontario"          ~ 3,
      ctry_order == "Saskatchewan"     ~ 3,
      ctry_order == "Manitoba"         ~ 3,
      ctry_order == "Prince Edward Island"  ~ 3,
      ctry_order == "New Brunswick"    ~ 3,
      ctry_order == "Newfoundland & Labrador"     ~ 3,
      ctry_order == "Nova Scotia"      ~ 3,
      ctry_order == "New South Wales"  ~ 4,
      ctry_order == "Victoria"         ~ 4
    )
  ) |>
    mutate(
      ctry_class = factor(
        ctry_class, 
        levels = c(1, 2, 3, 4),
        labels = c("United Kingdom", "Norway", "Canada", "Australia")
      )
    )
  
}

group_canada <- function(data) {
  
  atlantic_canada <- c("Prince Edward Island", "Nova Scotia", "New Brunswick", "Newfoundland & Labrador")
  skmb <- c("Saskatchewan", "Manitoba")
  
  data |>
    mutate(jdn_class = case_when(
      as.character(ctry_order) %in% atlantic_canada ~ "Atlantic Canada",
      as.character(ctry_order) %in% skmb ~ "Saskatchewan-Manitoba",
      TRUE ~ as.character(ctry_order)
    )) |>
    mutate(chk = as.character(ctry_order)) |>
    mutate(
      ctry_order = factor(
        jdn_class, 
        levels = c(
          "England",
          "Northern Ireland",
          "Scotland",
          "Wales",
          "Norway",
          "Alberta",
          "British Columbia",
          "Ontario",
          "Saskatchewan-Manitoba",
          "Atlantic Canada",
          "New South Wales",
          "Victoria"
        )
      )
    )
  
}

# Rationalise stage
rationalise_stage <- function(data) {
  
  data |>
    mutate(variable_value = case_when(
      variable != "stage" ~ variable_value,
      variable_value == "I"   | variable_value == "1" ~ "1",
      variable_value == "II"  | variable_value == "2" ~ "2",
      variable_value == "III" | variable_value == "3" ~ "3",
      variable_value == "IV"  | variable_value == "4" ~ "4",
      variable_value == "Localised"     ~ "L",
      variable_value == "L"             ~ "L",
      variable_value == "Regional spread (adjacent organs)" ~ "R",
      variable_value == "Regional spread (lymph nodes)"     ~ "R",
      variable_value == "R"                                 ~ "R",
      variable_value == "Distant"         ~ "D",
      variable_value == "Distant&Unknown" ~ "D",
      variable_value == "M"               ~ "D",
      variable_value == "X" ~ "X",
      TRUE ~ variable_value
    ))
}

# Add variables for plotting
plotting_age <- function(data) {
  # add age factor
  data |>
    mutate(age = ifelse(variable == "age", variable_value, NA)) |>
    mutate(age = ifelse(age == "84&99", "84", age)) |>
    mutate(age = factor(age, levels = c(64, 74, 84, 99), labels = c("15-64", "65-74", "75-84", "85-99")))
}

plotting_year <- function(data) {
  # add diagnosis year
  data |>
    mutate(year = ifelse(variable == "diagnosis_year", variable_value, NA)) |>
    mutate(year = as.numeric(year))
}

plotting_stage <- function(data) { 
  # add stages
  data |>
    mutate(stage_x = ifelse(variable == "stage", variable_value, NA)) |>
    mutate(stage_x = ifelse(
      stage_x == "X" & (ctry_order == "Norway" | ctry_order == "New South Wales" | ctry_order == "Victoria"), 
      "X SEER",
      stage_x)
    ) |>
    mutate(stage_x = factor(stage_x, levels = c("1", "2", "3", "4", "X", "L", "R", "D", "X SEER")))
}

plotting_sex  <- function(data) {
  # add sex
  data |>
    mutate(sex = ifelse(variable == "sex", str_to_upper(variable_value), NA)) |>
    mutate(sex = factor(sex, levels = c("F", "M")))
}

plotting_variables <- function(data) {
  data |>
    plotting_age() |>
    plotting_year() |>
    plotting_stage() |>
    plotting_sex()
}

# Exclusion criteria
exclusion <- function(data) {
  data |> 
    filter(!(ctry_order == "Newfoundland & Labrador" & trt == "Chemo"))
}



# Very common plot options -----------------------------------
theme_dotplots <- function(labfill = "grey", x_angle = 0) {
  
  theme(
    strip.text.x=element_text(size=6, color="black", hjust=0 ),
    strip.background = element_rect(fill = labfill) ,
    axis.text.x=element_text(size=8, color="black", vjust = 0.5, hjust=0.5, angle = x_angle),
    axis.text.y=element_text(size=8, color="black"),
    axis.ticks.length = unit(0, "cm"),
    panel.spacing.x = unit(0.5, "lines") ,
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    panel.spacing.y = unit(0.3,"line"),
    panel.grid.major.y = element_blank()
  )
  
} 

theme_lineplots <- function(labfill = "grey", x_angle = 0) {
  
  theme(
    strip.text.x=element_text(size=8, color="black", hjust=0 ),
    strip.background = element_rect(fill = labfill),
    axis.text.x=element_text(size=7, color="black", vjust = 0.5, angle = x_angle),
    axis.text.y=element_text(size=7, color="black"),
    axis.ticks.length = unit(0, "cm"),
    panel.spacing.x = unit(0.5, "lines"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    panel.spacing.y = unit(0.3,"line")
  )
  
}


# Colour scheme
icbp_colours <- c("#524A6F", "#003566", "#18ACBD", "#BC1C80")

icbp_colour_manual <- function() {
  scale_color_manual(values = icbp_colours)
}

icbp_fill_manual <- function() {
  scale_fill_manual(values = icbp_colours)
}


# Function for trend plot ------------------------------------------------------
trend_plot <- function(
      df, 
      panel_title,
      variable_touse, 
      trt_touse, 
      stages, 
      x_var, 
      x_name, 
      x_break, 
      x_limit, 
      x_label, 
      x_angle, 
      y_name, 
      labfill = "grey", 
      facet_scales = "fixed", 
      stage_plot = FALSE
    ) {
  
  # Some specific fudges
  if(stage_plot) {
    background_lines <- df |> filter(stage == "All stages") |> filter(trt == trt_touse) |> filter(variable == "stage") |> mutate(grouping = ctry_order)
    countries  <- background_lines |> select(ctry_order) |> unique() |> pull()
    countries2 <- background_lines |> select(ctry_order) |> unique() |> pull()
    crossmap <- expand.grid(countries, countries2)
    crossmap <- crossmap |>
      rename(grouping = Var1) |>
      rename(ctry_order = Var2) |> 
      mutate(keep = case_when(
        ctry_order %in% c("Norway", "New South Wales", "Victoria") & grouping %in% c("Norway", "New South Wales", "Victoria") ~ 1,
        !(ctry_order %in% c("Norway", "New South Wales", "Victoria")) & !(grouping %in% c("Norway", "New South Wales", "Victoria")) ~ 1,
        TRUE ~ 0
      )) |>
      filter(keep == 1) |>
      select(grouping, ctry_order)
    
    background_lines <- background_lines |> select(-ctry_order)
    background_lines <- background_lines |> inner_join(crossmap)

  }
  if(!stage_plot) {
    background_lines <- df |> filter(stage == stages) |> filter(trt == trt_touse) |> filter(variable == variable_touse) |> mutate(grouping = ctry_order) |> mutate(ctry_order = NULL)
  }

  p <- ggplot()
  p <- p + geom_line(
    data = background_lines, 
    colour = "gray", 
    aes(
      y = prop, 
      x = as.numeric({{x_var}}),
      group = grouping
    )
  )
  p <- p + geom_line(
    data = df |> filter(stage == stages) |> filter(trt == trt_touse) |> filter(variable == variable_touse), 
    size = 1, 
    aes(
      y = prop, 
      x = as.numeric({{x_var}}),
      group = ctry_order,
      colour = ctry_class
    )
  )
  p <- p + geom_line(
    data = df |> filter(stage == stages) |> filter(trt == trt_touse) |> filter(variable == variable_touse), 
    linewidth = 0.5, 
    linetype = 2,
    aes(
      y = lower, 
      x = as.numeric({{x_var}}),
      group = ctry_order,
      colour = ctry_class
    )
  )
  p <- p + geom_line(
    data = df |> filter(stage == stages) |> filter(trt == trt_touse) |> filter(variable == variable_touse), 
    linewidth = 0.5, 
    linetype = 2,
    aes(
      y = upper, 
      x = as.numeric({{x_var}}),
      group = ctry_order,
      colour = ctry_class
    )
  )
  p <- p + geom_point(
    data = df |> filter(stage == stages) |> filter(trt == trt_touse) |> filter(variable == variable_touse), 
    size = 1,
    aes(
      y = prop, 
      x = as.numeric({{x_var}}),
      group = ctry_order,
      colour = ctry_class
    )
  )
  p <- p + scale_x_continuous(
    name = x_name,
    limits = x_limit,
    breaks = x_break,
    minor_breaks = NULL,
    labels = x_label
  )
  p <- p + scale_y_continuous(
    limits = c(0, 0.75),
    breaks = c(0, 0.25, 0.5), 
    labels = c("0%", "25%", "50%"), 
    minor_breaks = c(0.125, 0.375, 0.625), 
    expand = c(0,0), 
    name = y_name
  )
  p <- p +
    theme_bw() +
    theme_lineplots(labfill = labfill, x_angle = x_angle) +
    guides(
      alpha = "none",
      size = "none",
      fill = "none",
      colour = "none",
      shape = "none"
    ) 
  p <- p + scale_colour_manual(values = c("#984ea3", "#e41a1c", "#4daf4a", "#377eb8"))
  p <- p + scale_fill_manual(values = c("#984ea3", "#e41a1c", "#4daf4a", "#377eb8"))
  p <- p + facet_wrap(ctry_order~., ncol = 4, scales = facet_scales)
  p <- p + labs(title = panel_title)
  p <- p + theme(plot.title = element_text(size = 12))
  p <- p + theme(strip.text.x = element_text(size = 6))
  p 
}



# Function for time-trend plot --------------------------------------------
dot_plot <- function(
    df, 
    trt_touse, 
    x_name, 
    labfill
) {
  
  df <- df |> filter(stage %in% c("All stages", "Stages 1-3", "Stage 4"))
  df <- df |> mutate(stage = as.character(stage)) |>
    mutate(stage = factor(stage, levels = c("All stages", "Stages 1-3", "Stage 4"))) |>
    mutate(stage = factor(stage, labels = c("All stages", "Stages 1-3 or L-R", "Stage 4 or distant")))
  
  
  p <- ggplot(data = df |> filter(trt == trt_touse))
  # 5th to 95th percentiles
  # p <- p + geom_errorbarh(
  #   colour = "black",
  #   linetype = 1,
  #   size = 0.3,
  #   height = 0,
  #   aes(
  #     xmin = pct_5,
  #     xmax = pct_95,
  #     y    = as.numeric(ctry_order)
  #   )
  # )
  # 25th to 75th percentiles
  p <- p + geom_errorbarh(
    linetype = 1,
    size = 1,
    height = 0.1,
    aes(
      xmin = pct_25,
      xmax = pct_75,
      y    = as.numeric(ctry_order),
      colour = ctry_class
    )
  )
  # median
  p <- p + geom_point(
    aes(
      x = pct_50,
      y = as.numeric(ctry_order),
      colour = ctry_class
    )
  )
  # X scale - days
  p <- p + scale_x_continuous(
    name = x_name,
    # limits = c(-5, 365),
    # breaks = c(0,90,180,270,365),
    # minor_breaks = NULL,
    # labels = c("0","90","180","270","365")
    limits = c(0, 165),
    breaks = c(0,30,60,90,120),
    minor_breaks = NULL
  )
  # Y scale - countries
  p <- p + scale_y_continuous(
    breaks = 1:16,
    labels = country_order,
    name = "Jurisdiction",
    trans = "reverse"
  )
  # Theme
  p <- p +
    theme_bw() +
    theme_dotplots(labfill = labfill) +
    guides(
      alpha = "none",
      size = "none",
      fill = "none",
      colour = "none",
      shape = "none"
    )
  p <- p + icbp_colour_manual()
  p <- p + icbp_fill_manual()
  # One column for each stage
  p <- p + facet_wrap(stage~., ncol = 4)
  p 
}

# Correlation plot --------------------------------------------------------
create_missing_data_frame <- function(df) {
  df |> 
    filter(stage == "All stages") |> 
    filter(variable == "stage") |>
    select(trt, ctry_class, ctry_order, variable_value, n, n_trt) |>
    mutate(variable_value = case_when(
      variable_value != "X" ~ "known",
      variable_value == "X" ~ "X"
    )) |>
    group_by(trt, ctry_class, ctry_order, variable_value) |>
    summarise_all(sum) |> 
    pivot_wider(id_cols = c(trt, ctry_class, ctry_order), values_from = c(n, n_trt), names_from = variable_value, values_fill = 0) |>
    mutate(X_pct = n_X/(n_X+n_known)) |>
    select(trt, ctry_class, ctry_order, n_X, n_trt_X, X_pct) |>
    mutate(X_trt_pct = n_trt_X/n_X) |>
    select(trt, ctry_class, ctry_order, X_trt_pct, X_pct) |>
    filter(ctry_order != "Victoria")
}


corr_plot <- function(
    df,
    trt_touse
) {
  
  X_corr2 <- df |> create_missing_data_frame()
  
  X_corr2 <- X_corr2 |> filter(trt == trt_touse)
  
  r <- lm(X_trt_pct ~ X_pct, data = X_corr2) |> glance() |> select(r.squared) |> sqrt() |> round(digits = 2) |> pull()
  yname <- str_to_lower(trt_touse)
  
  p <- ggplot(data = X_corr2, aes(x = X_pct, y = X_trt_pct, label = ctry_order, colour = ctry_class))
  p <- p + geom_point()
  p <- p + geom_text(size = 2, nudge_y = 0.01)
  p <- p + annotate(
    "text", label = paste("R =", r),
    x = 0.03, y = 0.48, size = 5, colour = "black"
  )
  p <- p + scale_x_continuous(
    name = "% Lung cancer with unknown stage",
    limits = c(0, 0.55),
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    minor_breaks = NULL,
    labels = c("0%", "10%", "20%", "30%", "40%", "50%")
  )
  p <- p + scale_y_continuous(
    name = paste0("% of unknown stage receiving ", yname),
    limits = c(0, 0.5),
    breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
    minor_breaks = NULL,
    labels = c("0%", "10%", "20%", "30%", "40%", "50%")
  )
  p <- p +
    theme_bw() +
    theme_dotplots() +
    guides(
      alpha = "none",
      size = "none",
      fill = "none",
      colour = "none",
      shape = "none"
    ) 
  p <- p + icbp_colour_manual()
  p <- p + icbp_fill_manual()
  p
}


# Treatment specific stage correlations -----------------------------------
trt_spec_stage_corrs <- function(df, trt_touse, trt_name, labfill) {
  
  p <- ggplot(data = df |> filter(trt == trt_touse), aes(x = stage13, y = stage4, label = ctry_order, colour = ctry_class))
  p <- p + geom_point()
  p <- p + scale_x_continuous(
      name = paste0("Median time to ", trt_name, ", stages 1-3 or L-R"),
      limits = c(20, 100),
      breaks = c(20, 40, 60, 80, 100),
      minor_breaks = NULL
    )
  p <- p + scale_y_continuous(
      name = paste0("Median time to ", trt_name, ", stage 4 or distant"),
      limits = c(20, 100),
      breaks = c(20, 40, 60, 80, 100),
      minor_breaks = NULL
    )
  p <- p +
    theme_bw() +
    theme_dotplots() +
    guides(
      alpha = "none",
      size = "none",
      fill = "none",
      colour = "none",
      shape = "none"
    ) +
    icbp_colour_manual() +
    icbp_fill_manual()
  p
}


# Use vs time -------------------------------------------------------------
plot_use_v_time <- function(data, x_var = dayChemo, y_var = trtChemo, trt_name = "chemotherapy", labfill = "grey") {
  p <- ggplot(data = data, aes(x = {{x_var}}, y = {{y_var}}, label = ctry_order, colour = ctry_class))
  p <- p + geom_point()
  p <- p + scale_y_continuous(
    name = paste0("Receiving ", trt_name, " (%)"),
    limits = c(0.2, 0.6),
    breaks = c(0.2, 0.3, 0.4, 0.5, 0.6),
    labels = c("20%", "30%", "40%", "50%", "60%"),
    minor_breaks = NULL
  )
  p <- p + scale_x_continuous(
    name = paste0("Median time to ", trt_name),
    limits = c(20, 100),
    breaks = c(20, 40, 60, 80, 100),
    minor_breaks = NULL
  )
  p <- p + facet_wrap(~stage)
  p <- p +
    theme_bw() +
    theme_dotplots(labfill = labfill) + 
    guides(
      alpha = "none",
      size = "none",
      fill = "none",
      colour = "none",
      shape = "none"
    ) +
    icbp_colour_manual() +
    icbp_fill_manual()
  p
}

# Included Jdns -----------------------------------------------------------
# Plot jurisdiction-specific results for the jdns used in the meta-analysis
incl_error <- function(yvar, data = meta_dat) {
  geom_errorbarh(
    data = data,
    height = .1,
    aes(y = {{yvar}}, xmin = lower, xmax = upper, colour = ctry_class)
  )
}
incl_point <- function(yvar, xvar, data = meta_dat) {
  geom_point(
    data = data,
    aes(y = {{yvar}}, x = {{xvar}}, colour = ctry_class, fill = include),
    shape = 21
  )
}

# Meta analysis diamond ---------------------------------------------------
# Plot the meta-analysis summaries
meta_error <- function(y_shift) {
  geom_errorbarh(
    data = poly_metan |> 
      filter( y == 0) |>
      mutate(y = y + {{y_shift}}) |> 
      mutate(
        jurisdiction = "RE meta-analysis",
        ctry_order = factor(jurisdiction, levels = country_order_chm)
      ),
    height=.3,
    colour = "darkgrey",
    linetype = "solid",
    aes(y = y, xmin = pred_lb, xmax = pred_ub)
  ) 
}
meta_diamond <- function(y_shift, groupvar) {
  geom_polygon(
    data = poly_metan |>
      mutate(y = y + {{y_shift}}) |>
      mutate(
        jurisdiction = "RE meta-analysis",
        ctry_order = factor(jurisdiction, levels = country_order_chm)
      ),
    colour = "darkgrey",
    fill = "darkgrey",
    aes(
      x = x, 
      y = y,
      group = trt,
      subgroup = {{groupvar}}
    )
  )
}

# Standard y-axis for jurisdictions
jdn_yscale <- function() {
  scale_y_continuous(
    name = "", 
    #limits = c(0. 20), 
    breaks = c(1:5, 7:15, 17:18, 20), 
    minor_breaks = NULL, 
    labels = c(
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
      "Victoria",
      "Pooled estimate"
    ), 
    #labels = blank_order,
    trans="reverse"
  ) 
}



theme_icbp <- function() {
  
  theme(
    strip.text.x = element_text(size=8, color="black", hjust=0 ), 
    axis.text.x  = element_text(size=7, color="black"),
    axis.text.y  = element_text(size=7, color="black"),
    axis.title   = element_text(size=8, color="black"),
    axis.ticks.length = unit(0, "cm"),
    panel.spacing.x = unit(0.5, "lines"),
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    panel.spacing.y = unit(0.3,"line") ,
    legend.position = "none",
    panel.grid.major.y = element_blank()
  ) 
  
}

# Colour scheme
icbp_colours <- c("#003566", "#524A6F", "#18ACBD", "#BC1C80")

icbp_colour_manual <- function() {
  scale_color_manual(values = icbp_colours)
}

icbp_fill_manual <- function() {
  scale_fill_manual(values = icbp_colours)
}



