#!/usr/bin/env Rscript

########################
## CONFIG + PACKAGES  ##
########################

##-------------------------------------------------------------------------
## Configs + packages

library(tidyverse)
library(data.table)
library(magrittr)
library(ggplotify)
library(ggforce)
library(packedcircles)
library(cowplot)

source("../config/config.R")

###########################################################################
###########################################################################
###                                                                     ###
###                            FUNCTIONS                                ###
###                                                                     ###
###########################################################################
###########################################################################

### Change variables for plotting / tables.
configure_plotting_variables <- function(meta, colours_order_df) {
  meta %>%
    mutate(#### For niceties when plotting the data
          family_group_type_broad = case_when(family_group_type == "Singleton" ~ "Singleton",
                                              family_group_type == "Duo with Mother or Father" ~ "Duo",
                                              family_group_type == "Trio with Mother and Father" ~ "Trio",
                                              TRUE ~ "Other"),
          family_group_type_broad = factor(family_group_type_broad, levels = c("Singleton", "Duo", "Trio", "Other")),
          family_group_type_broad_long = case_when(family_group_type_broad == "Singleton" ~ "Singleton probands",
                                                   family_group_type_broad == "Duo" ~ "Duos (with Mother or Father)",
                                                   family_group_type_broad == "Trio" ~ "Trios (with Mother and Father)",
                                                   TRUE ~ "Other"),
          family_group_type_broad_long = factor(family_group_type_broad_long, levels = c("Singleton probands",
                                                                                    "Duos (with Mother or Father)",
                                                                                    "Trios (with Mother and Father)",
                                                                                    "Other")),
          penetrance = ifelse(family_group_type == "Singleton", "n/a", penetrance),
          Penetrance = str_to_title(penetrance),
          Penetrance = factor(Penetrance, levels = c("N/A", "Complete", "Incomplete")),
          pen_show = case_when(penetrance == "complete" ~ "Comp.",
                               penetrance == "incomplete" ~ "Incomp.",
                               TRUE ~ "n/a"),
          participant_genetic_category = factor(participant_genetic_category,
                                                levels = colours_order_df$participant_genetic_category))
}

## Plotting function. Will work for any grouping label if desired.
## But make sure to have the correct colour/order in an aesthetics dataframe...
gmc_lat_long_plot <- function(df, gmc_lat_long, group_label, colours_order_df) {
  group_label_sym <- rlang::sym(group_label)
  gmc_prop_df <- df %>%
    group_by({{group_label_sym}}, handling_gmc_trust) %>%
    tally() %>%
    left_join(gmc_lat_long, by = "handling_gmc_trust") %>%
    group_by({{group_label_sym}}) %>%
    mutate(`Proportion (%)` = (n / sum(n)) * 100,
           {{group_label_sym}} := factor({{group_label_sym}},
                                         levels = pull(colours_order_df,
                                                       {{group_label_sym}})))

  world_sf <- ne_countries(scale = "medium", returnclass = "sf")
  uk_sf <- world_sf[which(world_sf$name == "United Kingdom"),]

  ggplot(data = uk_sf) +
    geom_sf(fill = "grey90", colour = "grey95") +
    theme_void(base_size = 6.5) +
    theme(legend.position = "top") +
    geom_point(data = gmc_prop_df,
               aes(x = longitude,
                   y = latitude,
                   size = `Proportion (%)`,
                   fill = {{group_label_sym}},
                   colour = {{group_label_sym}}),
               alpha = 0.7) +
    scale_size_area() +
    guides(colour = "none", fill = "none") +
    coord_sf(clip = 'off') +
    facet_wrap(as.formula(paste("~", group_label)), nrow = 3) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}}))
}

#### For Supplementary Figure 1
age_ethn_trend_plot <- function(df, group_label,
                                datasets = c("ONS 2021", "100KGP Probands"),
                                colours_order_df) {

  group_label_sym <- rlang::sym(group_label)
  ons_rd_ethn_age_props <- df %>%
    group_by(age_bin, group, {{group_label_sym}}) %>%
    summarise(num_inds = sum(num_inds)) %>%
    group_by(group, age_bin) %>%
    mutate(prop = num_inds / sum(num_inds)) %>%
    mutate({{group_label_sym}} := factor({{group_label_sym}}, levels = unique(pull(colours_order_df,
                                                                              {{group_label_sym}}))),
           group = factor(group, levels = datasets))

    ons_rd_ethn_age_props %>%
      ggplot(aes(x = age_bin, y = prop * 100,
             colour = {{group_label_sym}},
             shape = group,
             lineype = group,
             group = group)) +
      geom_point() +
      geom_line(aes(linetype = group)) +
      facet_wrap(as.formula(paste0("~", group_label)), scale = "free_y") +
      theme_classic(base_size = 6.5) +
      guides(colour = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks = element_blank(),
            axis.line = element_blank(),
            strip.background = element_blank(),
            legend.position = "top",
            legend.title = element_blank(),
            panel.grid.major.y = element_line(color = "grey90", linewidth = .6)) +
      xlab("") + ylab("Proportion (%)") +
      coord_cartesian(clip = 'off') +
      scale_colour_manual(values = pull(colours_order_df, colour),
                          limits = pull(colours_order_df, {{group_label_sym}}))
}

# For Figure 2 etc.
errorbar_plot <- function(df, group_label, reference_group = NULL,
                          colours_order_df, p_label = "p.value") {
    group_label_sym <- rlang::sym(group_label)
    p_label_sym <- rlang::sym(p_label)
    ### Add stuff
    df %<>%
    mutate(sig_stars = case_when(
        {{p_label_sym}} < 0.1  & {{p_label_sym}} >= 0.05  ~ ".",
        {{p_label_sym}} < 0.05 & {{p_label_sym}} >= 0.01 ~ "*",
        {{p_label_sym}} < 0.01 & {{p_label_sym}} >= 0.001 ~ "**",
        {{p_label_sym}} < 0.001 ~ "***",
        is.na({{p_label_sym}}) ~ NA_character_,
        TRUE ~ NA_character_),
        {{group_label_sym}} := factor({{group_label_sym}},
            levels = rev(pull(colours_order_df, {{group_label_sym}}))))
    
    df %>% filter({{p_label_sym}} < 0.05) %>%
      pull({{group_label_sym}}) -> sig_pops
    colours_order_sig <- colours_order_df %>%
      mutate(colour = ifelse({{group_label_sym}} %in% sig_pops,
                             colour, "white"))

    ## Adjust top group if requested
    if(!is.null(reference_group)) {
      df %<>%
        mutate({{group_label_sym}} := relevel({{group_label_sym}},
              ref = reference_group))
    }
    ### Plot
    ggplot(df) +
    aes(x = estimate, y = {{group_label_sym}},
        colour = {{group_label_sym}},
        fill = {{group_label_sym}}) +
    geom_vline(xintercept = 0, linewidth = 0.7, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(xmin = lower, xmax = upper), width = 0.1) +
    geom_point(shape = 21) +
    geom_text(data = df,
        aes(x = estimate, y = {{group_label_sym}},
            label = sig_stars), vjust = -0.05) +
    theme_minimal() +
    coord_cartesian(clip = 'off') +
    ylab("") +
    theme(legend.position = "NONE",
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    scale_fill_manual(values = pull(colours_order_sig, colour),
                      limits = pull(colours_order_sig, {{group_label_sym}}))
}

admixture_plot <- function(df, group_label,
                           order_rows_df = NULL, colours_order_df) {
  group_label_sym <- rlang::sym(group_label)
  adm_prop <- df %>%
    arrange(desc(`Proportional similarity to Europe North West`)) %>%
    dplyr::select(participant_id,
                  platekey,
                  {{group_label_sym}},
                  starts_with("Proportional similarity to")) %>%
    pivot_longer(cols = starts_with("Proportional similarity to"),
                 names_to = "Assignment",
                 values_to = "alpha") %>%
    mutate(Assignment = gsub("Proportional similarity to ",
                           "",
                           Assignment)) %>%
    adjust_genetic_labels(collapse_small = FALSE) %>%
    mutate(participant_genetic_category := factor(participant_genetic_category,
                                             levels = unique(pull(colours_order_df,
                                                             participant_genetic_category))),
           participant_id = as.character(participant_id),
           participant_id = factor(participant_id, levels = unique((.$participant_id))))
  
  if (!is.null(order_rows_df)) {
    adm_prop %<>% mutate({{group_label_sym}} := factor({{group_label_sym}},
                                             levels = unique(pull(order_rows_df,
                                                                  {{group_label_sym}}))))
  }

  ggplot(adm_prop) +
    aes(x = participant_id,
        y = alpha,
        group = Assignment,
        colour = Assignment,
        fill = Assignment) +
    geom_bar(stat = "identity") +
    facet_wrap(as.formula(paste0("~", {{group_label_sym}})),
               scales = "free_x",
               ncol = 1, strip.position = "left") +
    theme_minimal() +
    theme(legend.position = "bottom",
          strip.text.y.left = element_text(angle = 0, hjust=1),
          axis.text.x = element_blank(),
          strip.background = element_blank(),
          strip.placement = "outside",
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, participant_genetic_category)) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, participant_genetic_category)) +
    xlab("") + ylab("Alpha") +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    theme(plot.margin = margin(t = 0, r = 100, b = 0, l = 0, "points"))
}

### Will also work for any concordance between two variables...
concordance_plot <- function(df, group_label, cat_var, order_rows_df,
                             colours_order_df, hide_val) {
    cat_var_sym <- rlang::sym(cat_var)
    group_label_sym <- rlang::sym(group_label)
    df %<>%
      mutate(proportion = ifelse(n <= hide_val, NA, proportion),
             {{group_label_sym}} := factor({{group_label_sym}},
                                               levels = unique(pull(order_rows_df,
                                                                   {{group_label_sym}}))),
             {{cat_var_sym}} := factor({{cat_var_sym}},
                                              levels = unique(pull(colours_order_df,
                                                                   {{cat_var_sym}}))))
    ggplot(df) +
    aes(x = {{cat_var_sym}}, y = proportion * 100,
        colour = {{cat_var_sym}}, fill = {{cat_var_sym}}) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste(round((proportion * 100), 1), "%", sep = "")),
                  vjust = -0.5, size = 2.0) +
    facet_wrap(as.formula(paste("~", group_label)),
        nrow = length(unique(pull(df, {{group_label_sym}}))),
        scales = "fixed", strip.position = "left") +
    theme_minimal() +
    coord_cartesian(clip = 'off') +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey40") +
    theme(legend.position = "NONE",
        strip.text.y.left = element_text(angle = 0, hjust=1),
        axis.text.y = element_text(size = 6.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.background = element_blank(),
        strip.placement = "outside",
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    ylab("") + xlab("") +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{cat_var_sym}})) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, {{cat_var_sym}}))
}

### For Figure 1a
zoomed_bar_plot <- function(df, group_label,
                           colours_order_df,
                           zoom = c(22300, 29100)) {

  group_label_sym <- rlang::sym(group_label)
  df_counts <- df %>%
    group_by({{group_label_sym}}) %>%
    tally() %>%
    arrange(n) %>%
    mutate({{group_label_sym}} := factor({{group_label_sym}},
                                         levels = unique(pull(., {{group_label_sym}})))) 

  df_counts %>%
    ggplot(aes(y = "", x = n,
               fill = {{group_label_sym}})) +
    geom_bar(stat = "identity", colour = "white") +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, {{group_label_sym}})) +
    ylab("") +
    facet_zoom(xlim = zoom)
}

#### For Figure 1b
pca_plot_supp <- function(
    x, axis1, axis2, group_label, colours_order_df
    ) {
    group_label_sym <- rlang::sym(group_label)
    axis1_sym <- rlang::sym(axis1)
    axis2_sym <- rlang::sym(axis2)
    ggplot(x) +
        aes(x = {{axis1_sym}}, y = {{axis2_sym}},
            colour = {{group_label_sym}}) +
        geom_point(size = 0.2, alpha = 0.5) +
        theme_classic(base_size = 20) +
        xlab(axis1) + ylab(axis2) +
        scale_colour_manual(values = pull(colours_order_df, colour),
                            limits = pull(colours_order_df, {{group_label_sym}}))
}

raincloud_plot <- function(df, cont_var, group_label, colours_order_df) {
  cont_var_sym <- rlang::sym(cont_var)
  group_label_sym <- rlang::sym(group_label)
  df %<>%
    mutate({{group_label_sym}} := factor({{group_label_sym}},
                                          levels = colours_order_df %>%
                                                    pull({{group_label_sym}})))
  df %>%
    ggplot(aes(x = {{group_label_sym}},
               y = {{cont_var_sym}},
               colour = {{group_label_sym}},
               fill = {{group_label_sym}})) +
    theme_minimal() +
    theme(axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "NONE",
          plot.title = element_text(hjust = 0.5, size = 10)) +
    geom_boxplot(width = .2, fill = "white",
                 outlier.shape = NA) +
    ggdist::stat_halfeye(adjust = .75, ## bandwidth
                         width = .67,
                         color = NA,
                         position = position_nudge(x = .15)) +
    gghalves::geom_half_point(side = "l", 
                              range_scale = .3, 
                              alpha = .07, size = 0.1) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, {{group_label_sym}}))
}

### Cumulative ROH
box_facet_plot <- function(df, cont_var, group_label, colours_order_df) {
  cont_var_sym <- rlang::sym(cont_var)
  group_label_sym <- rlang::sym(group_label)
  df %>%
    mutate(row = "one") %>%
    ggplot(aes(x = {{cont_var_sym}},
               y = row,
               colour = {{group_label_sym}})) +
    theme_classic() +
    facet_wrap(as.formula(paste0("~", {{group_label_sym}})), ncol = 1, strip.position = "left") +
    theme(#axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          strip.text.y.left = element_text(angle = 0, hjust=1, size = 13),
          strip.background = element_blank(),
          strip.placement = "outside",
          plot.title = element_text(hjust = 0.5),
          legend.position = "NONE") +
    geom_boxplot(outlier.size = 0.1) +
    scale_colour_manual(values = pull(prive_colours, colour),
                        limits = pull(prive_colours, {{group_label_sym}}))
}

distribution_facet_plot <- function(
    df, stats_df, group_label, cont_var, type = "histogram", colours_order_df
    ) {
    cont_var_sym <- rlang::sym(cont_var)
    group_label_sym <- rlang::sym(group_label)
    ggObj <- ggplot(df) +
    aes(x = {{cont_var_sym}},
        fill = {{group_label_sym}},
        colour = {{group_label_sym}})
    if (type == "histogram") {
        ggObj <- ggObj + geom_histogram(aes(y = ..density.., group = 1), bins = 30)
    } else {
        ggObj <- ggObj + geom_density()
    }
    yinter <- median(layer_scales(ggObj)$y$range$range)
    ggObj <- ggObj + facet_wrap(as.formula(paste("~", group_label)),
        nrow = length(unique(pull(df, {{group_label_sym}})))) +
    geom_vline(data = stats_df,
        aes(xintercept = mean), linetype = "dashed") +
    theme_classic() +
    theme(legend.position = "NONE",
        axis.line = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5))
    if (type == "date") {
        ggObj <- ggObj + geom_text(data = stats_df,
            aes(x = mean,
            y = yinter,
            label = as.Date(mean, format = "%m/%y")),
            size = 3.0, hjust = -0.1,
            colour = "grey8")
    } else {
        ggObj <- ggObj + geom_text(data = stats_df,
            aes(x = mean,
            y = yinter,
            label = round(mean, 3),
            hjust = -0.1), colour = "grey8", size = 3.0)
    }
    ggObj <- ggObj +
      scale_colour_manual(values = pull(colours_order_df, colour),
                          limits = pull(colours_order_df, {{group_label_sym}})) +
      scale_fill_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}}))
}

proportion_facet_plot <- function(df, cat_var, group_label, colours_order_df) {
    group_label_sym <- rlang::sym(group_label)
    cat_var_sym <- rlang::sym(cat_var)
    paste(cat_var_sym)
    ggplot(df) +
    aes(x = {{cat_var_sym}}, y = proportion * 100,
        colour = {{group_label_sym}}, fill = {{group_label_sym}}) +
    geom_bar(stat = "identity") +
        geom_text(aes(label = paste(
            round((proportion * 100), 2), "%", sep = "")),
            vjust = -0.5, size = 3.0, colour = "grey8") +
    facet_wrap(as.formula(paste("~", group_label)),
        nrow = length(unique(pull(df, {{group_label_sym}}))),
        scales = "fixed") +
    theme_classic() +
    theme(legend.position = "NONE",
          axis.line = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 50, 100)) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, {{group_label_sym}}))
}

## Only valid for two level hierarchies
hierarchical_circles_plot <- function(df, hierarchy = c("phenotype_broad", "phenotype"),
                                     group_label = "participant_genetic_category",
                                     label_level = 1, label_as = "number",
                                     colours_order_df) {
  
  level1_sym <- rlang::sym(hierarchy[1])
  level2_sym <- rlang::sym(hierarchy[2])
  group_label_sym <- rlang::sym(group_label)
  
  hierarchies_df <- df %>%
    group_by({{level1_sym}}, {{level2_sym}}, {{group_label_sym}}) %>%
    summarise(size = n()) %>%
    ungroup() %>%
    mutate(name = paste0({{level2_sym}}, " (", {{group_label_sym}}, ")"),
           number = rownames(.))
  
  fill_categories <- hierarchies_df %>%
    left_join(colours_order_df, by = group_label) %>%
    pull(colour)
  
  vertices_df <- rbind(hierarchies_df %>%
                       dplyr::select(name, number, size),
                     hierarchies_df %>%
                       distinct({{level2_sym}}) %>%
                       rename(name = {{level2_sym}}) %>%
                       mutate(size = 0, number = rownames(.)),
                     hierarchies_df %>%
                       distinct({{level1_sym}}) %>%
                       rename(name = {{level1_sym}}) %>%
                       mutate(size = 0, number = rownames(.)) %>%
                     add_row(name = "all", size = 0, number = "1"))
  
  edges_df <- rbind(hierarchies_df %>%
                    distinct({{level2_sym}}, name) %>%
                    rename(from = {{level2_sym}}, to = name),
                  hierarchies_df %>%
                    distinct({{level1_sym}}, {{level2_sym}}) %>%
                    rename(from = {{level1_sym}}, to = {{level2_sym}}),
                  hierarchies_df %>%
                    distinct({{level1_sym}}) %>%
                    mutate(from = "all") %>%
                    rename(to = {{level1_sym}}))
  
  tree <- data.tree::FromDataFrameNetwork(edges_df)
  mylevels <- data.frame(name = tree$Get('name'),
                       level = tree$Get("level"))
  vertices_df %<>% left_join(mylevels, by=c("name" = "name")) 
  if (label_as == "number") {
    vertices_df %<>% mutate(new_label = ifelse(level == 2, number, NA))
  } else if (label_as == "name"){
    vertices_df %<>% mutate(new_label = ifelse(level == 2, name, NA))
  }

  hierarchies_graph <- graph_from_data_frame(edges_df, vertices = vertices_df)
  fill_all <- c(fill_categories,
                rep("white", nrow(vertices_df) - length(fill_categories)))
  
  circles_gg <- ggraph(hierarchies_graph, layout = 'circlepack', weight = size) + 
    geom_node_circle(aes(fill = as.factor(name), color = as.factor(depth), group = depth), size = 0.1) +
    scale_fill_manual(values = fill_all,
                      limits = as.factor(V(hierarchies_graph)$name)) +
    scale_color_manual(values=c("0" = "white", "1" = "black", "2" = "black", "3" = "white", "4" = "white")) +
    geom_node_label(aes(label = new_label), size = 2) +
    theme_void() +
    coord_fixed() +
    theme(legend.position="FALSE")
  return(circles_gg)
}

relate_popsizes_plot <- function(df, colours_order_df) {
  relate_popsizes <- df %>%
    filter(participant_genetic_category == participant_genetic_category_two) %>%
    mutate(participant_genetic_category = factor(participant_genetic_category,
                                                 levels = colours_order_df$participant_genetic_category))
  relate_popsize_background <- relate_popsizes %>%
    rename(participant_genetic_category2 = participant_genetic_category)

  relate_popsizes %>%
    ggplot(aes(x = log10(time), y = log10(pop_size),
               colour = participant_genetic_category,
               group = participant_genetic_category)) +
    geom_step(data = relate_popsize_background, aes(x = log10(time), y = log10(pop_size),
              group = participant_genetic_category2),
              colour = "grey80", alpha = 0.4) +
    geom_step() +
    theme_minimal() +
    theme(axis.ticks = element_blank(),
          strip.background = element_blank(),
          legend.position = "NONE") +
    facet_wrap(~participant_genetic_category, nrow = 3) +
    theme(legend.position = "NONE") +
    scale_x_continuous(limits = c(3, 6),
                       breaks = c(3, 4, 5, 6),
                       labels = scales::math_format(10^.x)) +
    scale_y_continuous(limits = c(3, 5),
                       breaks = c(3, 4, 5),
                       labels = scales::math_format(10^.x)) +
    annotation_logticks(size = 0.5,
                        colour = "grey70",
                        short = unit(0.025, "cm"),
                        mid = unit(0.05, "cm"),
                        long = unit(0.075, "cm")) +
    xlab("Years Before Present") +
    ylab("Ne") +
    coord_cartesian(clip = 'off') +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, participant_genetic_category))
}

compare_coeffs_plot <- function(df, group_label, corr_p, colours_order_df) {
  group_label_sym <- rlang::sym(group_label)
  gg <- df %>%
    mutate(variable = "dummy") %>%
    ggplot(aes(x = model1, y = model2, colour = {{group_label_sym}}, group = variable)) +
    geom_abline(slope=1, linetype = "dotted", colour = "grey80", linewidth = 0.75) +
    geom_point() + theme_minimal() +
    geom_errorbar(aes(xmax = upper951, xmin = lower951), alpha = 0.3) +
    geom_errorbar(aes(ymax = upper952, ymin = lower952), alpha = 0.3) +
    theme(panel.grid.minor = element_blank(), legend.title = element_blank()) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(linetype = 0)))
  if (corr_p < 0.01) {
    gg <- gg +
      geom_smooth(method = "lm", linetype = "dashed", colour = "grey60", se = FALSE, linewidth = 0.65)
  }
  return(gg)
}

raincloud_plot2 <- function(df, cont_var, group_label, colours_order_df) {
  cont_var_sym <- rlang::sym(cont_var)
  group_label_sym <- rlang::sym(group_label)
  df %<>%
    mutate({{group_label_sym}} := factor({{group_label_sym}},
                                          levels = rev(colours_order_df %>%
                                                    pull({{group_label_sym}}))))
  df %>%
    ggplot(aes(x = {{group_label_sym}},
               y = {{cont_var_sym}},
               colour = {{group_label_sym}},
               fill = {{group_label_sym}})) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "NONE") +
    geom_boxplot(width = .2, fill = "white",
                 outlier.shape = NA) +
    ggdist::stat_halfeye(adjust = .33, ## bandwidth
                         width = .67,
                         color = NA,
                         position = position_nudge(x = .15)) +
    gghalves::geom_half_point(side = "l", 
                              range_scale = .3, 
                              alpha = .07, size = 0.1) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, {{group_label_sym}})) +
    coord_flip()
}

individuals_points_plot <- function(df, x_var, y_var, strat_var, shape_var, group_label, colours_order_df) {
  x_var_sym <- rlang::sym(x_var)
  y_var_sym <- rlang::sym(y_var)
  strat_var_sym <- rlang::sym(strat_var)
  shape_var_sym <- rlang::sym(shape_var)
  group_label_sym <- rlang::sym(group_label)
  df %>%
    ggplot(aes(x = {{x_var_sym}},
               y = {{y_var_sym}},
               colour = {{group_label_sym}},
               shape = {{shape_var_sym}},
               group = {{shape_var_sym}})) +
    geom_point(size = 1) +
    guides(shape = guide_legend(override.aes = list(size = 2.5, alpha = 1)),
           colour = "none") +
    theme_minimal() +
    facet_wrap(as.formula(paste0("~", strat_var)), ncol = 1) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10),
          legend.position = "top",
          legend.title = element_text(size = 10)) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}}))
}


linear_points_plot <- function(model1, model2, group_label,
                               colours_order_df, y_pos = 1.5,
                               x_pos = 0.5, add_text = "yes",
                               p_label = "p.adjust") {
  
  df <- model1 %>%
    right_join(model2, by = "participant_genetic_category") %>%
    rename(model1_estimate = estimate.x,
           model1_upper = upper.x,
           model1_lower = lower.x,
           model2_estimate = estimate.y,
           model2_upper = upper.y,
           model2_lower = lower.y,
           response = response.y)

  if (nrow(df) > 1) {
    corr_test <- cor.test(df$model1_estimate, df$model2_estimate)
    corr_p <- corr_test$p.value
    print(corr_p)
    corr_r <- corr_test[[4]]
    corr_test2 <- cor.test(df$model1_estimate, df$model2_estimate, method = "spearman")
    corr_p2 <- corr_test2$p.value
    if (corr_p2 == 0) {
      corr_p2 = 2.2e-16
    }
    corr_r2 <- corr_test2[[4]]
  } else {
    corr_p <- 1
  }
  
  group_label_sym <- rlang::sym(group_label)
  p_label_sym <- rlang::sym(p_label)
  sig_pops <- model2 %>% filter({{p_label_sym}} < 0.05) %>%
    pull({{group_label_sym}})
  colours_order_sig <- colours_order_df %>%
    mutate(colour = ifelse({{group_label_sym}} %in% sig_pops,
                           colour, "white"))

  gg <- df %>%
    ggplot(aes(x = model1_estimate, y = model2_estimate, colour = {{group_label_sym}},
              group = response, fill = {{group_label_sym}})) +
    geom_errorbar(aes(ymax = model2_upper, ymin = model2_lower), width = 0) +
    geom_errorbar(aes(xmax = model1_upper, xmin = model1_lower), width = 0) +
    geom_point(shape = 21) + theme_minimal() +
    theme(panel.grid.minor = element_blank(), legend.title = element_blank()) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    scale_fill_manual(values = pull(colours_order_sig, colour),
                      limits = pull(colours_order_sig, {{group_label_sym}})) +
    guides(colour = guide_legend(ncol = 2, override.aes = list(linetype = 0))) +
      geom_smooth(method = "lm", linetype = "dotted", colour = "grey50", se = FALSE, linewidth = 0.7)
  
  if (add_text == "yes") {
    gg <- gg + geom_text(aes(x = x_pos, y = y_pos,
                  label = paste0("Pr: ", round(corr_r, 2), ", p: ",
                                 format(signif(corr_p, digits = 3), scientific = TRUE),
                                 "\n", "Sr: ", round(corr_r2, 2), ", p: ",
                                 format(signif(corr_p2, digits = 3), scientific = TRUE))),
                  colour = "grey50", size = 2.75)
  }
  return(gg)
}

group_bar_plot <- function(df, y_axis, group_label, colours_order_df) {
  y_axis_sym <- rlang::sym(y_axis)
  group_label_sym <- rlang::sym(group_label)
  df %<>%
    mutate({{group_label_sym}} := factor({{group_label_sym}},
                                          levels = colours_order_df %>%
                                                    pull({{group_label_sym}})))
  df %>%
    ggplot(aes(x = {{group_label_sym}}, y = {{y_axis_sym}},
               fill = {{group_label_sym}})) + 
    geom_bar(stat = "identity", colour = "white") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_manual(values = pull(colours_order_df,
                                    colour),
                      limits = pull(colours_order_df,
                                    {{group_label_sym}})) +
    scale_colour_manual(values = pull(colours_order_df,
                                    colour),
                      limits = pull(colours_order_df,
                                    {{group_label_sym}})) +
    coord_cartesian(clip = "off")
}

pairwise_rccr <- function(df, pop1 = "Africa South") {
  
  relate_coal_shared <- relate_coal_all %>%
    filter(Assignment == participant_genetic_category_two)
  relate_coal_pair <- relate_coal_all %>%
    filter(!groups %in% relate_coal_shared$groups & grepl(pop1, groups)) %>%
    mutate(participant_genetic_category_two = ifelse(Assignment != pop1, Assignment,
                                                     participant_genetic_category_two),
           Assignment = pop1)
  coal_pop1_vec <- relate_coal_shared %>% filter(Assignment == pop1) %>% drop_na() %>%
    pull(coal)
  time <- relate_coal_shared %>% filter(Assignment == pop1) %>% drop_na() %>%
    pull(time)

  other_pops <- relate_coal_pair %>% pull(participant_genetic_category_two) %>%unique()
  rccr_frame_list <- list()
  for (index in 1:length(other_pops)) {
    pop2 <- other_pops[index]
    rccr_frame <- data.frame(time = time, rccr = 0, participant_genetic_category = pop2)
    coal_pop2_vec <- relate_coal_shared %>%
      filter(participant_genetic_category_two == pop2) %>%
      drop_na() %>%
      pull(coal)
    coal_pair_vec <- relate_coal_pair %>%
      filter(participant_genetic_category_two == pop2) %>%
      drop_na() %>%
      pull(coal)
    rccr_vec <- coal_pair_vec / ((coal_pop1_vec + coal_pop2_vec) / 2)
    rccr_frame$rccr <- rccr_vec 
    rccr_frame_list[[index]] <- rccr_frame
  }
  rccr_frames <- do.call(rbind, rccr_frame_list)
  return(rccr_frames)
}

pred_errorbar_plot <- function(df, estimator, algos = c("AlphaMissense", "PolyPhen-2"),
                               tier = "TIERED", tier_to = "Ultra-rare cPAVs") {
  #### Only some have intermediate
  if ("AlphaMissense" %in% algos | "PolyPhen-2" %in% algos) {
    colour_palette <- c("#4CAF50", "orange", "red")
  } else {
    colour_palette <- c("#4CAF50", "red")
  }
  df %>%
    filter(term == estimator) %>%
    mutate(Method = case_when(grepl("AM", response) ~ "AlphaMissense",
                              grepl("Polyphen", response) ~ "PolyPhen-2",
                              grepl("SIFT", response) ~ "SIFT",
                              grepl("PrimateAI", response) ~ "PrimateAI"),
        tiered = ifelse(grepl(tier, response), tier_to, "other"),
        Prediction = case_when(grepl("Tolerated", response) ~ "Tolerated",
                               grepl("Deleterious", response) ~ "Deleterious",
                               grepl("Ambiguous", response) ~ "Ambiguous"),
        Prediction = factor(Prediction, levels = c("Tolerated",
                                                   "Ambiguous",
                                                   "Deleterious")),
        Method = factor(Method, levels = c("SIFT",
                                           "PrimateAI",
                                           "PolyPhen-2",
                                           "AlphaMissense")),
        sig_stars = case_when(
          p.value < 0.1  & p.value >= 0.05  ~ ".",
          p.value < 0.05 & p.value >= 0.01 ~ "*",
          p.value < 0.01 & p.value >= 0.001 ~ "**",
          p.value < 0.001 ~ "***",
          is.na(p.value) ~ NA_character_,
          TRUE ~ NA_character_)) %>%
    filter(Method %in% algos & tiered == tier_to) %>%
    ggplot(aes(y = Method, x = estimate, colour = Prediction)) +
      geom_errorbar(aes(xmax = upper, xmin = lower), width = 0.1,
                        position = position_dodge(width = 0.8)) +
      theme_minimal() +
      facet_wrap(~tiered) +
      geom_text(aes(label = sig_stars), vjust = -0.045, position = position_dodge(width = 0.8)) +
      theme(panel.grid.major.y = element_blank(),
            panel.grid.minor = element_blank(),
            legend.position = "top",
            strip.text = element_text(size = 10)) +
      geom_point(position = position_dodge(width = 0.8)) +
      scale_colour_manual(values = colour_palette) +
      guides(colour = guide_legend(override.aes = list(linetype = 0, label = 0, size = 2.5)))
}

pred_miss_points_plot <- function(df, algo = "AlphaMissense", tier = "TIERED",
                                  zoom_del = TRUE, y_zoom = 70, xvar = "missing_variants_count_tho") {
  xvar_sym <- rlang::sym(xvar)
  #### Only some have intermediate
  if (algo %in% c("SIFT", "PrimateAI")) {
    colour_palette <- c("#4CAF50", "red")
  } else {
    colour_palette <- c("#4CAF50", "orange", "red")
  }
  df %<>%
    dplyr::select(c("participant_id", contains(c("missing", "_SIFT_", "_AM_",
                                                 "_Polyphen2_HVAR_", "_PrimateAI_")))) %>%
    pivot_longer(contains(c("SIFT_", "AM_", "Polyphen2_HVAR_", "PrimateAI_")),
                 names_to = "response", values_to = "count") %>%
    mutate(Method = case_when(grepl("AM", response) ~ "AlphaMissense",
                              grepl("Polyphen", response) ~ "PolyPhen-2",
                              grepl("SIFT", response) ~ "SIFT",
                              grepl("PrimateAI", response) ~ "PrimateAI"),
           tiered = case_when(grepl("TIERED_", response) ~ "TIERED",
                             grepl("TIER12_", response) ~ "TIER12",
                             TRUE ~ NA_character_),
           Prediction = case_when(grepl("Tolerated", response) ~ "Tolerated",
                                 grepl("Deleterious", response) ~ "Deleterious",
                                 grepl("Ambiguous", response) ~ "Ambiguous"),
           Prediction = factor(Prediction, levels = c("Tolerated",
                                                     "Ambiguous",
                                                     "Deleterious"))) %>%
    drop_na() %>%
    filter(tiered == tier)
    df2 <- df %>%
      filter(Method == algo)

    gg <- df2 %>%
      ggplot(aes(x = {{xvar_sym}}, y = count,
                 colour = Prediction, fill = Prediction)) +
      geom_point(alpha = 0.05, size = 0.25) +
      geom_smooth(method = "lm", linewidth = 0.8) +
      theme_minimal() +
      facet_wrap(~Method) +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            strip.text = element_text(size = 10)) +
      scale_colour_manual(values = colour_palette) +
      scale_fill_manual(values = colour_palette) +
      theme(legend.position = "NONE") +
      coord_cartesian(ylim = c(0, max(df$count) + max(df$count)/5))
    if (zoom_del) {
      
      gg_zoom <- df2 %>%
        filter(Prediction != "Tolerated") %>%
        ggplot(aes(x = {{xvar_sym}}, y = count,
               colour = Prediction, fill = Prediction)) +
          geom_point(alpha = 0.05, size = 0.25) +
          geom_smooth(method = "lm", linewidth = 0.8) +
          theme_minimal() +
          coord_cartesian(ylim = c(0, y_zoom)) +
          theme(panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                legend.position = "NONE",
                axis.title = element_blank()) +
          scale_colour_manual(values = colour_palette[-1]) +
          scale_fill_manual(values = colour_palette[-1])

      gg <- gg +
        annotation_custom(ggplotGrob(gg_zoom), 
                           xmin = min(df %>% pull({{xvar_sym}})) + 2.75 ,
                           xmax = max(df %>% pull({{xvar_sym}})) / 1.8,
                           ymin = max(df$count) - max(df$count)/3,
                           ymax = max(df$count) + max(df$count)/9) +
        geom_rect(aes(xmin = min(df %>% pull({{xvar_sym}})) + 2.75 ,
                      xmax = max(df$missing_variants_count_tho) / 1.8,
                  ymin = max(df$count) - max(df$count)/3,
                  ymax = max(df$count) + max(df$count)/9),
            colour = "grey20", linetype='dotted', alpha=0, size = 0.5) +
        geom_rect(aes(xmin = min(df %>% pull({{xvar_sym}})) - 0.25,
                      xmax = max(df %>% pull({{xvar_sym}})) + 0.25,
                      ymin = max(df$count) / 100, ymax = y_zoom),
                  color = "grey20", alpha = 0, size = 0.5, linetype = "dotted")
    }
  return(gg)
}

tier_pred_prop_table <- function(df, tier, algo,
                                 encoding, group_label) {
  group_label_sym <- rlang::sym(group_label)
  tier_all <- rlang::sym(paste0(tier, "_", algo))
  tier_effect <- rlang::sym(paste0(tier, "_", algo, "_", encoding))
  df %>%
    group_by({{group_label_sym}}) %>%
    summarise({{tier_all}} := sum({{tier_all}}),
              {{tier_effect}} := sum({{tier_effect}}),
              tiered_prop = ({{tier_effect}} / {{tier_all}}) * 100) %>%
    dplyr::select({{group_label_sym}}, tiered_prop) %>%
    mutate(method = algo)
}

exit_pred_var_boxplot <- function(df, exit_annot, acmg_class = "vus",
                                  xvar = "missing_variants_count_tho",
                                  algos = c("AlphaMissense", "PrimateAI"),
                                  tiers = c("TIER1", "TIER2"),
                                  group_label, colours_order_df) {

  group_label_sym <- rlang::sym(group_label)
  xvar_sym <- rlang::sym(xvar)
  tiers_grep <- paste(tiers, collapse = "|")
  algo_acmg_frame <- exit_annot %>%
    filter(participant_id %in% df$participant_id,
           acmg_classification_broad == acmg_class,
           grepl(tiers_grep, tier)) %>%
    pivot_longer(!c(varkey, participant_id, tier, acmg_classification_broad),
                 names_to = "algo", values_to = "prediction") %>%
    mutate(algo = case_when(algo == "AM_pred" ~ "AlphaMissense",
                            algo == "PrimateAI_pred" ~ "PrimateAI",
                            algo == "SIFT_pred" ~ "SIFT",
                            algo == "Polyphen2_pred" ~ "PolyPhen-2")) %>% 
    filter(algo %in% algos) %>%
    drop_na(prediction) %>%
    left_join(df %>% dplyr::select(participant_id,
                                   {{xvar_sym}},
                                   {{group_label_sym}}))
  cat(paste0("Number of VUS: ", nrow(algo_acmg_frame), "\n"))
  if ("AlphaMissense" %in% algos | "PolyPhen-2" %in% algos) {
    algo_acmg_frame %<>%
      mutate(prediction := factor(prediction, levels = c("Tolerated", "Deleterious", "Ambiguous")))
  } else {
    algo_acmg_frame %<>%
      mutate(prediction := factor(prediction, levels = c("Tolerated", "Deleterious")))
  }

  algo_acmg_frame %>%
    ggplot(aes(x = prediction, y = {{xvar_sym}},
               colour = {{group_label_sym}}, group = prediction)) +
    geom_jitter(size = 1) +
    geom_boxplot(colour = "grey50") +
    theme_minimal() +
    theme(legend.position = "NONE",
          strip.background = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.text = element_text(size = 11)) +
    facet_grid(~algo, scales = "free", space = "free") +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    stat_pvalue_manual(algo_acmg_frame %>%
                         group_by(algo) %>%
                         pairwise_wilcox_test(as.formula(paste0({{xvar_sym}}, "~", "prediction"))) %>%
                         add_sig_stars(label = "p.adj.signif", p_label = "p.adj") %>%
                         add_xy_position(scales = "free_y"))
}

add_sig_stars <- function(df, label = "sig_stats", p_label) {
    p_label_sym <- rlang::sym(p_label)
    label_sym <- rlang::sym(label)
    ### Add stuff
    df %>%
    mutate({{label_sym}} := case_when(
        {{p_label_sym}} < 0.1  & {{p_label_sym}} >= 0.05  ~ ".",
        {{p_label_sym}} < 0.05 & {{p_label_sym}} >= 0.01 ~ "*",
        {{p_label_sym}} < 0.01 & {{p_label_sym}} >= 0.001 ~ "**",
        {{p_label_sym}} < 0.001 ~ "***",
        is.na({{p_label_sym}}) ~ NA_character_,
        TRUE ~ NA_character_))
}

## For making the figure...
# Given a vector of true AFs (aka. probability of observing allele on a single haplotype, assuming HWE)..
# .. what is the probability that a variant at that true AF will appear above some AC threshold (given an AF threshold) in a sample composed of a certain number of haplotypes (pop_size) drawn from that population.
filtering_confidence <- function(pop_size = 1000, threshold = 0.01, max_ac = "AF", step = 0.001,
                                 start_freq = 0.001, end_freq = 0.1) {
  if (max_ac == "AF95") {
    max_ac <- find_max_ac(threshold, pop_size)
  } else if (max_ac == "AF") {
    max_ac <- floor(pop_size * threshold)
  } else {
    max_ac <- max_ac
  }
  freq_vec <- seq(start_freq, end_freq, by = step)
  filt_likelihood <- pbinom(max_ac, pop_size, freq_vec, lower.tail = FALSE)
  return(data.frame(
    true_freq = freq_vec,
    filt_freq = threshold,
    filt_like = filt_likelihood
  ))
}

find_max_ac = function(af,an,ci=.95) {
  if (af == 0) {
    return (0)
  } else {
    quantile_limit = ci # ci for one-sided, 1-(1-ci)/2 for two-sided
    max_ac = qpois(quantile_limit,an*af)
    return (max_ac)
  }
}

filtering_confidence <- function(pop_size = 1000, threshold = 0.01, max_ac = "AF", step = 0.001,
                                 start_freq = 0.001, end_freq = 0.1) {
  if (max_ac == "FAF95") {
    max_ac <- find_max_ac(threshold, pop_size)
  } else if (max_ac == "AF") {
    max_ac <- floor(pop_size * threshold)
  } else {
    max_ac <- max_ac
  }
  freq_vec <- seq(start_freq, end_freq, by = step)
  filt_likelihood <- pbinom(max_ac, pop_size, freq_vec, lower.tail = FALSE)
  return(data.frame(
    true_freq = freq_vec,
    filt_freq = threshold,
    filt_like = filt_likelihood
  ))
}

pca_plot <- function(
    x, axis1, axis2, group_label, colours_order_df
    ) {
    group_label_sym <- rlang::sym(group_label)
    axis1_sym <- rlang::sym(axis1)
    axis2_sym <- rlang::sym(axis2)
    ggplot(x) +
        aes(x = {{axis1_sym}}, y = {{axis2_sym}},
            colour = {{group_label_sym}}) +
        geom_point(size = 0.5) +
        theme_classic() +
        xlab(axis1) + ylab(axis2) +
        scale_colour_manual(values = pull(colours_order_df, colour),
                            limits = pull(colours_order_df, {{group_label_sym}}))
}

zoomed_bar_plot <- function(df, group_label,
                           colours_order_df,
                           zoom = c(22300, 29100)) {

  group_label_sym <- rlang::sym(group_label)
  df_counts <- df %>%
    group_by({{group_label_sym}}) %>%
    tally() %>%
    arrange(n) %>%
    mutate({{group_label_sym}} := factor({{group_label_sym}},
                                         levels = unique(pull(., {{group_label_sym}})))) 

  df_counts %>%
    ggplot(aes(y = "", x = n,
               fill = {{group_label_sym}})) +
    geom_bar(stat = "identity", colour = "white") +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid = element_blank()) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, {{group_label_sym}})) +
    ylab("") +
    facet_zoom(xlim = zoom)
}

errorbar_plot2 <- function(df, group_label, reference_group = NULL,
                           colours_order_df, p_label = "p.value") {
    group_label_sym <- rlang::sym(group_label)
    p_label_sym <- rlang::sym(p_label)
    ### Add stuff
    df %<>%
    mutate(sig_stars = case_when(
        {{p_label_sym}} < 0.1  & {{p_label_sym}} >= 0.05  ~ ".",
        {{p_label_sym}} < 0.05 & {{p_label_sym}} >= 0.01 ~ "*",
        {{p_label_sym}} < 0.01 & {{p_label_sym}} >= 0.001 ~ "**",
        {{p_label_sym}} < 0.001 ~ "***",
        is.na({{p_label_sym}}) ~ NA_character_,
        TRUE ~ NA_character_),
        {{group_label_sym}} := factor({{group_label_sym}},
            levels = pull(colours_order_df, {{group_label_sym}})))
    
    df %>% filter({{p_label_sym}} < 0.05) %>%
      pull({{group_label_sym}}) -> sig_pops
    colours_order_sig <- colours_order_df %>%
      mutate(colour = ifelse({{group_label_sym}} %in% sig_pops,
                             colour, "white"))

    ## Adjust top group if requested
    if(!is.null(reference_group)) {
      df %<>%
        mutate({{group_label_sym}} := relevel({{group_label_sym}},
              ref = reference_group))
    }
    ### Plot
    ggplot(df) +
    aes(y = estimate, x = {{group_label_sym}},
        colour = {{group_label_sym}},
        fill = {{group_label_sym}}) +
    geom_line(yintercept = 0, linewidth = 0.7, linetype = "dashed", colour = "grey50") +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
    geom_point(shape = 21, size = 2.5) +
    geom_text(data = df,
        aes(y = estimate, x = {{group_label_sym}},
            label = sig_stars), hjust = -0.05) +
    theme_minimal() +
    coord_cartesian(clip = 'off') +
    ylab("") +
    theme(legend.position = "NONE",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    scale_fill_manual(values = pull(colours_order_sig, colour),
                      limits = pull(colours_order_sig, {{group_label_sym}}))
}


raincloud_plot2 <- function(df, cont_var, group_label, colours_order_df) {
  cont_var_sym <- rlang::sym(cont_var)
  group_label_sym <- rlang::sym(group_label)
  df %<>%
    mutate({{group_label_sym}} := factor({{group_label_sym}},
                                          levels = rev(colours_order_df %>%
                                                    pull({{group_label_sym}}))))
  df %>%
    ggplot(aes(x = {{group_label_sym}},
               y = {{cont_var_sym}},
               colour = {{group_label_sym}},
               fill = {{group_label_sym}})) +
    theme_minimal() +
    theme(panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "NONE") +
    geom_boxplot(width = .2, fill = "white",
                 outlier.shape = NA) +
    ggdist::stat_halfeye(adjust = .75, ## bandwidth
                         width = .67,
                         color = NA,
                         position = position_nudge(x = .15)) +
    gghalves::geom_half_point(side = "l", 
                              range_scale = .3, 
                              alpha = .07, size = 0.1) +
    scale_colour_manual(values = pull(colours_order_df, colour),
                        limits = pull(colours_order_df, {{group_label_sym}})) +
    scale_fill_manual(values = pull(colours_order_df, colour),
                      limits = pull(colours_order_df, {{group_label_sym}})) +
    coord_flip()
}

filtering_conf_af_faf95_plot <- function(covid_sizes, step = 0.0001, threshold = 0.01,
                                         start_freq = 0.001, end_freq = 0.1, colours_order_df) {

  ## Loop through pops and create simulations
  pop_list_filter <- list()
  covid_pops <- covid_sizes$participant_genetic_category
  for (index in 1:length(covid_pops)) {
    pop <- covid_pops[index]
    haplotypes <- as.numeric(covid_sizes[index, 3])

    sim_adjust <- filtering_confidence(pop_size = haplotypes, step = step, threshold = threshold,
                                       start_freq = start_freq, end_freq = end_freq,
                                       max_ac = "FAF95") %>%
      mutate(participant_genetic_category = pop, max_ac = "FAF95")

    sim_noadjust <- filtering_confidence(pop_size = haplotypes, step = step, threshold = threshold,
                                         start_freq = start_freq, end_freq = end_freq) %>%
      mutate(participant_genetic_category = pop, max_ac = "AF")

    pop_list_filter[[index]] <- rbind(sim_adjust, sim_noadjust)
  }

  ## Combine simulations
  sims_all <- do.call(rbind, pop_list_filter) %>%
    mutate(participant_genetic_category = factor(participant_genetic_category,
                                                 levels = colours_order_df$participant_genetic_category))
  
  ## Plot
  ggplot(sims_all) +
    aes(x = true_freq, y = filt_like, group = participant_genetic_category,
        colour = participant_genetic_category) +
    geom_line(linewidth = 0.75) +
    facet_wrap(~max_ac, ncol = 1) +
    theme_minimal() +
    annotation_logticks(base = 10, side = "b") +
    scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "NONE",
          strip.text = element_text(size = 10)) +
    scale_colour_manual(values = pull(colours_order_df %>% filter(participant_genetic_category %in%
                                                                  covid_sizes$participant_genetic_category),
                                      colour),
                        limits = pull(colours_order_df %>% filter(participant_genetic_category %in%
                                                                  covid_sizes$participant_genetic_category),
                                      participant_genetic_category))
}

faf_filter_barplot <- function(covid19_faf95_filter_tiering, filter_test = "rare_cPAV_common_FAF95",
                               colours_order_df) {
  # Create plot                            
  covid19_faf95_filter_tiering %>%
    filter(test == filter_test) %>%
    mutate(participant_genetic_category = factor(participant_genetic_category,
           levels = colours_order_df$participant_genetic_category)) %>%
    ggplot(aes(x = participant_genetic_category, y = num_instances,
               colour = participant_genetic_category,
               fill = participant_genetic_category,
               group = participant_genetic_category)) +
    geom_bar(stat = "identity", colour = "white", width = 0.65) + theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "NONE") +
    scale_colour_manual(values = pull(colours_order_df %>% filter(participant_genetic_category %in% covid19_faf95_filter_tiering$participant_genetic_category),
                        colour),
                        limits = pull(colours_order_df %>% filter(participant_genetic_category %in% covid19_faf95_filter_tiering$participant_genetic_category),
                                      participant_genetic_category)) +
    scale_fill_manual(values = pull(colours_order_df %>% filter(participant_genetic_category %in% covid19_faf95_filter_tiering$participant_genetic_category),
                      colour),
                      limits = pull(colours_order_df %>% filter(participant_genetic_category %in% covid19_faf95_filter_tiering$participant_genetic_category),
                                    participant_genetic_category)) +
    geom_hline(yintercept = 0, linewidth = 0.5, colour = "grey90") +
    coord_cartesian(clip = "off") +
    geom_text(aes(x = participant_genetic_category, y = num_instances,
                label = num_instances), vjust = -0.2, size = 3.5)
}

###########################################################################
###########################################################################
###                                                                     ###
###                          DATA INGEST                                ###
###                                                                     ###
###########################################################################
###########################################################################

### Next, we will read in the relevant aesthetics data from /data/aesthetics/
### This provides information on the colours + order of the cohort labels utilised throughout the paper.

##### Ethnicity Aesthetics Data
ethn_colours <- read_csv("aesthetics_data/ethnicity_colours.csv") %>%
  configure_ethnicity_labels()

##### Ancestry Aesthetics Data
prive_colours <- read_csv("aesthetics_data/prive_colours.csv")

### gnomAD order
gnomad_order <- read_csv("aesthetics_data/gnomad_order.csv")

### PCA UMAP
pca_umap <- read_csv(paste0(analysis_dir, "pca_analysis/pca_umap.csv"))

### Ancestry coeffiencts data
ukbb_assignments_data <- read_csv(paste0(out_dir, "ancestry_data/ukbb_assignments.csv")) %>%
  dplyr::select(-Assignment)

##### GLM analysis data (unrelated etc.)
analysis_data <- read_csv(paste0(analysis_dir, "glm_analysis/analysis_data.csv.gz"))
analysis_data_plots <- analysis_data %>%
  configure_plotting_variables(colours_order_df = prive_colours) %>%
  left_join(pca_umap, by = "platekey") %>%
  left_join(ukbb_assignments_data, by = "platekey")

### cPAVs counts
cpavs_counts <- read_csv(paste0(analysis_dir, "glm_analysis/cpavs_counts_rd.csv.gz"))

#### GLM tables for plotting
ukbb_cpav_glm <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_cpavs_glm.csv"))
ukbb_pavs_miss_glm <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_missing_glm.csv"))
ukbb_patho_pred_glm <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_patho_prop_glm.csv"))

### Needs some editing for plotting niceties
ukbb_path_count <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_patho_count_glm.csv")) %>%
  mutate(method = case_when(grepl("AM", response) ~ "AlphaMissense > 0.564",
                            grepl("PAI3D", response) ~ "PrimateAI-3D > 0.8",
                            grepl("CADD", response) ~ "CADD > 30"),
         method = factor(method, levels = c("CADD > 30", "PrimateAI-3D > 0.8", "AlphaMissense > 0.564")))

ukbb_vus_glm <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_vus_glm.csv")) %>%
  mutate(variable = ifelse(variable == "Yes", "Diagnosed cohort", "Undiagnosed cohort"))
####

ukbb_ppv_glm <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_ppv_glm.csv"))
ukbb_diagnostic_yield_glm <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_diagnostic_yield_glm.csv"))
gnomad_diagnostic_yield_glm <- read_csv(paste0(analysis_dir, "glm_analysis/gnomad_diagnostic_yield_glm.csv"))

### ANOVA for Yield GLM
lrt_tiered_yield <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_diagnostic_yield_anova.csv"))
## Full model GLM diagnostic yield from GPcPAVs
model_tidy <- read_csv(paste0(analysis_dir, "glm_analysis/ukbb_diagnostic_yield_full_model_glm.csv"))

### COVID-19 AF filter stuff
covid19_faf95_filter_tiering <- read_csv(paste0(analysis_dir, "covid_filter_analysis/aggCOVID_faf95_filter_tiering.csv"))

### Age, sex ethnicity matrix.
ons_rd_age_sex_ethn_matrix <- fread(paste0(analysis_dir, "ons_ethnicity_analysis/ons_100kGP_probands_age_sex_ethnicity.csv.gz"))

### Missing PAVs count (from the data directory...)
missing_gnomad_pavs_data <- read_csv(paste0(out_dir, "combined_data/combined_missing_pavs_gnomad.csv.gz")) %>%
  mutate(missing_pavs_count = HOM_ALT_CT + HET_CT)

###########################################################################
###########################################################################
###                                                                     ###
###                        PLOTTING (MAIN FIGURES)                      ###
###                                                                     ###
###########################################################################
###########################################################################

## Plots create (for outputting plots)
dir.create(plots_dir)

####################
#---- FIGURE 1 ----#
####################

figure_1a_gg <- analysis_data_plots %>%
  filter(participant_type == "Proband") %>%
  zoomed_bar_plot(group_label = "participant_genetic_category",
                  colours_order_df = prive_colours %>%
                   filter(participant_genetic_category %in% analysis_data_plots$participant_genetic_category)) +
  xlab("100kGP Probands") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(title="GIA group"))

figure_1b_gg <- analysis_data_plots %>%
  filter(participant_type == "Proband") %>%
  pca_plot(axis1 = "PC1", axis2 = "PC2",
           group_label = "participant_genetic_category",
           colours_order_df = prive_colours %>%
            filter(participant_genetic_category %in% analysis_data_plots$participant_genetic_category)) +
  theme(legend.position = "NONE")

figure_1_gg <- plot_grid(figure_1a_gg, figure_1b_gg, ncol = 1, rel_heights = c(0.45, 1))
ggsave(paste0(plots_dir, "figure1.pdf"), width = 8.5, height = 10)

####################
#---- FIGURE 2 ----#
####################

figure_2_gg <- errorbar_plot(ukbb_cpav_glm %>%
                filter(response == "cPAV"),
              "participant_genetic_category", NULL,
              prive_colours, "p.value") +
  ggtitle("cPAVs, n: 1,951,659") +
  xlab("log(RR) vs Europe North West") +
  coord_cartesian(xlim = c(-0.55, 1.25)) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("GIA group") +
  theme_bw() + 
  theme(legend.position = "NONE",
        plot.title = element_text(size = 11),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

#### Plot out (PDF)
ggsave(paste0(plots_dir, "figure2.pdf"), width = 5, height = 5)

####################
#---- FIGURE 3 ----#
####################

figure_3a_gg <- analysis_data_plots %>%
  raincloud_plot2("missing_pavs_count",
                 "participant_genetic_category",
                 prive_colours) +
  xlab("GIA group") + ylab("Protein-altering variants\nmissing from gnomAD\n") +
  scale_y_continuous(limits = c(0, 250),
                     breaks = c(0, 50, 100, 150, 200, 250))

figure_3b_gg <-
  linear_points_plot(model1 = ukbb_pavs_miss_glm,
                     model2 = ukbb_cpav_glm %>%
                       filter(response == "cPAV"),
                     group_label = "participant_genetic_category",
                     colours_order_df = prive_colours,
                     y_pos = 1.1, x_pos = 0.35) +
  ylab("log(RR) cPAVs\nvs Europe North West") +
  xlab("log(RR) Protein-altering variants\nmissing from gnomAD\nvs Europe North West") +
  theme(legend.position = "NONE") +
  geom_abline(slope = 1, size = 0.5, colour = "grey92") +
  coord_obs_pred()

figure_3c_gg <- analysis_data_plots %>%
  filter(family_group_type_broad != "Other" & penetrance %in% c("complete", "incomplete", "n/a") & cPAV > 0) %>%
  individuals_points_plot("missing_pavs_count",
                          "cPAV",
                          "family_group_type_broad_long",
                          "Penetrance",
                          "participant_genetic_category",
                          prive_colours) +
  xlab("Protein-altering variants\nmissing from gnomAD\n") +
  ylab("\ncPAVs") +
  scale_x_continuous(limits = c(0, 250),
                     breaks = c(0, 50, 100, 150, 200, 250)) +
  scale_y_continuous(limits = c(min(analysis_data_plots$cPAV),
                                max(analysis_data_plots$cPAV)))

figure3_gg <- plot_grid(plot_grid(figure_3a_gg, figure_3b_gg,
                                  nrow = 2, rel_heights = c(2.65, 1)),
                                  figure_3c_gg, ncol = 2, rel_widths = c(1, 1.35))
ggsave(paste0(plots_dir, "figure3.pdf"), width = 8.5, height = 9.5)

####################
#---- FIGURE 4 ----#
####################

figure_4a_gg <- analysis_data_plots %>%
  tier_pred_prop_table(
  "rare_cPAV",
  "CADD",
  "Deleterious",
  "participant_genetic_category"
) %>%
  group_bar_plot("tiered_prop", "participant_genetic_category", prive_colours) +
  ylab("Score >30\nDeleterious (%)") +
  xlab(" ") +
  ggtitle(paste0("CADD\n",
          "1,828,418",
          " cPAVs (popmax AF < 0.1%)")) +
  theme(plot.title = element_text(size = 11)) +
  coord_cartesian(ylim = c(0, 15))

figure_4b_gg <- linear_points_plot(
  model1 = ukbb_cpav_glm %>%
    filter(response == "rare_cPAV"),
  model2 = ukbb_patho_pred_glm %>%
    filter(response == "rare_cPAV_CADD_Deleterious"),
  group_label = "participant_genetic_category",
  colours_order_df = prive_colours,
  y_pos = -0.45,
  x_pos = 0.05,
  p_label = "p.value"
) +
  ylab("\nlog(OR) Deleterious (%)\nvs Europe North West") +
  xlab("log(RR) cPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") + ggtitle(" ")

figure_4ab_gg <- plot_grid(figure_4a_gg, figure_4b_gg, ncol = 2, rel_widths = c(1.8, 1))

figure_4c_gg <- analysis_data_plots %>%
  tier_pred_prop_table(
  "rare_cPAV",
  "PAI3D",
  "Deleterious",
  "participant_genetic_category"
) %>%
  group_bar_plot("tiered_prop", "participant_genetic_category", prive_colours) +
  ylab("Score > 0.8\nDeleterious (%)") +
  xlab(" ") +
  ggtitle(paste0("PrimateAI-3D\n",
          "1,340,758",
          " missense cPAVs (popmax AF < 0.1%)")) +
  theme(plot.title = element_text(size = 11))  +
  coord_cartesian(ylim = c(0, 15))

figure_4d_gg <- linear_points_plot(
  model1 = ukbb_cpav_glm %>%
    filter(response == "rare_cPAV"),
  model2 = ukbb_patho_pred_glm %>%
    filter(response == "rare_cPAV_PAI3D_Deleterious"),
  group_label = "participant_genetic_category",
  colours_order_df = prive_colours,
  y_pos = -0.45,
  x_pos = 0.05,
  p_label = "p.value"
) +
  ylab("\nlog(OR) Deleterious (%)\nvs Europe North West") +
  xlab("log(RR) cPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") + ggtitle(" ")

figure_4cd_gg <- plot_grid(figure_4c_gg, figure_4d_gg, ncol = 2, rel_widths = c(1.8, 1))

figure_4e_gg <- analysis_data_plots %>%
  tier_pred_prop_table(
  "rare_cPAV",
  "AM",
  "Deleterious",
  "participant_genetic_category"
) %>%
  group_bar_plot("tiered_prop", "participant_genetic_category", prive_colours) +
  ylab("Score > 0.564\nDeleterious (%)") +
  xlab(" ") +
  ggtitle(paste0("AlphaMissense\n",
          "1,369,524",
          " missense cPAVs (popmax AF < 0.1%)")) +
  theme(plot.title = element_text(size = 11))  +
  coord_cartesian(ylim = c(0, 15))

figure_4f_gg <- linear_points_plot(
  model1 = ukbb_cpav_glm %>%
    filter(response == "rare_cPAV"),
  model2 = ukbb_patho_pred_glm %>%
    filter(response == "rare_cPAV_AM_Deleterious"),
  group_label = "participant_genetic_category",
  colours_order_df = prive_colours,
  y_pos = -0.45,
  x_pos = 0.05,
  p_label = "p.value"
) +
  ylab("\nlog(OR) Deleterious (%)\nvs Europe North West") +
  xlab("log(RR) cPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") +
  ggtitle(" ")

figure_4ef_gg <- plot_grid(figure_4e_gg, figure_4f_gg, ncol = 2, rel_widths = c(1.8, 1))

### Full figure
figure_4_gg <- plot_grid(figure_4ab_gg, figure_4cd_gg, figure_4ef_gg, ncol = 1)
ggsave(paste0(plots_dir, "figure4.pdf"), width = 8.75, height = 9.75)

####################
#---- FIGURE 5 ----#
####################

### COVID-19 sizes circle plot (needs editing in post...)
covid_sizes <- missing_gnomad_pavs_data %>%
  filter(dataset == "COVID-19 Cohort") %>%
  group_by(participant_genetic_category) %>%
  tally() %>%
  filter(n > 100 & participant_genetic_category != "Remaining participants") %>%
  mutate(haplotypes = n * 2)

### Circle packing
packing <- circleProgressiveLayout(covid_sizes$n, sizetype='area')
covid_sizes_packed <- cbind(covid_sizes, packing) %>%
  mutate(id = as.integer(rownames(.)), participant_genetic_category = factor(participant_genetic_category,
                                                                             levels = prive_colours$participant_genetic_category)) %>%
  dplyr::select(id, participant_genetic_category)
covid_sizes_packed.gg <- circleLayoutVertices(packing, npoints = 100) %>%
  left_join(covid_sizes_packed, by = "id")

### Figure 5a (also used for) Supp Figure 19a
figure5a_gg <- ggplot() +
  geom_polygon(data = covid_sizes_packed.gg, aes(x, y, group = as.factor(id),
               fill = participant_genetic_category), colour = "white") +
  scale_size_continuous(range = c(1,5)) + theme_void() + coord_equal() +
  theme(legend.position = "bottom", legend.title = element_text()) +
  guides(fill = guide_legend(ncol = 3, title = "1000 participants")) +
  scale_fill_manual(values = pull(prive_colours %>% filter(participant_genetic_category %in% covid_sizes_packed.gg$participant_genetic_category),
                                    colour),
                    limits = pull(prive_colours %>% filter(participant_genetic_category %in% covid_sizes_packed.gg$participant_genetic_category),
                                  participant_genetic_category)) +
  ggtitle("COVID-19 cohort")
ggsave(paste0(plots_dir, "figure5a.pdf"), height = 5, width = 5)

### Figure 5b
figure5b_gg <- filtering_conf_af_faf95_plot(covid_sizes = covid_sizes,
                                            colours_order_df = prive_colours) +
  ylab("\n\nProbability observed >1%") +
  xlab("True AF")
ggsave(paste0(plots_dir, "figure5b.pdf"), height = 4.5, width = 4.5)

### Figure 5c
figure5c_gg <- faf_filter_barplot(covid19_faf95_filter_tiering = covid19_faf95_filter_tiering,
                                  colours_order_df = prive_colours) +
  scale_y_continuous(breaks = c(0, 500, 1000)) +
  xlab("GIA group\nCOVID-19 cohort") +
  ylab("cPAVs (popmax AF < 0.1%)\nCOVID-19 cohort FAF95 > 1%") +
  ggtitle("")
ggsave(paste0(plots_dir, "figure5c.pdf"), width = 8, height = 4)

####################
#---- FIGURE 6 ----#
####################

figure_6a_gg <- analysis_data %>%
  group_by(participant_genetic_category) %>%
  summarise(GPcPAV_ppv = (sum(GPcPAV_ACMG_Pathogenic) / sum(GPcPAV)) * 100) %>%
  group_bar_plot("GPcPAV_ppv", "participant_genetic_category", prive_colours) +
  ylab("P/LP (%)") +
  ggtitle("Positive predictive value (PPV) of GPcPAVs") +
  xlab("GIA group") +
  theme(plot.title = element_text(size = 11))

figure_6b_gg <-
  linear_points_plot(model1 = ukbb_cpav_glm %>%
                      filter(response == "GPcPAV"),
                     model2 = ukbb_ppv_glm %>%
                       filter(response == "GPcPAV_ACMG_Pathogenic"),
                     group_label = "participant_genetic_category",
                     colours_order_df = prive_colours,
                     y_pos = -1.25, x_pos = 0.15, p_label = "p.adjust") +
  ylab("\nlog(RR) PPV\nvs Europe North West") +
  xlab("log(RR) GPcPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") +
  ggtitle(" ")

figure_6ab_gg <- plot_grid(figure_6a_gg, figure_6b_gg, ncol = 2, rel_widths = c(1.9, 1))

### Stats including gnomAD + prive assignments
gnomad_prive_cont_stats <- analysis_data %>%
mutate(participant_genetic_category = factor(participant_genetic_category,
       levels = prive_colours$participant_genetic_category),
       participant_gnomad_category = factor(participant_gnomad_category,
       levels = rev(c("nfe", "sas", "afr", "eas")))) %>%
filter(participant_gnomad_category %in% c("nfe", "afr", "eas", "sas")) %>%
cat_var_stats("participant_gnomad_category", "participant_genetic_category")

figure6_c_gg <- gnomad_prive_cont_stats %>%
ggplot(aes(y = participant_gnomad_category, x = proportion * 100,
           fill = participant_genetic_category,
           colour = participant_genetic_category)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(legend.position = "none",
          axis.line = element_blank(),
          panel.border = element_rect(linewidth = 0.2, colour = "grey50"),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank()) +
    scale_fill_manual(values = pull(prive_colours,
                                    colour),
                      limits = pull(prive_colours,
                                    participant_genetic_category)) +
     scale_colour_manual(values = pull(prive_colours,
                                    colour),
                      limits = pull(prive_colours,
                                    participant_genetic_category)) +
    coord_cartesian(clip = "off") +
    xlab("GIA group (%)") +
    ylab("Continental GIA group") + ggtitle("\n")

figure6_d_gg <- gnomad_diagnostic_yield_glm %>%
  filter(!participant_gnomad_category %in% c("oth") & response == "DIAGNOSED_GPcPAV_ACMG_Pathogenic") %>%
  mutate(participant_gnomad_category = factor(participant_gnomad_category,
         levels = rev(c("nfe", "sas", "afr", "eas")))) %>%
  ggplot(aes(y = participant_gnomad_category, x = estimate)) +
  geom_point(colour = "grey25") +
  geom_errorbar(aes(xmax = upper, xmin = lower), colour = "grey25", width = 0.1) +
  geom_vline(xintercept = 0, linewidth = 0.7, linetype = "dashed", colour = "grey50") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(linewidth = 0.2, colour = "grey50"),
        panel.grid.major.y = element_blank()) +
  xlab("\nlog(OR) vs nfe") +
  ylab("") + ggtitle("\nCase solved with GPcPAV(s)") +
  theme(plot.title = element_text(size = 11))

figure_6cd_gg <- plot_grid(figure6_c_gg, figure6_d_gg, align = "h", rel_widths = c(1, 1))

### Full figure
figure_6_gg <- plot_grid(figure_6ab_gg, figure_6cd_gg, ncol = 1, rel_heights = c(1, 0.75))
ggsave(paste0(plots_dir, "figure6.pdf"), width = 8.75, height = 5.5)

###########################################################################
###########################################################################
###                                                                     ###
###                  PLOTTING (SUPPLEMENTARY FIGURES)                   ###
###                                                                     ###
###########################################################################
###########################################################################

############################
##### Supplementary Figure 1

supp_figure1_gg <- ons_rd_age_sex_ethn_matrix %>%
  configure_ethnicity_labels() %>%
  age_ethn_trend_plot(group_label = "ethn_label",
                      colours_order_df = ethn_colours)
ggsave(paste0(plots_dir, "supp_figure1.pdf"), width = 6, height = 6)

############################
##### Supplementary Figure 2

### Only get odd numbers as list indices
justodd <- function(x) x[x %% 2 == 1]
pcs_vec <- justodd(1:16)

## Empty list for plots
pca_list <- list()
### Loop through vector of odd numbers and create plot i vs i + 1
for (i in 1:length(pcs_vec)) {

  pcX = paste0("PC", pcs_vec[i])
  pcY = paste0("PC", pcs_vec[i] + 1)
  pca_list[[i]] <- analysis_data_plots %>%
    pca_plot_supp(axis1 = pcX, axis2 = pcY,
             group_label = "participant_genetic_category",
             colours_order_df = prive_colours %>%
              filter(participant_genetic_category %in% analysis_data_plots$participant_genetic_category)) +
    theme(legend.position = "NONE")

  ### Also create plot specifically for the legend to add to the final grid space
  if (i == length(pcs_vec)) {
    pca_plot_legend <- analysis_data_plots %>%
      pca_plot_supp(axis1 = pcX, axis2 = pcY,
               group_label = "participant_genetic_category",
               colours_order_df = prive_colours %>%
                 filter(participant_genetic_category %in% analysis_data_plots$participant_genetic_category)) +
      guides(colour = guide_legend(ncol = 2,
                                   title="GIA group",
                                   override.aes = list(size = 3,
                                                       alpha = 1)))
    pca_legend <- get_legend(pca_plot_legend)
    pca_list[[i + 1]] <- as.ggplot(pca_legend)
  }
}
supp_figure_2_gg <- cowplot::plot_grid(plotlist = pca_list)
ggsave(paste0(plots_dir, "supp_figure2.png"), width = 15, height = 15)

############################
##### Supplementary Figure 3

supp_figure3a_gg <- analysis_data_plots %>%
  zoomed_bar_plot(group_label = "participant_genetic_category",
                  colours_order_df = prive_colours%>%
                 filter(participant_genetic_category %in% analysis_data_plots$participant_genetic_category)) +
  xlab("100kGP Probands") +
  theme(legend.position = "top") +
  guides(fill = guide_legend(title="GIA group"))

supp_figure3b_gg <- analysis_data_plots %>%
  pca_plot_supp(axis1 = "UMAP1", axis2 = "UMAP2",
           group_label = "participant_genetic_category",
           colours_order_df = prive_colours %>%
                 filter(participant_genetic_category %in% analysis_data_plots$participant_genetic_category)) +
  theme_void() +
  theme(legend.position = "NONE")

### Combine and output for Supplementary Figure 3
supp_figure_3_gg <- plot_grid(supp_figure3a_gg, supp_figure3b_gg, ncol = 1, rel_heights = c(0.5, 1))
ggsave(paste0(plots_dir, "supp_figure3.png"), width = 9, height = 10)

############################
##### Supplementary Figure 4

supp_figure4_gg <- analysis_data_plots %>%
  filter(participant_genetic_category == "Remaining participants") %>%
  admixture_plot("participant_type",
                 colours_order_df = prive_colours %>% filter(participant_genetic_category != "Remaining participants")) +
  ggtitle("Remaining participants") +
  guides(fill = guide_legend(title="Reference population"), colour = "none") +
  theme(strip.text = element_text(colour = "white"))

### Output for Supplementary Figure 4 as PDF.
ggsave(paste0(plots_dir, "supp_figure4.pdf"), width = 12, height = 3)

############################
##### Supplementary Figure 5

supp_figure5_gg <- analysis_data_plots %>%
  admixture_plot("participant_ethnic_category",
                 ethn_colours,
                 prive_colours %>% filter(participant_genetic_category != "Remaining participants")) +
  guides(fill = guide_legend(title="Reference population"), colour = "none") +
  ylab("Ethnic group") + theme(plot.margin = unit(c(0.25,5,0.25,0.25), "cm"))
ggsave(paste0(plots_dir, "supp_figure5.png"), width = 10, height = 8)

############################
##### Supplementary Figure 6

ethn_prive_stats <- analysis_data_plots %>%
  cat_var_stats("participant_ethnic_category",
                "participant_genetic_category")

supp_figure_4_gg <- ethn_prive_stats %>%
  concordance_plot("participant_ethnic_category", "participant_genetic_category", ethn_colours, prive_colours, 5) +
  xlab("GIA group") + ylab("Ethnic group")
ggsave(paste0(plots_dir, "supp_figure6.pdf"), width = 8, height = 8)

############################
##### Supplementary Figure 7

supp_figure7_gg <- analysis_data_plots %>%
  mutate(participant_gnomad_category = ifelse(participant_gnomad_category == "fin", "oth", participant_gnomad_category)) %>%
  admixture_plot("participant_gnomad_category",
                 order_rows_df = gnomad_order,
                 colours_order_df = prive_colours %>% filter(participant_genetic_category != "Remaining participants")) +
  guides(fill = guide_legend(title = "Reference population"), colour = "none") +
  ylab("gnomAD group")
ggsave(paste0(plots_dir, "supp_figure7.png"), width = 8, height = 8)

############################
##### Supplementary Figure 8

gnomad_prive_stats <- analysis_data_plots %>%
  mutate(participant_gnomad_category = ifelse(participant_gnomad_category == "fin", "oth", participant_gnomad_category)) %>%
  cat_var_stats("participant_gnomad_category",
                "participant_genetic_category")

supp_figure_8_gg <- gnomad_prive_stats %>%
  concordance_plot("participant_gnomad_category", "participant_genetic_category", gnomad_order, prive_colours, 5) +
  xlab("GIA group") + ylab("gnomAD group")
ggsave(paste0(plots_dir, "supp_figure8.pdf"), width = 6, height = 8)

################################################
### Supplementary Figure 9 created in PowerPoint...

############################
##### Supplementary Figure 10

supp_figure10a_gg <- analysis_data_plots %>%
  ggplot(aes(x = cPAV)) +
  geom_histogram(bins = 100, colour = "grey8", fill = "grey70") +
  theme_classic() +
  xlab("cPAVs per proband") + ylab("Count")

### Count of the number of rare / uncommon cPAVs
plot_frame <- data.frame(ylab = c(" ", " "),
                         label = c("popmax AF < 0.1%", "popmax 1% > AF > 0.1%\nBiallelic MOI only"),
                         count = c(1865089, 86570))

supp_figure10b_gg <- plot_frame %>%
  ggplot(aes(x = ylab, y = count, fill = label)) +
  geom_bar(stat = "identity", colour = "grey8") +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "bottom", axis.ticks.y = element_blank(),
        axis.line.y = element_blank(), axis.title.y = element_text(colour = "white")) +
  ylab("cPAVs (all)") + xlab("")

supp_figure10ab_gg <- plot_grid(supp_figure10a_gg, supp_figure10b_gg, ncol = 1, rel_heights = c(1, 0.4))

### cPAvs counts across probands
supp_figure10cd_gg <- cpavs_counts %>%
  mutate(n = ifelse(n > 10, ">10", as.character(n))) %>%
  group_by(n, freq) %>%
  summarise(counts = n()) %>%
  group_by(freq) %>%
  mutate(counts_prop = counts / sum(counts),
         n = factor(n, levels = c(as.character(seq(1,10)), ">10"))) %>%
  ggplot(aes(x = n, y = counts_prop)) +
  geom_bar(stat = "identity", fill = "grey60", colour = "grey8") +
  facet_wrap(~freq, ncol = 2) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.text = element_text(size = 11, hjust = 0)) +
  xlab("Count among 100kGP probands") + ylab("Proportion (%)")

plot_grid(supp_figure10ab_gg, supp_figure10cd_gg, ncol = 1, rel_heights = c(1.5, 1))
ggsave(paste0(plots_dir, "supp_figure10.pdf"), width = 8, height = 8)


############################
##### Supplementary Figure 11

### Create a dummy C/O for running the raincloud plot function...
dummy_colours_order_fgt <- data.frame(family_group_type = c("Singleton",
                                                            "Duo with other Biological Relative",
                                                            "Duo with Mother or Father",
                                                            "Trio with other Biological Relatives",
                                                            "Trio with Mother or Father and other Biological Relationship",
                                                            "Trio with Mother and Father",
                                                            "Families with more than three Participants"),
                                      colour = "grey35")

supp_figure_11a_gg <- analysis_data_plots %>%
  drop_na(family_group_type) %>%
  raincloud_plot("cPAV",
                 "family_group_type",
                 dummy_colours_order_fgt) +
  facet_grid(~Penetrance, scales = "free", space = "free") +
  theme(strip.background = element_blank(),
        plot.title = element_text(size = 8),
        axis.text.x = element_text(size = 6)) +
  ylab("cPAVs") + xlab("Family type") +
  ggtitle("Penetrance")

ggsave(paste0(plots_dir, "supp_figure11.png"), height = 6, width = 6)

#############################
##### Supplementary Figure 12

### CROH
roh_facet <- analysis_data_plots %>%
  box_facet_plot("c_roh_mb",
                 "participant_genetic_category",
                 prive_colours) +
  xlab("Cumulative >1Mb (Mb)") + ggtitle("ROH")

### Age
age_stats <- analysis_data_plots %>%
  cont_var_stats("participant_genetic_category", "age")
age_facet <- analysis_data_plots %>%
  distribution_facet_plot(age_stats,
                         "participant_genetic_category",
                         "age",
                         type = "density",
                         prive_colours) +
  xlab("Years") + ggtitle("Age") + ylab("")

### Karyotype
sex_stats <- analysis_data_plots %>%
  drop_na(karyotype) %>%
  cat_var_stats("participant_genetic_category", "karyotype")
sex_facet <- sex_stats %>%
  proportion_facet_plot("karyotype", "participant_genetic_category", prive_colours) +
  ylab("") + xlab("")  + ggtitle("Karyotype")

### Family group type
fgt_stats <- analysis_data_plots %>%
  cat_var_stats("participant_genetic_category", "family_group_type_broad")
fgt_facet <- fgt_stats %>%
  proportion_facet_plot("family_group_type_broad", "participant_genetic_category", prive_colours) +
  ylab("") + xlab("") + ggtitle("Family type")

### Penetrance
penetrance_stats <- analysis_data_plots %>%
  cat_var_stats("participant_genetic_category", "pen_show") %>%
  filter(pen_show != "n/a")
pen_facet <- penetrance_stats %>%
  proportion_facet_plot("pen_show", "participant_genetic_category", prive_colours) +
  ylab("") + xlab("") + ggtitle("Penetrance*")

### Parent affected
fma_stats <- analysis_data_plots %>%
  cat_var_stats("participant_genetic_category", "only_proband_affected")
fma_facet <- fma_stats %>%
  proportion_facet_plot("only_proband_affected", "participant_genetic_category", prive_colours) +
  ylab("") + xlab("") + ggtitle("Only FM affected")

### Panels size
panels_size_facet <- analysis_data_plots %>%
  box_facet_plot("unique_genes_panels_applied_mb",
                 "participant_genetic_category",
                 prive_colours) +
  theme(strip.text.y.left = element_blank()) +
  xlab("Size (Mb)") + ggtitle("Panel(s)")

######## Combine
supp_figure_12_gg <- cowplot::plot_grid(roh_facet, age_facet, sex_facet, fgt_facet, pen_facet, fma_facet, panels_size_facet,
          rel_rel_widths = c(2.3, 1, 0.75, 1.30, 0.75, 1, 1), nrow = 1)
ggsave(paste0(plots_dir, "supp_figure12.pdf"), width = 15, height = 11)

#############################
##### Supplementary Figure 13

## Adapted from Speidel et al. 2019
years_per_gen <- 28
filename <- paste0(analysis_dir, "relate_analysis/phased_panel_relate_subset_anc.pairwise.coal")
groups <- read.table(filename, nrow = 1)
t <- years_per_gen * t(as.matrix(read.table(filename, skip = 1, nrow = 1)))
num_pops <- round(sqrt(dim(read.table(filename, skip = 2))[1]))
pop_size <- data.frame(time = numeric(0), coal = numeric(0), pop_size = numeric(0), groups = numeric(0))

for (p1 in 1:num_pops) {
  for (p2 in 1:p1) {
    c <- as.matrix(read.table(filename, skip = (p1-1) * num_pops + p2 + 1, nrow = 1))[-c(1:2)]
    str <- rep(paste(groups[p1], " - ", groups[p2], sep = ""), length(c))
    pop_size <- rbind(pop_size, data.frame(time = t, coal = c, pop_size = 0.5/c, groups = str))
  }
}

#### Relate coal into pop-sizes...
relate_coal_all <- pop_size %>%
  separate(groups,
        into = c("Assignment", "participant_genetic_category_two"),
        sep = " - ", remove = FALSE,
        extra = "drop", fill = "right") %>%
  mutate(Assignment = gsub("_", " ", Assignment),
         participant_genetic_category_two = gsub("_", " ", participant_genetic_category_two),
         groups = gsub("_", " ", groups)) %>%
  adjust_genetic_labels()

supp_figure_13_gg <- relate_popsizes_plot(relate_coal_all, prive_colours)
ggsave(paste0(plots_dir, "supp_figure13.pdf"), height = 7, width = 8)

#############################
##### Supplementary Figure 14

supp_figure_14a_gg <- errorbar_plot(ukbb_path_count %>%
                                    filter(grepl("rare_cPAV", response) & 
                                    !grepl("nDeleterious", response)),
                                    "participant_genetic_category", NULL, prive_colours, "p.value") +
  ggtitle(paste0("cPAVs (Deleterious)\npopmax AF < 0.1%")) + 
  xlab("log(RR) vs Europe North West") +
  coord_cartesian(xlim = c(-0.75, 1.75)) + scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("") +
  theme_bw() +
  facet_wrap(~method) +
  theme(legend.position = "NONE",
        plot.title = element_text(size = 11, hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_blank())

supp_figure_14b_gg <- errorbar_plot(ukbb_path_count %>%
                                    filter(grepl("rare_cPAV", response) &
                                    grepl("nDeleterious", response)),
                                    "participant_genetic_category", NULL, prive_colours, "p.value") +
  ggtitle(paste0("cPAVs (Non-deleterious)\npopmax AF < 0.1%")) + 
  xlab("log(RR) vs Europe North West") +
  coord_cartesian(xlim = c(-0.75, 1.75)) + scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("") +
  theme_bw() +
  facet_wrap(~method) +
  theme(legend.position = "NONE",
        plot.title = element_text(size = 11, hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_blank())

plot_grid(supp_figure_14a_gg, supp_figure_14b_gg, ncol = 1)
ggsave(paste0(plots_dir, "supp_figure14.pdf"), width = 10, height = 9)

#############################
##### Supplementary Figure 15

supp_figure_15_gg <- errorbar_plot(ukbb_cpav_glm %>%
                                    filter(response == "GPcPAV"),
                                  "participant_genetic_category", NULL, prive_colours, "p.value") +
  ggtitle("GPcPAVs, n: 51,838") + 
  xlab("log(RR) vs Europe North West") +
  coord_cartesian(xlim = c(-0.55, 1.25)) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("GIA group") +
  theme_bw() + 
  theme(legend.position = "NONE",
        plot.title = element_text(size = 11),
        axis.ticks.y = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(paste0(plots_dir, "supp_figure15.pdf"), width = 5, height = 5)

#############################
##### Supplementary Figure 16

supp_figure_16a_gg <- analysis_data_plots %>%
  tier_pred_prop_table(
  "rare_GPcPAV",
  "CADD",
  "Deleterious",
  "participant_genetic_category"
) %>%
  group_bar_plot("tiered_prop", "participant_genetic_category", prive_colours) +
  ylab("Score >30\nDeleterious (%)") +
  xlab(" ") +
  ggtitle(paste0("CADD\n",
          "44,238",
          " GPcPAVs (popmax AF < 0.1%)")) +
  theme(plot.title = element_text(size = 11)) +
  coord_cartesian(ylim = c(0, 30))

supp_figure_16b_gg <- linear_points_plot(
  model1 = ukbb_cpav_glm %>%
    filter(response == "rare_GPcPAV"),
  model2 = ukbb_patho_pred_glm %>%
    filter(response == "rare_GPcPAV_CADD_Deleterious"),
  group_label = "participant_genetic_category",
  colours_order_df = prive_colours,
  y_pos = 0.6,
  x_pos = 0.65,
  p_label = "p.value"
) +
  ylab("\nlog(OR) Deleterious (%)\nvs Europe North West") +
  xlab("log(RR) GPcPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") + ggtitle(" ")

supp_figure_16ab_gg <- plot_grid(supp_figure_16a_gg, supp_figure_16b_gg, ncol = 2, rel_widths = c(1.8, 1))

supp_figure_16c_gg <- analysis_data_plots %>%
  tier_pred_prop_table(
  "rare_GPcPAV",
  "PAI3D",
  "Deleterious",
  "participant_genetic_category"
) %>%
  group_bar_plot("tiered_prop", "participant_genetic_category", prive_colours) +
  ylab("Score > 0.8\nDeleterious (%)") +
  xlab(" ") +
  ggtitle(paste0("PrimateAI-3D\n",
          "30,071",
          " missense GPcPAVs (popmax AF < 0.1%)")) +
  theme(plot.title = element_text(size = 11))  +
  coord_cartesian(ylim = c(0, 30))

supp_figure_16d_gg <- linear_points_plot(
  model1 = ukbb_cpav_glm %>%
    filter(response == "rare_GPcPAV"),
  model2 = ukbb_patho_pred_glm %>%
    filter(response == "rare_GPcPAV_PAI3D_Deleterious"),
  group_label = "participant_genetic_category",
  colours_order_df = prive_colours,
  y_pos = 0.3,
  x_pos = 0.65,
  p_label = "p.value"
) +
  ylab("\nlog(OR) Deleterious (%)\nvs Europe North West") +
  xlab("log(RR) GPcPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") + ggtitle(" ")

supp_figure_16cd_gg <- plot_grid(supp_figure_16c_gg, supp_figure_16d_gg, ncol = 2, rel_widths = c(1.8, 1))

supp_figure_16e_gg <- analysis_data_plots %>%
  tier_pred_prop_table(
  "rare_GPcPAV",
  "AM",
  "Deleterious",
  "participant_genetic_category"
) %>%
  group_bar_plot("tiered_prop", "participant_genetic_category", prive_colours) +
  ylab("Score > 0.564\nDeleterious (%)") +
  xlab(" ") +
  ggtitle(paste0("AlphaMissense\n",
          "32,170",
          " missense GPcPAVs (popmax AF < 0.1%)")) +
  theme(plot.title = element_text(size = 11))  +
  coord_cartesian(ylim = c(0, 30))

supp_figure_16f_gg <- linear_points_plot(
  model1 = ukbb_cpav_glm %>%
    filter(response == "rare_GPcPAV"),
  model2 = ukbb_patho_pred_glm %>%
    filter(response == "rare_GPcPAV_AM_Deleterious"),
  group_label = "participant_genetic_category",
  colours_order_df = prive_colours,
  y_pos = 0.3,
  x_pos = 0.65,
  p_label = "p.value"
) +
  ylab("\nlog(OR) Deleterious (%)\nvs Europe North West") +
  xlab("log(RR) GPcPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") +
  ggtitle(" ")

supp_figure_16ef_gg <- plot_grid(supp_figure_16e_gg, supp_figure_16f_gg, ncol = 2, rel_widths = c(1.8, 1))

### Full figure
supp_figure_16_gg <- plot_grid(supp_figure_16ab_gg, supp_figure_16cd_gg, supp_figure_16ef_gg, ncol = 1)
ggsave(paste0(plots_dir, "supp_figure16.pdf"), width = 8.75, height = 9.75)

#############################
##### Supplementary Figure 17

supp_figure_17a_gg <- errorbar_plot(ukbb_path_count %>%
                                    filter(grepl("rare_GPcPAV", response) & 
                                    !grepl("nDeleterious", response)),
                                    "participant_genetic_category", NULL, prive_colours, "p.value") +
  ggtitle(paste0("GPcPAVs (Deleterious)\npopmax AF < 0.1%")) + 
  xlab("log(RR) vs Europe North West") +
  coord_cartesian(xlim = c(-1, 1.5)) + scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("") +
  theme_bw() +
  facet_wrap(~method) +
  theme(legend.position = "NONE",
        plot.title = element_text(size = 11, hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_blank())

supp_figure_17b_gg <- errorbar_plot(ukbb_path_count %>%
                                    filter(grepl("rare_GPcPAV", response) &
                                    grepl("nDeleterious", response)),
                                    "participant_genetic_category", NULL, prive_colours, "p.value") +
  ggtitle(paste0("GPcPAVs (Non-deleterious)\npopmax AF < 0.1%")) + 
  xlab("log(RR) vs Europe North West") +
  coord_cartesian(xlim = c(-1, 1.5)) + scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("") +
  theme_bw() +
  facet_wrap(~method) +
  theme(legend.position = "NONE",
        plot.title = element_text(size = 11, hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_blank())

plot_grid(supp_figure_17a_gg, supp_figure_17b_gg, ncol = 1)
ggsave(paste0(plots_dir, "supp_figure17b.pdf"), width = 10, height = 9)

############################
#### Supplementary Figure 18

supp_figure_18_gg <- missing_gnomad_pavs_data %>%
  raincloud_plot2("missing_pavs_count",
                 "participant_genetic_category",
                 prive_colours) +
  facet_wrap(~dataset, nrow = 1) +
  xlab("GIA group") + ylab("Protein-altering variants\nmissing from gnomAD") +
  theme(panel.border = element_rect(linewidth = 0.5, colour = "grey80", fill = NA),
        strip.text = element_text(size = 10))
ggsave(paste0(plots_dir, "supp_figure18.png"), width = 8, height = 8)

###########################
#### Supplmentary Figure 19

supp_figure19b_gg <- filtering_conf_af_faf95_plot(covid_sizes = covid_sizes, step = 0.00001,
                                                  threshold = 0.001, start_freq = 0.00001, end_freq = 0.1,
                                                  colours_order_df = prive_colours) +
  ylab("\n\nProbability observed >0.1%") +
  xlab("True AF")
ggsave(paste0(plots_dir, "supp_figure19b.pdf"), height = 4.5, width = 4.5)

supp_figure19c_gg <- faf_filter_barplot(covid19_faf95_filter_tiering = covid19_faf95_filter_tiering,
                                        filter_test = "missing_cPAV_uncommon_FAF95",
                                        colours_order_df = prive_colours) +
  scale_y_continuous(breaks = c(0, 500, 1000)) +
  xlab("GIA group\nCOVID-19 cohort") +
  ylab("cPAVs (popmax AF = 0%)\nCOVID-19 cohort FAF95 > 0.1%") +
  ggtitle("")

## Write to plots
ggsave(paste0(plots_dir, "supp_figure19c.pdf"), width = 8, height = 4)

#############################
##### Supplementary Figure 20

## This code is a bit rough...
covariates_cat <- c("participant_genetic_category", "karyotype", "family_group_type", "penetrance", "only_proband_affected")
covariates_cont <- c("c_roh_mb", "n_roh", "age", "date_months", "unique_genes_panels_applied_mb", "panels_applied_count", "samtools_error_rate + samtools_reads_mapped_percentage + illumina_autosome_mean_coverage")

model_tidy %>%
  filter(grepl(paste(c(covariates_cat, covariates_cont), collapse = "|"), term)) -> model_tidy_collapse

model_tidy_collapse$covar <- NA_character_
for (covariate in covariates_cat) {
  model_tidy_collapse%<>% mutate(covar = ifelse(grepl(covariate, term), covariate, covar))
  model_tidy_collapse$term <- gsub(covariate, "", model_tidy_collapse$term)
}

### Take the tidied model and make the names nicer...
model_tidy_collapse %>%
  mutate(covar = ifelse(is.na(covar), term, covar)) %>%
  filter(!covar %in% c("phenotype_narrow", "handling_gmc_trust", "samtools_reads_mapped_percentage")) %>%
  mutate(covar_nice = case_when(covar == "consanguinity" ~ "Consanguinity vs Noncongsanguineous",
                           covar == "only_proband_affected" ~ "Only prob. affected\nvs No",
                           covar == "family_group_type" ~ "Family type\nvs Singleton", 
                           covar == "penetrance" ~ "Penetr.\nvs incompl.", 
                           covar == "karyotype" ~ "Karyotype\nvs XX",
                           covar == "n_roh" ~ "No. ROH",
                           covar == "c_roh_mb" ~ "cROH\nin Mb",
                           covar == "age" ~ "Age\nin years",
                           covar == "participant_genetic_category" ~ "GIA group\nvs Europe North West",
                           covar == "date_months" ~ "Date",
                           covar == "unique_genes_panels_applied_mb" ~  "Panel(s)\nin Mb",
                           covar == "panels_applied_count" ~  "No. Panels\napplied",
                           covar == "samtools_error_rate" ~ "Map error rate",
                           covar == "illumina_autosome_mean_coverage" ~ "Read depth",
                           covar == "samtools_reads_mapped_percentage" ~ "Percentage aligned reads (Samtools)"),
        #covar = ifelse(term == "", covar, term),
        cont_cat = ifelse(covar %in% covariates_cat, "Categorical", "Continuous")) -> model_plot_names

model_plot_names %>%
  filter(cont_cat == "Categorical") %>%
  drop_na(estimate) %>%
  add_sig_stars(p_label = "p.value") %>%
  ggplot(aes(y = term, x = estimate)) +
  geom_point() +
  geom_errorbar(aes(xmax = upper, xmin = lower), width = 0.1) +
  theme_bw() + 
  geom_text(aes(label = sig_stats), vjust = -0.045, position = position_dodge(width = 0.8)) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(rows = vars(covar_nice), scales = "free", space = "free") +
  coord_cartesian(clip = "off") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  xlab("log(Odds Ratio)") + ggtitle(" ") +
  ylab("") -> supp_figre16a_gg

model_plot_names %>%
  filter(cont_cat == "Continuous") %>%
  add_sig_stars(p_label = "p.value") %>%
  ggplot(aes(y = covar_nice, x = estimate)) +
  geom_point() +
  geom_errorbar(aes(xmax = upper, xmin = lower), width = 0.1) +
  theme_bw() + 
  geom_text(aes(label = sig_stats), vjust = -0.045, position = position_dodge(width = 0.8)) +
  theme(strip.background = element_blank(),
        panel.grid.major.y = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_grid(rows = vars(covar_nice), scales = "free", space = "free") +
  geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  coord_cartesian(xlim = c(-0.08, 0.08), clip = "off") +
  xlab("Change in log(Odds) per unit* increase") + ggtitle(" ") +
  ylab("") -> supp_figure16b_gg

lrt_tiered_yield %>%
  add_sig_stars(p_label = "Pr(>Chisq)") %>%
  mutate(sig_stats = ifelse(is.na(sig_stats), "", sig_stats)) %>%
  arrange(desc(`LR Chisq`)) %>%
  mutate(covar = case_when(covar == "c_roh_mb" ~ "cROH (Mb)",
                           covar == "n_roh" ~ "No. ROH",
                           covar == "only_proband_affected" ~ "Only prob. affected",
                           covar == "family_group_type" ~ "Family type", 
                           covar == "penetrance" ~ "Penetrance",
                           covar == "karyotype" ~ "Karyotype",
                           covar == "age" ~ "Age",
                           covar == "samtools_error_rate" ~ "Map error rate",
                           covar == "illumina_autosome_mean_coverage" ~ "Read depth",
                           covar == "phenotype_narrow" ~ "Phenotype",
                           covar == "participant_genetic_category" ~ "GIA group",
                           covar == "date_months" ~ "Date",
                           covar == "panels_applied_count" ~ "No. panels applied",
                           covar == "unique_genes_panels_applied_mb" ~ "Panel(s) Mb"),
          covar = factor(covar, levels = rev(covar)),
          facet_null = "null\nnull",
          shown = ifelse(covar %in% c("GIA group", "Phenotype"), "No", "Yes")) %>%
  drop_na(covar) %>%
  ggplot(aes(y = covar, x = `LR Chisq`, colour = shown, fill = shown)) +
  geom_bar(stat = "identity", width = 0.25) +
  theme_bw() +
  facet_grid(rows = vars(facet_null)) +
  theme(panel.grid.major.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), legend.position = "NONE",
        plot.title = element_text(size = 11), strip.text = element_text(colour = "white"),
        axis.line = element_blank(),
        panel.border = element_rect(size = 0.3, colour = "grey40"),
        strip.background = element_blank()) +
  geom_text(aes(label = paste0(covar, "\n    ", sig_stats), x = `LR Chisq`, y = covar),
                hjust = -0.075, size = 2.8, colour = "black") +
  scale_colour_manual(values = c("darkred", "grey40")) +
  scale_fill_manual(values = c("darkred", "grey40")) + ylab("") +
  ggtitle("Analysis of Deviance (Type II)") + xlab("Likelihood Ratio") +
  coord_cartesian(xlim = c(0, 3250), clip = "off") -> supp_figure16c_gg

### Arrange them all together + add a title
plot_grid(supp_figure16b_gg, supp_figure16c_gg, ncol = 1, rel_heights = c(0.6, 1)) -> supp_figure16bc_gg
plot_grid(supp_figre16a_gg, supp_figure16bc_gg, ncol = 2, rel_widths = c(1.8, 1)) -> supp_figure16_gg
ggpubr::annotate_figure(supp_figure16_gg, top = ggpubr::text_grob("Case solved with GPcPAV(s)", 
                        size = 14))
ggsave(paste0(plots_dir, "supp_figure20.pdf"), width = 12, height = 15)

############################
#### Supplementary Figure 21

supp_figure_21a_gg <- analysis_data %>%
  cat_var_stats("participant_genetic_category", "DIAGNOSED_GPcPAV_ACMG_Pathogenic") %>%
  filter(DIAGNOSED_GPcPAV_ACMG_Pathogenic == "Yes") %>%
  mutate(proportion = proportion * 100) %>%
  group_bar_plot("proportion", "participant_genetic_category", prive_colours) +
  ylab("Yield (%)") +
  xlab(" ") +
  ggtitle("Case solved using GPcPAV(s)") +
  theme(plot.title = element_text(size = 11))

supp_figure_21b_gg <-
  linear_points_plot(model1 = ukbb_cpav_glm %>%
                      filter(response == "GPcPAV"),
                     model2 = ukbb_diagnostic_yield_glm %>%
                       filter(response == "DIAGNOSED_GPcPAV_ACMG_Pathogenic"),
                     group_label = "participant_genetic_category",
                     colours_order_df = prive_colours,
                     y_pos = -0.75, x_pos = 0.45, p_label = "p.adjust") +
  ylab("\nlog(OR) Yield (%)\nvs Europe North West") +
  xlab("log(RR) GPcPAVs\nvs Europe North West") +
  theme(legend.position = "NONE") +
  ggtitle(" ")

supp_figure_21_gg <- plot_grid(supp_figure_21a_gg, supp_figure_21b_gg, ncol = 2, rel_widths = c(1.9, 1))
ggsave(paste0(plots_dir, "supp_figure21.pdf"), width = 8.75, height = 3)

##########################
## Supplementary Figure 22

supp_figure_22_gg <- errorbar_plot(ukbb_vus_glm %>%
                                      filter(response == "GPcPAV_VUS_or_noclass"),
                                    "participant_genetic_category", NULL, prive_colours, "p.value") +
  ggtitle("GPcPAVs (unclassified or VUS), n: 47,021") + 
  xlab("log(RR) vs Europe North West") +
  coord_cartesian(xlim = c(-0.85, 1.45)) +
  scale_x_continuous(breaks = c(-2, -1, 0, 1, 2)) +
  ylab("GIA group") +
  theme_bw() +
  facet_wrap(~variable, nrow = 1) +
  theme(legend.position = "NONE",
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 11, hjust = 0.5),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(paste0(plots_dir, "supp_figure22.pdf"), width = 8, height = 5)

###########################################################################
###########################################################################
###                                                                     ###
###                         END OF SCRIPT                               ###
###                                                                     ###
###########################################################################
###########################################################################