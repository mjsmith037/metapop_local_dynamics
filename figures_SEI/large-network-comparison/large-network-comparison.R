library(igraph)
library(magrittr)
library(tidyverse)
library(tidygraph)
library(JuliaCall)
library(kableExtra)
library(assertthat)
library(patchwork)
library(e1071)
library(ggraph)
library(janitor)
library(ggbeeswarm)
library(ggforce)

my_cols <- c("er"="#016394", "sb"="#f74700", "tr"="#b6003b", "ba"="#005342", "ws"="#1b264f",
             "Random"="#016394", "Modular"="#f74700", "Tree"="#b6003b", "Scale-free"="#005342", "Small-world"="#1b264f",
             "Erdős-Rényi-Gilbert"="#016394", "Stochastic Block"="#f74700", "Tree"="#b6003b", "Barabasi-Albert"="#005342", "Watts-Strogatz"="#1b264f")

filler_plot <- ggplot() + theme(panel.background=element_blank())

## load the networks
make_network_plots <- function() {
  source("../../Code/random_tree_networks.R")
  source("../../Code/directed_watts_strogatz_networks.R")
  set.seed(0)
  n_patches <- 25
  example_er <- play_erdos_renyi(n_patches, m=90) %>% mutate(net_type = "er") %>%
    activate(edges) %>% mutate(net_type = "er")
  example_sb <- play_blocks(n_patches, c(floor(n_patches / 2), ceiling(n_patches / 2)),
                            matrix(c(0.28, 0.03, 0.03, 0.28), 2, 2)) %>% mutate(net_type = "sb") %>%
    activate(edges) %>% mutate(net_type = "sb")
  example_ba <- play_barabasi_albert(n_patches, 0.5, 4) %>% mutate(net_type = "ba") %>%
    activate(edges) %>% distinct(to, from) %>% filter(to != from) %>%
    activate(edges) %>% mutate(net_type = "ba")
  example_tr <- simple_tree(n_patches, 1, 6, 12) %>% mutate(net_type = "tr") %>%
    activate(edges) %>% mutate(net_type = "tr")
  example_ws <- directed_watts_strogatz(n_patches, 4, 0.5, 90) %>% mutate(net_type = "ws") %>%
    activate(edges) %>% mutate(net_type = "ws")

  layout_matrix <- rbind(example_er %>% as.igraph() %>% layout_nicely(),
                         example_sb %>% as.igraph() %>% layout_nicely(),
                         example_ba %>% as.igraph() %>% layout_nicely(),
                         example_ws %>% as.igraph() %>% layout_nicely(),
                         example_tr %>% as.igraph() %>% layout_as_tree()) %>%
    set_colnames(letters[24:25]) %>%
    as_tibble()

  nets <- bind_graphs(example_er, example_sb, example_ba, example_ws, example_tr) %>%
    mutate(net_type = factor(net_type, levels=c("er", "sb", "ws", "tr", "ba"),
                             labels=c("Random", "Modular", "Small-world", "Tree", "Scale-free"))) %>%
    activate(nodes) %>%
    mutate(net_type = factor(net_type, levels=c("er", "sb", "ws", "tr", "ba"),
                             labels=c("Random", "Modular", "Small-world", "Tree", "Scale-free"))) %>%
    {ggraph(., layout=layout_matrix) +
        geom_edge_link(aes(colour=net_type),
                       end_cap=circle(5, "pt"),
                       edge_width=0.33,
                       edge_alpha=0.5,
                       arrow=arrow(angle=30, length=unit(3, "pt"), type="closed")) +
        geom_node_point(aes(colour=net_type), size=1) +
        facet_nodes(~net_type, scales="free", nrow=1) +
        scale_colour_manual(values=my_cols) +
        scale_edge_colour_manual(values=my_cols) +
        theme_bw() +
        theme(panel.spacing.x=unit(3, "pt"),
              axis.text=element_blank(),
              axis.title=element_blank(),
              axis.ticks=element_blank(),
              panel.grid=element_blank(),
              legend.position="none")}
  return(nets)
}

make_metrics_plots <- function(file_name) {
  all_sims <- read_csv(file_name) %>% rename_with(~str_remove(., "prevalence_"))

  plot_data <- all_sims %>%
    mutate(var = ifelse(summarise_dynamics <= 1, NA, var),
           ave_cycle = ifelse(summarise_dynamics > 1, ave_cycle, NA),
           var_cycle = ifelse(summarise_dynamics > 1, var_cycle, NA)) %>%
    filter(summarise_dynamics != 1) %>% # remove unconverged
    pivot_longer(!any_of(c("patch", "summarise_dynamics", "replicate", "net_type", "strain")),
                 names_to="metric", values_to="value") %>%
    filter(across(any_of("strain"), ~equals(., 1))) %>%
    mutate(net_type = factor(net_type, levels=rev(c("er", "sb", "ws", "tr", "ba")),
                             labels=rev(c("Random", "Modular", "Small-world", "Tree", "Scale-free"))),
           metric = factor(metric, levels=c("sum", "var", "ave_cycle", "var_cycle")))

  met1 <- ggplot(plot_data %>% filter(metric == "sum") %>% mutate(metric = "Cumulative Prevalence")) +
    aes(x=value, y=net_type, colour=net_type) +
    geom_quasirandom(width=0.25, alpha=0.15, varwidth=TRUE, groupOnX=FALSE) +
    scale_colour_manual(values=my_cols) +
    scale_x_log10() +
    xlab("Cumulative Prevalence") +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y=element_blank())

  met2 <- ggplot(plot_data %>% filter(metric == "var") %>% mutate(metric = "Variance in Prevalence")) +
    aes(x=value, y=net_type, colour=net_type) +
    geom_quasirandom(width=0.25, alpha=0.15, varwidth=TRUE, groupOnX=FALSE) +
    scale_colour_manual(values=my_cols) +
    scale_x_log10() +
    xlab("Variance in Prevalence") +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  met3 <- ggplot(plot_data %>% filter(metric == "ave_cycle") %>% mutate(metric = "Average Cycle Length")) +
    aes(x=value, y=net_type, colour=net_type) +
    geom_quasirandom(width=0.25, alpha=0.15, varwidth=TRUE, groupOnX=FALSE) +
    scale_colour_manual(values=my_cols) +
    scale_x_log10() +
    xlab("Average Cycle Length") +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  met4 <- ggplot(plot_data %>% filter(metric == "var_cycle") %>% mutate(metric = "Variance in Cycle Length")) +
    aes(x=value, y=net_type, colour=net_type) +
    geom_quasirandom(width=0.25, alpha=0.15, varwidth=TRUE, groupOnX=FALSE) +
    scale_colour_manual(values=my_cols) +
    scale_x_log10() +
    xlab("Variance in Cycle Length") +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

  return(list(met1, met2, met3, met4))
}

make_propdyn_plot <- function(file_name) {
  all_sims <- read_csv(file_name) %>% rename_with(~str_remove(., "prevalence_"))

  prop_dynamics <- all_sims %>%
    mutate(var = ifelse(summarise_dynamics <= 1, 0, var),
           ave_cycle = ifelse(summarise_dynamics > 1, ave_cycle, NA),
           var_cycle = ifelse(summarise_dynamics > 1, var_cycle, NA),
           dynamics = case_when(summarise_dynamics == -1 ~ "extinct",
                                summarise_dynamics ==  0 ~ "stable",
                                summarise_dynamics ==  1 ~ "unconverged",
                                TRUE ~ "cycles")) %>%
    mutate(net_type = factor(net_type, levels=c("er", "sb", "ws", "tr", "ba"),
                             labels=c("Random", "Modular", "Small-world", "Tree", "Scale-free"))) %>%
    group_by(net_type, replicate) %>%
    group_modify(~tabyl(., dynamics, )) %>%
    pivot_wider(-n, names_from=dynamics, values_from=percent) %>%
    mutate(across(everything(), ~replace_na(., 0))) %>%
    arrange(desc(cycles), desc(stable), desc(extinct), desc(unconverged)) %>%
    group_by(net_type) %>%
    mutate(replicate = 1:n()) %>%
    ungroup() %>%
    pivot_longer(!c(net_type, replicate), names_to="dynamics", values_to="percent") %>%
    mutate(dynamics = dynamics %>% str_to_title() %>% str_replace_all("Cycles", "Cycles or Chaos") %>%
             factor(levels=c("Cycles or Chaos", "Stable", "Extinct", "Unconverged") %>% rev()))

  dynamics <- ggplot(prop_dynamics) +
    aes(x=replicate, y=percent, fill=dynamics, colour=dynamics) +
    geom_col(position="stack") +
    facet_wrap(~net_type, nrow=1, scales="free_x") +
    scale_fill_manual(values=my_cols %>% set_names(NULL) %>% .[1:4] %>% rev(), aesthetics=c("colour", "fill")) +
    ylab("Proportion of Patches") +
    coord_cartesian(expand=FALSE) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          legend.title=element_blank(),
          plot.margin=margin(0, 5.5, 5.5, 5.5))

  return(dynamics)
}

# single strain main parameterization
metrics_plots <- make_metrics_plots("../../Results/LNS_sigma=13/large-network-simulations.csv")
ggsave(gridExtra::arrangeGrob(metrics_plots[[1]] + facet_zoom(xlim = c(0.25, 0.65), zoom.size=1),
                              metrics_plots[[2]],
                              metrics_plots[[3]] + facet_zoom(xlim = c(0.38, 0.54), zoom.size=1),
                              metrics_plots[[4]] + facet_zoom(xlim = c(-1.75, 0), zoom.size=1),
                              filler_plot,
                              layout_matrix=matrix(c(2,3,4,5,2,6,4,5), 2, 4, byrow=TRUE),
                              widths=c(1.6,1,1,1), heights=c(4,4)),
       file="main_metrics.pdf", width=9, height=5, device=cairo_pdf)

ggsave((make_network_plots() + scale_colour_manual(values=rep("black", 5), aesthetics=c("colour", "edge_colour"))) /
         make_propdyn_plot("../../Results/LNS_sigma=13/large-network-simulations.csv") +
         theme(strip.text=element_blank(), strip.background=element_blank()) +
         plot_layout(heights=c(2, 2.25)),
       filename="main_propdyn.pdf", width=8, height=4.25, device=cairo_pdf)

# single strain alternate parameterization
metrics_plots <- make_metrics_plots("../../Results/LNS_sigma=40/large-network-simulations.csv")
ggsave(gridExtra::arrangeGrob(metrics_plots[[1]] + facet_zoom(xlim = c(0.15, 0.55), zoom.size=1),
                              metrics_plots[[2]],
                              metrics_plots[[3]] + facet_zoom(xlim = c(0.303, 0.43), zoom.size=1),
                              metrics_plots[[4]] + facet_zoom(xlim = c(-1.6, -0.1), zoom.size=1),
                              filler_plot,
                              layout_matrix=matrix(c(2,3,4,5,2,6,4,5), 2, 4, byrow=TRUE),
                              widths=c(1.6,1,1,1), heights=c(4,4)),
       file="alternate_metrics.pdf", width=9, height=5, device=cairo_pdf)
ggsave(make_propdyn_plot("../../Results/LNS_sigma=40/large-network-simulations.csv"),
       filename="alternate_propdyn.pdf", width=8, height=2.5, device=cairo_pdf)
