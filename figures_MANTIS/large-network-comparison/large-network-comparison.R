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

make_figs <- function(data_file, fig_id) {
  load(data_file)
  tmp <- all_sims %>%
    lapply(function(x) lapply(x, function(y) {
      if (length(y) == 0) return(NULL)
      y$timeseries %>% mutate(net_type=y$network$net_type)
    }) %>% bind_rows()) %>% bind_rows() %>%
    rename(summarise_dynamics=num_nonzero) %>%
    mutate(net_type = factor(net_type, levels=c("er", "sb", "ws", "tr", "ba"),
                             labels=c("Random", "Modular", "Small-world", "Tree", "Scale-free")))
  met1 <- ggplot(tmp %>%
           pivot_longer(c(ave_inf, ave_res, ave_cycle)) %>%
           filter(name == "ave_inf")) +
    aes(x=value, y=net_type, colour=net_type) +
    geom_quasirandom(width=0.25, alpha=0.15, varwidth=TRUE, groupOnX=FALSE) +
    scale_colour_manual(values=my_cols) +
    xlab("Average number infectious") +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y=element_blank())
  met2 <- ggplot(tmp %>%
                   pivot_longer(c(ave_inf, ave_res, ave_cycle)) %>%
                   filter(name == "ave_res")) +
    aes(x=value, y=net_type, colour=net_type) +
    geom_quasirandom(width=0.25, alpha=0.15, varwidth=TRUE, groupOnX=FALSE) +
    scale_colour_manual(values=my_cols) +
    xlab("Average number partially immune") +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y=element_blank())
  met3 <- ggplot(tmp %>%
                   pivot_longer(c(ave_inf, ave_res, ave_cycle)) %>%
                   filter(name == "ave_cycle")) +
    aes(x=value, y=net_type, colour=net_type) +
    geom_quasirandom(width=0.25, alpha=0.15, varwidth=TRUE, groupOnX=FALSE) +
    scale_colour_manual(values=my_cols) +
    scale_x_log10() +
    xlab("Average time between peaks") +
    theme_bw() +
    theme(legend.position="none",
          axis.title.y=element_blank())

  prop_dynamics <- tmp %>%
    group_by(net_type) %>%
    mutate(replicate = rep(1:(n() / 100), each=100)) %>%
    ungroup() %>%
    mutate(dynamics = case_when(summarise_dynamics == -1 ~ "extinct",
                                summarise_dynamics ==  0 ~ "stable",
                                summarise_dynamics ==  1 ~ "unconverged",
                                TRUE ~ "cycles")) %>%
    group_by(net_type, replicate) %>%
    group_modify(~tabyl(., dynamics, )) %>%
    pivot_wider(-n, names_from=dynamics, values_from=percent) %>%
    mutate(across(everything(), ~replace_na(., 0))) %>%
    arrange(desc(cycles), desc(stable), desc(unconverged)) %>%
    group_by(net_type) %>%
    mutate(replicate = 1:n()) %>%
    ungroup() %>%
    pivot_longer(!c(net_type, replicate), names_to="dynamics", values_to="percent") %>%
    mutate(dynamics = dynamics %>% str_to_title() %>% str_replace_all("Cycles", "Cycles or Chaos") %>%
             factor(levels=c("Cycles or Chaos", "Stable", "Extinct", "Unconverged") %>% rev()))

  ggplot(prop_dynamics) +
    aes(x=replicate, y=percent, fill=dynamics, colour=dynamics) +
    geom_col(position="stack") +
    facet_wrap(~net_type, nrow=1, scales="free_x") +
    scale_fill_manual(values=my_cols %>% set_names(NULL) %>% .[1:4] %>% rev(), aesthetics=c("colour", "fill"),
                      breaks=rev(c("Cycles or Chaos", "Stable", "Extinct", "Unconverged"))) +
    ylab("Proportion of Patches") +
    coord_cartesian(expand=FALSE) +
    theme_bw() +
    theme(axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          legend.position="bottom",
          legend.title=element_blank(),
          plot.margin=margin(0, 5.5, 5.5, 5.5))
  ggsave(filename=str_glue("{fig_id}_propdyn.pdf"), width=8, height=2.5)

  return(list(met1, met2, met3))
}
metrics_plots <- make_figs("large-network-simulations_g0.66,s8,m0.05,d0.01.RData", "main")
ggsave(gridExtra::arrangeGrob(metrics_plots[[1]],
                              metrics_plots[[2]] +
                                theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                              metrics_plots[[3]] +
                                theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                              widths=c(1.3,1,1)),
       file="main_metrics.pdf", width=9, height=3, device=cairo_pdf)

metrics_plots <- make_figs("large-network-simulations_g0.55,s32,m0.05,d0.01.RData", "alternate")
ggsave(gridExtra::arrangeGrob(metrics_plots[[1]],
                              metrics_plots[[2]] +
                                theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                              metrics_plots[[3]] +
                                theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                              widths=c(1.3,1,1)),
       file="alternate_metrics.pdf", width=9, height=3, device=cairo_pdf)
