library(magrittr)
library(tidyverse)
library(tidygraph)
library(JuliaCall)
library(ggraph)

julia_setup()
julia_source("../../code/multipop_MANTIS.jl")

## plotting parameters
my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")

## structural parameters
set.seed(0)
struct <- c(2,2,2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))
n_patches <- 2
initial_conditions <- runif(prod(struct) * n_patches) %>%
  matrix(n_patches, prod(struct)) %>%
  apply(1, . %>% {. / (5 * sum(.))}) %>%
  matrix(n_patches, prod(struct), byrow=TRUE)

## dynamical parameters
gamma <- 0.66           # partial cross-protective immunity (cpi)
sigma <- 16             # recovery rate
mu <- 0.125             # death rate
beta <- c(48, 80)       # transmission rate
movement_rate <- 0.025

## integration parameters
maxtime <- 1000
tsteps <- round(0.85 * maxtime):maxtime

## run the simulation
timeseries <- bind_rows(
  julia_call("runMANTIS", strainstructure=struct, tstep=tsteps, tmax=maxtime,
             beta=beta, gamma=gamma, sigma=sigma, mu=mu,
             chi=diag(0, 2), initvals=initial_conditions)$timeseries %>%
    set_colnames(c(expand.grid(str_c("Patch_", 1:n_patches),
                               str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                               c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("patch", "strain", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("z", "w", "y"),
                             labels=c("specific immunity (z)",
                                      "cross-reactive immunity (w)",
                                      "currently infectious (y)")),
           network = "A      B"),
  julia_call("runMANTIS", strainstructure=struct, tstep=tsteps, tmax=maxtime,
             beta=beta, gamma=gamma, sigma=sigma, mu=mu,
             chi=matrix(c(-movement_rate, 0, movement_rate, 0), 2, 2),
             initvals=initial_conditions)$timeseries %>%
    set_colnames(c(expand.grid(str_c("Patch_", 1:n_patches),
                               str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                               c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("patch", "strain", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("z", "w", "y"),
                             labels=c("specific immunity (z)",
                                      "cross-reactive immunity (w)",
                                      "currently infectious (y)")),
           network = "A -> B"),
  julia_call("runMANTIS", strainstructure=struct, tstep=tsteps, tmax=maxtime,
             beta=beta, gamma=gamma, sigma=sigma, mu=mu,
             chi=matrix(c(0, movement_rate, 0, -movement_rate), 2, 2),
             initvals=initial_conditions)$timeseries %>%
    set_colnames(c(expand.grid(str_c("Patch_", 1:n_patches),
                               str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                               c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("patch", "strain", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("z", "w", "y"),
                             labels=c("specific immunity (z)",
                                      "cross-reactive immunity (w)",
                                      "currently infectious (y)")),
           network = "A <- B")
) %>%
  mutate(network = factor(network, levels=c("A      B", "A -> B", "A <- B")),
         variable = factor(variable, levels=c("currently infectious (y)",
                                              "specific immunity (z)",
                                              "cross-reactive immunity (w)")),
         patch = factor(patch, levels=c("Patch_1", "Patch_2"),
                             labels=c("Patch A", "Patch B"))) %>%
  as_tibble()

ggplot(timeseries %>% filter(strain == "Strain_111", equation=="y")) +
  aes(colour=patch, y=prevalence, x=time) +
  geom_line(size=0.75) +
  facet_wrap(~network) +
  scale_x_continuous(expand=c(0,0), breaks=c(860, 900, 940, 980)) +
  scale_colour_manual(values=my_cols) +
  ylab("Proportion infectious (y)") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave("two-pop-comparison.pdf", width=8, height=2.5)
