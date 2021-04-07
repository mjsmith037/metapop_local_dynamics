library(magrittr)
library(tidyverse)
library(JuliaCall)

julia_setup()
julia_source("../../Code/multipop_MANTIS.jl")

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")

## structural parameters
set.seed(0)
struct <- c(2,2,2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))
n_patches <- 3
# initial_conditions <- rep(0.00001, prod(struct) * n_patches)
initial_conditions <- runif(prod(struct) * n_patches) %>%
  matrix(n_patches, prod(struct)) %>%
  apply(1, . %>% {. / (5 * sum(.))}) %>%
  t()
## dynamical parameters
gamma <- 0.66                            # partial cross-protective immunity (cpi)
sigma <- 16                              # recovery rate
movement_rate <- 0.025
mu <- 0.125
chi <- matrix(c(-movement_rate, 0, 0,
                0, -movement_rate, 0,
                movement_rate, movement_rate, 0), 3, 3)
beta <- c(5, 3, 5) * (sigma + mu - diag(chi))     # Infection rate


## integration parameters
maxtime <- 1000
tsteps <- round(0.85 * maxtime):maxtime

## run the simulation
timeseries <- julia_call("runMANTIS", strainstructure=struct, tstep=tsteps, tmax=maxtime,
                         beta=beta, gamma=gamma, sigma=sigma, mu=mu, chi=chi,
                         initvals=initial_conditions)$timeseries %>%
  set_colnames(c(expand.grid(str_c("Patch_", 1:n_patches),
                             str_c("Strain_", apply(strains, 1, str_c, collapse="")),
                             c("y", "z", "w")) %>% apply(1, str_c, collapse="."),
                 "time")) %>%
  gather("details", "prevalence", -time) %>%
  separate(details, c("patch", "strain", "equation"), "\\.") %>%
  mutate(variable = factor(equation, levels=c("y", "z", "w"),
                           labels=c("currently infectious (y)",
                                    "specific immunity (z)",
                                    "cross-reactive immunity (w)")),
         patch = factor(patch, levels=c("Patch_1", "Patch_2", "Patch_3"),
                             labels=c("Patch A", "Patch B", "Patch C"))) %>%
  as_tibble()

ggplot(timeseries %>% filter(strain == "Strain_122", equation=="y")) +
  aes(colour=patch, y=prevalence, x=time) +
  geom_line(size=0.75) +
  facet_grid_sc(~patch) +
  scale_x_continuous(expand=c(0,0), breaks=c(860, 900, 940, 980)) +
  scale_colour_manual(values=my_cols) +
  ylab("Proportion infectious (y)") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave("dynamics-hierarchy.pdf", width=8, height=2.5)
