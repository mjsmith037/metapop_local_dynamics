library(magrittr)
library(tidyverse)
library(JuliaCall)

julia_setup()
julia_source("../../Code/multipop_MANTIS.jl")

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")

## structural parameters
set.seed(0)
struct <- c(2, 2, 2)
strains <- expand.grid(lapply(struct, seq, from=1, by=1) %>%
                         setNames(str_c("locus", 1:length(.))))
n_patches <- 4

# initial_conditions <- rep(0.00001, prod(struct) * n_patches)
initial_conditions <- runif(prod(struct) * n_patches) %>%
  matrix(n_patches, prod(struct)) %>%
  apply(1, . %>% {. / (5 * sum(.))}) %>%
  matrix(n_patches, prod(struct) , byrow=TRUE)

## dynamical parameters
gamma <- 0.66
mu <- 0.125 #c(0, 0.05, 0.05, 0.05, 0.05, 0.1)
movement_rate <- 0.025
chi <- matrix(c(rep(c(-movement_rate, rep(0, n_patches-1), movement_rate), n_patches - 1), 0), n_patches, n_patches)
# chi <- matrix(c(rep(c(-movement_rate, movement_rate, rep(0, n_patches-2), movement_rate), n_patches - 1), 0), n_patches, n_patches)
# chi[n_patches,1] <- movement_rate; chi[n_patches,n_patches] <- -movement_rate
# chi <- matrix(0, n_patches, n_patches)
sigma <- 16
beta <- 48#sample(15:80, 26)

r_0 <- function(beta, sigma, gamma, mu, w) beta * (1 - gamma * w) / (sigma + mu)

## integration parameters
maxtime <- 1000
tsteps <- round(0.8 * maxtime):maxtime

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
         # beta = beta[as.integer(str_extract(patch, "\\d+"))], sigma=sigma, gamma=gamma, mu=mu,
         patch = factor(patch,
                             levels=str_c("Patch_", 1:n_patches),
                             labels=str_c("Patch ", LETTERS[1:n_patches]))) %>%
  as_tibble()

ggplot(timeseries %>% filter(strain == "Strain_111", equation=="y")) +
  aes(colour=patch, y=prevalence, x=time) +
  geom_line(size=0.75) +
  facet_wrap(~patch, nrow=1) +
  scale_x_continuous(expand=c(0,0), breaks=c(860, 900, 940, 980)) +
  scale_colour_manual(values=my_cols) +
  ylab("Proportion infectious (y)") +
  theme_bw() +
  theme(legend.position="bottom",
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave("dampening-cycles.pdf", width=10.5, height=2.5)

# timeseries %>%
#   select(-variable) %>%
#   pivot_wider(names_from=equation, values_from=prevalence) %>%
#   mutate(R0 = r_0(beta, sigma, gamma, mu, w)) %>%
#   group_by(patch) %>%
#   summarise(min(R0)) %>%
#   print(n=26)

timeseries %>%
  filter(strain == "Strain_111") %>%
  group_by(patch, equation) %>%
  summarise(prevalence = mean(prevalence), .groups="drop") %>%
  pivot_wider(names_from=equation, values_from=prevalence)
