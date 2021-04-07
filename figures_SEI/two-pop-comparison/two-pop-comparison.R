library(magrittr)
library(tidyverse)
library(JuliaCall)
library(scales)

julia_setup()
julia_source("../../Code/multipop_SI.jl")

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342")

mysqrt_trans <- function() {
  trans_new("mysqrt",
            transform = base::sqrt,
            inverse = function(x) ifelse(x<0, 0, x^2),
            domain = c(0, Inf))
}

## structural parameters
set.seed(0)
n_patches <- 2
# initial_conditions <- rep(0.00001, prod(struct) * n_patches)
initial_conditions <- matrix(c(5, 1, 1), n_patches, 3, byrow=TRUE)

## dynamical parameters
r <- 0.5             # intrinsic growth rate
K <- c(15, 5)       # carrying capacity
beta <- 80           # infection rate
sigma <- 13          # 1 / sigma = latency period
mu <- 0.5            # intrinsic death rate
nu <- 73             # disease induced mortality
movement_rate <- 0.1
chi <- matrix(c(rep(c(-movement_rate, rep(0, n_patches-1), movement_rate), n_patches - 1), 0),
              n_patches, n_patches)

## integration parameters
maxtime <- 10000
tsteps <- round(0.9 * maxtime):maxtime

## run the simulation
timeseries <- bind_rows(
  julia_call("runsimulation", tstep=tsteps, tmax=maxtime,
             r=r, K=K, beta=beta, sigma=sigma, mu=mu, nu=nu,
             chi=diag(0, 2),
             initvals=initial_conditions) %>%
    set_colnames(c(expand.grid(str_c("Patch_", 1:n_patches),
                               c("S", "E", "I")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("patch", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("S", "E", "I"),
                             labels=c("Susceptible", "Exposed", "Infectious")),
           patch = factor(patch, levels=str_c("Patch_", 1:n_patches),
                               labels=str_c("Patch ", LETTERS[1:n_patches]))) %>%
    as_tibble() %>%
    mutate(network = "A      B"),
  julia_call("runsimulation", tstep=tsteps, tmax=maxtime,
             r=r, K=K, beta=beta, sigma=sigma, mu=mu, nu=nu,
             chi=matrix(c(-movement_rate, 0, movement_rate, 0), 2, 2),
             initvals=initial_conditions) %>%
    set_colnames(c(expand.grid(str_c("Patch_", 1:n_patches),
                               c("S", "E", "I")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("patch", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("S", "E", "I"),
                             labels=c("Susceptible", "Exposed", "Infectious")),
           patch = factor(patch, levels=str_c("Patch_", 1:n_patches),
                               labels=str_c("Patch ", LETTERS[1:n_patches]))) %>%
    as_tibble() %>%
    mutate(network = "A -> B"),
  julia_call("runsimulation", tstep=tsteps, tmax=maxtime,
             r=r, K=K, beta=beta, sigma=sigma, mu=mu, nu=nu,
             chi=matrix(c(0, movement_rate, 0, -movement_rate), 2, 2),
             initvals=initial_conditions) %>%
    set_colnames(c(expand.grid(str_c("Patch_", 1:n_patches),
                               c("S", "E", "I")) %>% apply(1, str_c, collapse="."),
                   "time")) %>%
    gather("details", "prevalence", -time) %>%
    separate(details, c("patch", "equation"), "\\.") %>%
    mutate(variable = factor(equation, levels=c("S", "E", "I"),
                             labels=c("Susceptible", "Exposed", "Infectious")),
           patch = factor(patch, levels=str_c("Patch_", 1:n_patches),
                               labels=str_c("Patch ", LETTERS[1:n_patches]))) %>%
    as_tibble() %>%
    mutate(network = "A <- B")
) %>%
  mutate(network = fct_inorder(network))

#R =
# sigma * beta * r / gamma / (sigma + r + mu) / (nu + nu)

ggplot(timeseries %>% filter(variable == "Infectious")) +
  aes(colour=patch, y=prevalence, x=time) +
  geom_line(size=0.75) +
  facet_grid(.~network, scales="free") +
  scale_x_continuous(expand=c(0,0), name="Time") +
  coord_trans(y="mysqrt") +
  scale_colour_manual(values=my_cols) +
  ylab("Prevalence") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="bottom",
        legend.background=element_blank(),
        legend.title=element_blank())
ggsave("two-pop-comparison.pdf", width=8, height=3)

## how many extrema (what are the dynamics for each patch in each network structure)
timeseries %>%
  group_by(equation, patch, network) %>%
  ## note that results are subject to rounding;
  ## this level corresponds to mimimum accuracy set in Julia integration
  summarise(n_extrema = prevalence %>% round(8) %>% unique() %>%
              {ifelse(tail(., 1) == 0, -1, (. > lead(., default=-1) & . > lag(., default=-1)) %>% sum())}) %>%
  pivot_wider(names_from=network, values_from=n_extrema) %>%
  arrange(equation, patch)
