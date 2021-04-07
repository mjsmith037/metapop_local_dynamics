library(magrittr)
library(tidyverse)
library(JuliaCall)
library(scales)

julia_setup()
julia_source("../../code/multipop_SEI.jl")

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342", "#eac435")

mysqrt_trans <- function() {
  trans_new("mysqrt",
            transform = base::sqrt,
            inverse = function(x) ifelse(x<0, 0, x^2),
            domain = c(0, Inf))
}

## structural parameters
set.seed(0)
n_patches <- 3
# initial_conditions <- rep(0.00001, prod(struct) * n_patches)
initial_conditions <- matrix(c(5, 1, 1), n_patches, 3, byrow=TRUE)

## dynamical parameters
r <- 0.5             # intrinsic growth rate
K <- c(5, 15, 5)       # carrying capacity
beta <- 80           # infection rate
sigma <- 13          # 1 / sigma = latency period
mu <- 0.5            # intrinsic death rate
nu <- 73             # disease induced mortality
movement_rate <- 0.1
chi <- matrix(c(-movement_rate, 0, 0,
                0, -movement_rate, 0,
                movement_rate, movement_rate, 0), 3, 3)

## integration parameters
maxtime <- 10000
tsteps <- round(0.9 * maxtime):maxtime

## run the simulation
timeseries <- julia_call("runsimulation", tstep=tsteps, tmax=maxtime,
                         r=r, K=K, beta=beta, sigma=sigma, mu=mu, nu=nu, chi=chi,
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
  as_tibble()

ggplot(timeseries %>% filter(variable == "Infectious")) +
  aes(colour=patch, y=prevalence, x=time) +
  geom_line(size=0.75) +
  facet_grid(.~patch, scales="free") +
  scale_x_continuous(expand=c(0,0), name="Time") +
  coord_trans(y="mysqrt") +
  scale_colour_manual(values=my_cols) +
  ylab("Prevalence") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")
ggsave("dynamics-hierarchy.pdf", width=8, height=2.5)
