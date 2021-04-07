library(magrittr)
library(tidyverse)
library(JuliaCall)
library(rstatix)
library(scales)

julia_setup()
julia_source("../../Code/multipop_SI.jl")

my_cols <- c("#f74700", "#016394", "#b6003b", "#005342", "#eac435")

mysqrt_trans <- function() {
  trans_new("mysqrt",
            transform = base::sqrt,
            inverse = function(x) ifelse(x<0, 0, x^2),
            domain = c(0, Inf))
}

## structural parameters
set.seed(0)
n_patches <- 4
# initial_conditions <- rep(0.00001, prod(struct) * n_patches)
initial_conditions <- matrix(c(5, 1, 1), n_patches, 3, byrow=TRUE) #+ runif(n_patches * 3, -0.5, 2)

## dynamical parameters
r <- 0.5             # intrinsic growth rate
K <- 15#runif(n_patches, 13, 17) # carrying capacity
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
# maxtime <- 15000
# tsteps <- round(0.3334 * maxtime):maxtime

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

ggplot(timeseries %>% filter(variable == "Infectious")) + # %>% filter(patch != "Patch A", patch != "Patch J")) +
  aes(colour=patch, y=prevalence, x=time) +
  geom_line(size=0.75) +
  facet_grid(.~patch, scales="free") +
  scale_colour_manual(values=my_cols) +
  scale_x_continuous(expand=c(0,0), name="Time") +
  coord_trans(y="mysqrt") +
  ylab("Prevalence") +
  theme_bw() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position="none")
ggsave("dampening-cycles.pdf", width=8, height=2.5)

# timeseries %>%
#   filter(patch != "Patch A") %>% #, patch != "Z") %>%
#   group_by(patch, equation) %>%
#   summarise(across(prevalence, list("mean"=mean, "max"=max, "min"=min)), .groups="drop") %>%
#   ggplot() +
#   aes(x=as.integer(patch),
#       y=prevalence_mean, ymax=prevalence_max, ymin=prevalence_min,
#       # label=str_remove(patch, "Patch "),
#       colour=patch) +
#   geom_linerange(size=1) + geom_point(size=2) + #geom_smooth(method="lm", se=FALSE) +
#   facet_grid(equation~., scales="free") +
#   theme_bw()

# timeseries %>%
#   filter(patch != "Patch A") %>% #, patch != "Z") %>%
#   group_by(patch, equation) %>%
#   summarise(across(prevalence, list("mean"=mean, "max"=max, "min"=min)), .groups="drop") %>%
#   ggplot() +
#   aes(x=as.integer(patch),
#       y=prevalence_max - prevalence_min,
#       # label=str_remove(patch, "Patch "),
#       colour=patch) +
#   geom_point(size=2) + #geom_smooth(method="lm", colour="black", se=FALSE) +
#   facet_grid(equation~., scales="free") +
#   theme_bw()

# cor_mat <- timeseries %>%
#   filter(equation == "I") %>%
#   select(time, patch, prevalence) %>%
#   pivot_wider(names_from="patch", values_from=prevalence) %>%
#   select(-time) %>%
#   cor_mat(method="spearman") %>%
#   rename_with(~str_remove(., "Patch "))
# pdf(str_glue("long_chain_corrplot_{movement_rate}.pdf"), width=8, height=8);
# cor_plot(cor_mat, method="shade", type="upper",
#          palette=ggpubr::get_palette(c(my_cols[3], "white", my_cols[2]), 200));
# dev.off()
# save.image(str_glue("long_chain_corrplot_{movement_rate}.RData"))
