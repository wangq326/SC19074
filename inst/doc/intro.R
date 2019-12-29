## -----------------------------------------------------------------------------
library(tidyverse)
library(rrpack)
library(knitr)
library(mvtnorm)
library(gridExtra)
library(MASS)
library(xtable)
set.seed(6419)

## -----------------------------------------------------------------------------
n = 38 # samples
q = 44 # number of response variables
rho = 0.5 # correlation between Xs
rX = 8 # rank of X
r = 5 # rank of C
p = 23 # number of explanatory variables

## -----------------------------------------------------------------------------
library(SC19074)
sim21 <- rrr.sim6(n, p, q, r, b = 0.05, rho, rX)
sim_ann21 <- rrr2(sim21$Y, sim21$X, penaltySVD = "ann",
                 modstr = list(gamma = 2))
sim_res21 <- stability_ann_sim(sim21$Y, sim21$X, sim_ann21$lambda, subset = 0.7, ntimes = 100)

## -----------------------------------------------------------------------------
sim22 <- rrr.sim6(n, p, q, r, b = 0.1, rho, rX)
sim_ann22 <- rrr2(sim22$Y, sim22$X, penaltySVD = "ann",
                 modstr = list(gamma = 2))
sim_res22 <- stability_ann_sim(sim22$Y, sim22$X, sim_ann22$lambda, subset = 0.7, ntimes = 100)

## -----------------------------------------------------------------------------
sim23 <- rrr.sim6(n, p, q, r, b = 0.3, rho, rX)
sim_ann23 <- rrr2(sim23$Y, sim23$X, penaltySVD = "ann",
                 modstr = list(gamma = 2))
sim_res23 <- stability_ann_sim(sim23$Y, sim23$X, sim_ann23$lambda, subset = 0.7, ntimes = 100)

## -----------------------------------------------------------------------------
# Rank 0
# Generate X
p = 46
sim4 <- rrr.sim6(n, p, q, r, b = 0, rho, rX)
sim_ann4 <- rrr2(sim4$Y, sim4$X, penaltySVD = "ann",
                 modstr = list(gamma = 2))
sim_res4 <- stability_ann_sim(sim4$Y, sim4$X, sim_ann4$lambda, subset = 0.7, ntimes = 100)

## -----------------------------------------------------------------------------
detach("package:rrpack", unload = TRUE)
detach("package:MASS", unload=TRUE)

# ---- stability-s2n-table ----
# Table of signal to noise ratios for each simulation
s2n_table <- data.frame(Simulation = 1:4, n = n, p = p, q = q, b = c(0, 0.05, 0.1, 0.3),
                       "Signal to noise ratio" = c(sim4$s2n, sim21$s2n, sim22$s2n, sim23$s2n),
                       check.names = FALSE)
print(xtable(s2n_table, label = "tab:stability-s2n-table",
       caption = "Parameters and signal to noise ratios for the four simulated datasets.",
       digits = c(0,0,0,0,0, 2, 3)),
      include.rownames = FALSE)

## -----------------------------------------------------------------------------
# ---- stability-sims-plots ----
stability_combine_sim <- function(stability_ann_results){
  st_plot <- data.frame(rank = c(stability_ann_results$st_res),
                        lambda = rep(stability_ann_results$lambda, each = nrow(stability_ann_results$st_res)))
  st_group <- group_by(st_plot, by = lambda)
  st_sum <- summarize(st_group, var = var(rank), avg = mean(rank)) %>% rename(lambda = by)
  st_sum <- full_join(st_sum, stability_ann_results$st_coef, by = "lambda")
}
stability_plots_sim <- function(sim_results, rank, title = "sim_results", subtitle = ""){
  st_combined <- stability_combine_sim(sim_results)
  # Plot mean vs variance
  varmeanplot <- ggplot(st_combined, aes(x = avg, y = sqrt(var))) + 
    geom_point() +
    geom_vline(aes(xintercept = rank)) +
    xlab("Mean estimated rank") + 
    ylab("Standard deviation of subspace distance") +
    ggtitle(title, subtitle = subtitle) +
    theme_bw()
  # Plot mean vs variance of coefficients
  coefplot <- ggplot(st_combined, aes(x = avg, y = sqrt(var.subspace.coef))) + 
    geom_point() +
    geom_vline(aes(xintercept = rank)) +
    xlab("Mean estimated rank") + 
    ylab(expression(sqrt("Mean of coefficient ranks"))) +
    ggtitle(title, subtitle = subtitle) +
    theme_bw()
  
  return(list(rplot = varmeanplot, coefplot = coefplot))
}

## -----------------------------------------------------------------------------
p21 <- stability_plots_sim(sim_res21, rank = r, title = "Simulation 2", subtitle = "b = 0.05")

p22 <- stability_plots_sim(sim_res22, rank = r, title = "Simulation 3", subtitle = "b = 0.1")

p23 <- stability_plots_sim(sim_res23, rank = r, title = "Simulation 4", subtitle = "b = 0.3")

p4 <- stability_plots_sim(sim_res4, rank = 0, title = "Simulation 1", subtitle = "b = 0")

## -----------------------------------------------------------------------------
# ---- stability-sims-rplots ----
grid.arrange(p4$rplot, p21$rplot, p22$rplot, p23$rplot, nrow = 2)

# ---- stability-sims-coefplots ----
grid.arrange(p4$coefplot, p21$coefplot, p22$coefplot, p23$coefplot, nrow = 2)

# ---- stability-sims-r0-plots ----
grid.arrange(p4$rplot, p4$coefplot, nrow = 1)

## -----------------------------------------------------------------------------
# ---- rank1-sims ----
# Simulations and plots for rank 1 matrices
set.seed(9988999)
n = 38 # samples
q = 44 # number of response variables
rho = 0.5 # correlation between Xs
rX = 8 # rank of X
p = 23 # number of explanatory variables
simr1 <- rrr.sim6(n, p, q, r=1, b = 0.3, rho, rX)
simr1$s2n
sim_annr1 <- rrr2(simr1$Y, simr1$X, penaltySVD = "ann",
                  modstr = list(gamma = 2))
sim_resr1 <- stability_ann_sim(simr1$Y, simr1$X, sim_annr1$lambda, subset = 0.7, ntimes = 100)

simr12 <- rrr.sim6(n, p, q, r=1, b = 0.05, rho, rX)
simr12$s2n
sim_annr12 <- rrr2(simr12$Y, simr12$X, penaltySVD = "ann",
                  modstr = list(gamma = 2))
sim_resr12 <- stability_ann_sim(simr12$Y, simr12$X, sim_annr12$lambda, subset = 0.7, ntimes = 100)

simr13 <- rrr.sim6(n, p, q, r=1, b = 0.1, rho, rX)
simr13$s2n
sim_annr13 <- rrr2(simr13$Y, simr13$X, penaltySVD = "ann",
                  modstr = list(gamma = 2))
sim_resr13 <- stability_ann_sim(simr13$Y, simr13$X, sim_annr13$lambda, subset = 0.7, ntimes = 100)

simr14 <- rrr.sim6(n, p, q, r=1, b = 0, rho, rX)
simr14$s2n
sim_annr14 <- rrr2(simr14$Y, simr14$X, penaltySVD = "ann",
                   modstr = list(gamma = 2))
sim_resr14 <- stability_ann_sim(simr14$Y, simr14$X, sim_annr14$lambda, subset = 0.7, ntimes = 100)

pr1 <- stability_plots_sim(sim_resr1, rank = 1, title = "", subtitle = paste("Signal-to-noise ratio: ", round(simr1$s2n, 2)))
pr12 <- stability_plots_sim(sim_resr12, rank = 1, title = "", subtitle = paste("Signal-to-noise ratio: ", round(simr12$s2n, 2)))
pr13 <- stability_plots_sim(sim_resr13, rank = 1, title = "", subtitle = paste("Signal-to-noise ratio: ", round(simr13$s2n, 2)))
pr14 <- stability_plots_sim(sim_resr14, rank = 1, title = "", subtitle = paste("Signal-to-noise ratio: ", round(simr14$s2n, 2)))

## -----------------------------------------------------------------------------
grid.arrange(pr1$coefplot, pr12$coefplot, pr13$coefplot, pr14$coefplot, nrow = 2)

