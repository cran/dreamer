## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(fig.width = 6, fig.height = 4)

## ---- message = FALSE---------------------------------------------------------
library(dreamer)
library(dplyr)
library(ggplot2)
set.seed(888)
data <- dreamer_data_quad(
  n_cohorts = c(20, 20, 20, 20, 20), # number of subjects in each cohort
  doses = c(0, .5, 1, 1.5, 2), # dose administered to each cohort
  b1 = 5, # intercept
  b2 = 5,
  b3 = - 2,
  sigma = 2 # standard deviation
)
head(data)
ggplot(data, aes(dose, response)) + geom_point()

## -----------------------------------------------------------------------------
model_quad(
  mu_b1 = 0,
  sigma_b1 = 10,
  mu_b2 = 0,
  sigma_b2 = 10,
  mu_b3 = 0,
  sigma_b3 = 10,
  shape = 1,
  rate = .001,
  w_prior = 1 / 2
)

## -----------------------------------------------------------------------------
# Bayesian model averaging
output <- dreamer_mcmc(
  data = data,
  n_adapt = 1e3,
  n_burn = 1e3,
  n_iter = 1e4,
  n_chains = 2,
  silent = TRUE,
  # the argument name "mod_quad" is arbitrary and is chosen by the user
  mod_quad = model_quad(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 10,
    shape = 1,
    rate = .001,
    w_prior = 1 / 2
  ),
  # the argument name "mod_emax" is arbitrary and is chosen by the user
  mod_emax = model_emax(
    mu_b1 = 0,
    sigma_b1 = 10,
    mu_b2 = 0,
    sigma_b2 = 10,
    mu_b3 = 0,
    sigma_b3 = 2,
    mu_b4 = 1,
    sigma_b4 = 5,
    shape = 1,
    rate = .001,
    w_prior = 1 / 2
  )
)

## -----------------------------------------------------------------------------
output

## -----------------------------------------------------------------------------
summary(output)

## -----------------------------------------------------------------------------
posterior(output)

# posterior Pr(effect of dose - effect of reference dose)
posterior(x = output, reference_dose = 0)

## -----------------------------------------------------------------------------
plot(output)

## -----------------------------------------------------------------------------
plot(output, data = data)

## -----------------------------------------------------------------------------
plot(output, reference_dose = 0) # adjust relative to dose = 0

## -----------------------------------------------------------------------------
plot(output, data = data, predictive = 10)

## -----------------------------------------------------------------------------
plot(output, reference_dose = 0, predictive = 10)

## -----------------------------------------------------------------------------
plot(output$mod_emax, data = data)

## -----------------------------------------------------------------------------
plot_comparison(output)

## -----------------------------------------------------------------------------
plot_comparison(
  mod_emax = output$mod_emax,
  mod_quad = output$mod_quad
)

## -----------------------------------------------------------------------------
# absolute: pr(mean at dose 0.50 > 1 | data)
pr_eoi(output, eoi = 1, dose = 0.50)

# relative: pr( (mean at dose 0.50) - (mean at dose 0.25) > 0.55 | data )
pr_eoi(output, eoi = 0.55, dose = 0.50, reference_dose = 0.25)

# vectorized
n_doses <- length(output$doses)
pr_eoi(
  output,
  eoi = rep(.55, n_doses),
  dose = output$doses,
  reference_dose = rep(0, n_doses)
)

## -----------------------------------------------------------------------------
# from the quadratic model
pr_eoi(
  output$mod_quad,
  eoi = 1.5,
  dose = .75,
  reference_dose = 0
)

## -----------------------------------------------------------------------------
pr_med(output, csd = 8)

# placebo adjusted
pr_med(
  output,
  csd = 3,
  reference_dose = 0
)

# relative to placebo grp (dose = 0) for just the EMAX model
pr_med(
  output$mod_emax,
  csd = 3,
  reference_dose = 0
)

## -----------------------------------------------------------------------------
# looking for smallest dose with 95% of the maximum efficacy
pr_medx(output, ed = 95)

## -----------------------------------------------------------------------------
post_medx(output, ed = 95)

## -----------------------------------------------------------------------------
post_perc_effect(
  output,
  dose = c(.05, .5)
)

## -----------------------------------------------------------------------------
diagnostics(output)
# single model
diagnostics(output$mod_emax)

## ---- fig.height = 8, fig.width = 6-------------------------------------------
# single model
plot_trace(output$mod_quad)

# traceplot for all parameters for each model using: plot_trace(output)

## -----------------------------------------------------------------------------
dreamer_plot_prior(
  doses = c(0, 2.5, 5),
  mod_linear_binary = model_linear_binary(
    mu_b1 = - 1,
    sigma_b1 = .1,
    mu_b2 = 1,
    sigma_b2 = .1,
    link = "logit",
    w_prior = 1
  )
)

## -----------------------------------------------------------------------------
dreamer_plot_prior(
  doses = seq(from = 0, to = 5, length.out = 50),
  n_samples = 100,
  plot_draws = TRUE,
  mod_quad_binary = model_quad_binary(
    mu_b1 = - .5,
    sigma_b1 = .2,
    mu_b2 = - .5,
    sigma_b2 = .2,
    mu_b3 = .5,
    sigma_b3 = .1,
    link = "logit",
    w_prior = 1
  )
)

## -----------------------------------------------------------------------------
dreamer_plot_prior(
  doses = c(0, 2.5, 5),
  mod_linear_binary = model_linear_binary(
    mu_b1 = - 1,
    sigma_b1 = .1,
    mu_b2 = 1,
    sigma_b2 = .1,
    link = "logit",
    w_prior = .75
  ),
  mod_quad_binary = model_quad_binary(
    mu_b1 = - .5,
    sigma_b1 = .2,
    mu_b2 = - .5,
    sigma_b2 = .2,
    mu_b3 = .5,
    sigma_b3 = .1,
    link = "logit",
    w_prior = .25
  )
)

## -----------------------------------------------------------------------------
output_independent <- dreamer_mcmc(
  data = data,
  n_adapt = 1e3,
  n_burn = 1e3,
  n_iter = 1e4,
  n_chains = 2,
  silent = TRUE, # make rjags be quiet,
  # this model has the same prior on the mean for each dose
  mod_indep = model_independent(
    mu_b1 = 0,
    sigma_b1 = 1,
    shape = 1,
    rate = .001
  )
)
# prior is too strong!
plot(output_independent, data = data)

output_independent2 <- dreamer_mcmc(
  data = data,
  n_adapt = 1e3,
  n_burn = 1e3,
  n_iter = 1e4,
  n_chains = 2,
  silent = TRUE, # make rjags be quiet,
  # this model has the different priors on the mean for each dose
  mod_indep = model_independent(
    mu_b1 = c(0, 1, 2, 3, 4),
    sigma_b1 = c(10, 10, 20, 20, 30),
    shape = 1,
    rate = .001,
    doses = c(0, 0.5, 1, 1.5, 2)
  )
)
plot(output_independent2, data = data)

## -----------------------------------------------------------------------------
set.seed(889)
data_long <- dreamer_data_linear(
  n_cohorts = c(10, 10, 10, 10), # number of subjects in each cohort
  doses = c(.25, .5, .75, 1.5), # dose administered to each cohort
  b1 = 0, # intercept
  b2 = 2, # slope
  sigma = .5, # standard deviation,
  longitudinal = "itp",
  times = c(0, 12, 24, 52),
  t_max = 52, # maximum time
  a = .5,
  c1 = .1
)

ggplot(data_long, aes(time, response, group = dose, color = factor(dose))) +
  geom_point()

## -----------------------------------------------------------------------------
# Bayesian model averaging
output_long <- dreamer_mcmc(
  data = data_long,
  n_adapt = 1e3,
  n_burn = 1e3,
  n_iter = 1e4,
  n_chains = 2,
  silent = TRUE, # make rjags be quiet,
  mod_linear = model_linear(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    shape = 1,
    rate = .001,
    w_prior = 1 / 2, # prior probability of the model
    longitudinal = model_longitudinal_itp(
      mu_a = 0,
      sigma_a = 1,
      a_c1 = 0,
      b_c1 = 1,
      t_max = 52
    )
  ),
  mod_quad = model_quad(
    mu_b1 = 0,
    sigma_b1 = 1,
    mu_b2 = 0,
    sigma_b2 = 1,
    mu_b3 = 0,
    sigma_b3 = 1,
    shape = 1,
    rate = .001,
    w_prior = 1 / 2,
    longitudinal = model_longitudinal_linear(
      mu_a = 0,
      sigma_a = 1,
      t_max = 52
    )
  )
)

## -----------------------------------------------------------------------------
plot(output_long, data = data_long)
# plot dose response at last time point
plot(output_long, times = 52, data = data_long)

## -----------------------------------------------------------------------------
posterior(output_long)

