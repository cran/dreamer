test_that("get_extreme.quad()", {
  data <- dreamer_data_linear(n_cohorts = c(10, 20, 30), c(1, 3, 5), 1, 2, 2)
  mcmc <- dreamer_mcmc(
    data,
    mod = model_quad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01
    ),
    n_iter = 4,
    n_chains = 1,
    silent = TRUE,
    convergence_warn = FALSE
  )

  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))
  b1 <- 1:4
  b2 <- c(- 1.25, - 1.5, 2, 2.5)
  b3 <- c(.5, .5, - .5, - .5)
  mcmc <- mcmc %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3)
  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = TRUE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(5, 5, 2, 2.5)) %>%
    dplyr::mutate(
      extreme_responses = b1 + b2 * doses + b3 * doses ^ 2,
      greater = TRUE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(1.25, 1.5, 5, 5)) %>%
    dplyr::mutate(
      extreme_responses = b1 + b2 * doses + b3 * doses ^ 2,
      greater = FALSE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = NULL,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = 2
  )
  exp <- tibble::tibble(doses = c(1.25, 1.5, 5, 5)) %>%
    dplyr::mutate(
      extreme_responses = b1 + b2 * doses + b3 * doses ^ 2,
      greater = FALSE
    ) %>%
    dplyr::slice(2)
  expect_equal(obs, exp)
})

test_that("get_extreme.quad() longitudinal", {
  times <- c(0, 10)
  t_max <- max(times)
  data <- dreamer_data_linear(
    n_cohorts = c(10, 25, 30),
    dose = c(0, 2, 4),
    b1 = .5,
    b2 = 3,
    sigma = .5,
    longitudinal = "linear",
    a = .5,
    times = times,
    t_max = t_max
  )
  mcmc <- dreamer_mcmc(
    data = data,
    n_iter = 4,
    n_chains = 1,
    convergence_warn = FALSE,
    silent = TRUE,
    mod = model_quad(
      mu_b1 = 0,
      sigma_b1 = 1,
      mu_b2 = 0,
      sigma_b2 = 1,
      mu_b3 = 0,
      sigma_b3 = 1,
      shape = 1,
      rate = .01,
      longitudinal = model_longitudinal_linear(0, 1, t_max)
    )
  )

  lower <- min(attr(mcmc, "doses"))
  upper <- max(attr(mcmc, "doses"))

  a <- c(.1, .2, .3, .4)
  b1 <- 1:4
  b2 <- c(- 1.25, - 1.5, 2, 2.5)
  b3 <- c(.5, .5, - .5, - .5)

  mcmc <- mcmc %>%
    replace_mcmc("mod", "a", a) %>%
    replace_mcmc("mod", "b1", b1) %>%
    replace_mcmc("mod", "b2", b2) %>%
    replace_mcmc("mod", "b3", b3)

  time <- 3
  obs <- get_extreme(
    mcmc$mod,
    time = time,
    greater = TRUE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(4, 4, 2, 2.5)) %>%
    dplyr::mutate(
      extreme_responses = a + (b1 + b2 * doses + b3 * doses ^ 2) * time / t_max,
      greater = TRUE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = time,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = NULL
  )
  exp <- tibble::tibble(doses = c(1.25, 1.5, 0, 0)) %>%
    dplyr::mutate(
      extreme_responses = a + (b1 + b2 * doses + b3 * doses ^ 2) * time / t_max,
      greater = FALSE
    )
  expect_equal(obs, exp)

  obs <- get_extreme(
    mcmc$mod,
    time = time,
    greater = FALSE,
    lower = lower,
    upper = upper,
    index = 2
  )
  exp <- tibble::tibble(doses = c(1.25, 1.5, 0, 0)) %>%
    dplyr::mutate(
      extreme_responses = a + (b1 + b2 * doses + b3 * doses ^ 2) * time / t_max,
      greater = FALSE
    ) %>%
    dplyr::slice(2)
  expect_equal(obs, exp)
})
