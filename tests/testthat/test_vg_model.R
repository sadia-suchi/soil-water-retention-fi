# Tests for R/vg_model.R
# Run from the repository root:  Rscript -e 'testthat::test_dir("tests/testthat")'

source(file.path("..", "..", "R", "vg_model.R"), chdir = TRUE)

test_that("water content at zero suction equals theta_s", {
  expect_equal(vg_theta(0, 0.1, 0.45, 0.02, 1.5), 0.45)
})

test_that("water content decreases as suction increases", {
  h <- pf_to_h(c(0, 1, 2, 3, 4.2))
  theta <- vg_theta(h, 0.1, 0.45, 0.02, 1.5)
  expect_true(all(diff(theta) < 0))
})

test_that("water content stays between theta_r and theta_s", {
  h <- pf_to_h(seq(0, 7, by = 0.5))
  theta <- vg_theta(h, 0.08, 0.40, 0.03, 1.4)
  expect_true(all(theta <= 0.40))
  expect_true(all(theta >= 0.08))
})

test_that("invalid parameters are rejected", {
  expect_error(vg_theta(100, 0.1, 0.45, 0.02, 0.9), "n must be greater than 1")
  expect_error(vg_theta(100, 0.5, 0.45, 0.02, 1.5), "theta_s must be greater")
})

test_that("fit recovers known parameters from noise-free data", {
  pF <- c(0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.7, 4.2)
  true <- list(theta_r = 0.09, theta_s = 0.48, alpha = 0.0084, n = 1.52)
  theta <- vg_theta(pf_to_h(pF), true$theta_r, true$theta_s, true$alpha, true$n)

  fit <- fit_vg(pF, theta)

  expect_equal(fit$theta_s, true$theta_s, tolerance = 0.01)
  expect_equal(fit$n, true$n, tolerance = 0.05)
  expect_lt(fit$rmse, 1e-4)
})

test_that("fit needs at least 5 points", {
  expect_error(fit_vg(c(0, 2, 4.2), c(0.45, 0.3, 0.15)), "at least 5")
})

test_that("plant-available water is field capacity minus wilting point", {
  params <- list(theta_r = 0.06, theta_s = 0.40, alpha = 0.011, n = 1.47)
  props <- hydraulic_properties(params)
  expect_equal(props$plant_available_water,
               props$field_capacity - props$wilting_point)
  expect_gt(props$plant_available_water, 0)
  expect_equal(props$paw_mm_per_m, props$plant_available_water * 1000)
})
