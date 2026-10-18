# vg_model.R
# Van Genuchten (1980) soil water retention model: fitting and derived
# soil hydraulic properties.
#
# Units used throughout:
#   h      suction head in cm (positive number)
#   pF     log10(h in cm)
#   theta  volumetric water content, m3/m3
#   alpha  1/cm
#   n      dimensionless (> 1); m = 1 - 1/n (Mualem constraint)

# ---------------------------------------------------------------------------
# Convert between pF and suction head
# ---------------------------------------------------------------------------
pf_to_h <- function(pF) {
  10^pF
}

h_to_pf <- function(h) {
  log10(h)
}

# ---------------------------------------------------------------------------
# Van Genuchten water content at suction head h
#   theta(h) = theta_r + (theta_s - theta_r) / (1 + (alpha*h)^n)^m
# ---------------------------------------------------------------------------
vg_theta <- function(h, theta_r, theta_s, alpha, n) {
  if (any(n <= 1)) {
    stop("n must be greater than 1")
  }
  if (any(theta_s <= theta_r)) {
    stop("theta_s must be greater than theta_r")
  }
  m <- 1 - 1 / n
  theta_r + (theta_s - theta_r) / (1 + (alpha * h)^n)^m
}

# ---------------------------------------------------------------------------
# Fit the van Genuchten model to measured (pF, theta) pairs for one sample.
#
# Uses bounded non-linear least squares (nls, algorithm = "port").
# Several starting values for alpha and n are tried, because nls is
# sensitive to where it starts. The fit with the lowest residual sum
# of squares is kept.
#
# Returns a list with parameters, RMSE, R2 and the nls object.
# ---------------------------------------------------------------------------
fit_vg <- function(pF, theta) {
  if (length(pF) != length(theta)) {
    stop("pF and theta must have the same length")
  }
  if (length(pF) < 5) {
    stop("at least 5 measurement points are needed to fit 4 parameters")
  }

  h <- pf_to_h(pF)
  dat <- data.frame(h = h, theta = theta)

  # Starting guesses: theta_s near the wettest measurement, theta_r near zero
  theta_s_start <- max(theta)
  theta_r_start <- 0.5 * min(theta)

  alpha_starts <- c(0.001, 0.01, 0.05)
  n_starts <- c(1.1, 1.5, 2.5)

  best_fit <- NULL
  best_rss <- -1

  for (a0 in alpha_starts) {
    for (n0 in n_starts) {
      fit <- tryCatch(
        nls(
          theta ~ vg_theta(h, theta_r, theta_s, alpha, n),
          data = dat,
          start = list(theta_r = theta_r_start, theta_s = theta_s_start,
                       alpha = a0, n = n0),
          lower = c(theta_r = 0, theta_s = 0.2, alpha = 1e-5, n = 1.01),
          upper = c(theta_r = min(theta), theta_s = 0.8, alpha = 1, n = 8),
          algorithm = "port"
        ),
        error = function(e) NULL
      )

      if (!is.null(fit)) {
        rss <- sum(residuals(fit)^2)
        if (best_rss < 0 || rss < best_rss) {
          best_fit <- fit
          best_rss <- rss
        }
      }
    }
  }

  if (is.null(best_fit)) {
    stop("van Genuchten fit did not converge for any starting values")
  }

  p <- coef(best_fit)
  rmse <- sqrt(best_rss / length(theta))
  r2 <- 1 - best_rss / sum((theta - mean(theta))^2)

  list(
    theta_r = unname(p["theta_r"]),
    theta_s = unname(p["theta_s"]),
    alpha = unname(p["alpha"]),
    n = unname(p["n"]),
    rmse = rmse,
    r2 = r2,
    model = best_fit
  )
}

# ---------------------------------------------------------------------------
# Derived hydraulic properties from fitted parameters.
#
#   Field capacity (FC):      theta at pF 2.0  (h = 100 cm, ~ -10 kPa)
#   Permanent wilting point:  theta at pF 4.2  (h ~ 15 850 cm, ~ -1500 kPa)
#   Plant-available water:    FC - PWP
#
# pF 2.0 is the convention used in Finnish and Nordic soil surveys;
# some studies use pF 2.5 (-33 kPa) instead, so it is an argument.
# ---------------------------------------------------------------------------
hydraulic_properties <- function(params, pf_fc = 2.0, pf_pwp = 4.2) {
  fc <- vg_theta(pf_to_h(pf_fc), params$theta_r, params$theta_s,
                 params$alpha, params$n)
  pwp <- vg_theta(pf_to_h(pf_pwp), params$theta_r, params$theta_s,
                  params$alpha, params$n)
  paw <- fc - pwp

  # Plant-available water in mm per 1 m soil depth (m3/m3 * 1000 mm)
  paw_mm_per_m <- paw * 1000

  # Air capacity: pore space drained at field capacity
  air_capacity <- params$theta_s - fc

  data.frame(
    field_capacity = fc,
    wilting_point = pwp,
    plant_available_water = paw,
    paw_mm_per_m = paw_mm_per_m,
    air_capacity = air_capacity
  )
}

# ---------------------------------------------------------------------------
# Bootstrap uncertainty for plant-available water.
#
# Residual bootstrap: refit the model to (fitted value + resampled
# residual) many times and take the 2.5 and 97.5 percentiles.
# With only 7-8 points per sample this is rough, but it shows how
# much the estimate can move, which a single number hides.
# ---------------------------------------------------------------------------
bootstrap_paw <- function(pF, theta, n_boot = 200, seed = 1) {
  set.seed(seed)
  base_fit <- fit_vg(pF, theta)
  fitted_vals <- fitted(base_fit$model)
  res <- residuals(base_fit$model)

  paw_values <- c()
  for (i in seq_len(n_boot)) {
    theta_star <- fitted_vals + sample(res, length(res), replace = TRUE)
    fit_star <- tryCatch(fit_vg(pF, theta_star), error = function(e) NULL)
    if (!is.null(fit_star)) {
      paw_star <- hydraulic_properties(fit_star)$plant_available_water
      paw_values <- c(paw_values, paw_star)
    }
  }

  if (length(paw_values) < 0.5 * n_boot) {
    warning("more than half of bootstrap refits failed; interval unreliable")
  }

  c(
    lower = unname(quantile(paw_values, 0.025)),
    upper = unname(quantile(paw_values, 0.975))
  )
}
