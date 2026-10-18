# make_example_data.R
# Generates the ILLUSTRATIVE example dataset in data/example_retention.csv.
#
# These are NOT field measurements. Each texture class uses the
# class-average van Genuchten parameters from the ROSETTA model
# (Schaap, Leij & van Genuchten 2001, J. Hydrology 251:163-176),
# with small random measurement noise added. Clay and silty clay loam
# are included because they dominate agricultural fields in southern
# Finland.
#
# The pF steps match a typical sandbox + pressure plate laboratory
# series. Replace this file with your own measurements to analyse
# real soils (same columns).

source("R/vg_model.R")

set.seed(42)

classes <- data.frame(
  texture = c("clay", "silty_clay_loam", "loam", "sandy_loam"),
  theta_r = c(0.098, 0.090, 0.061, 0.039),
  theta_s = c(0.459, 0.482, 0.399, 0.387),
  alpha = 10^c(-1.825, -2.076, -1.954, -1.574),
  n = 10^c(0.098, 0.182, 0.168, 0.161)
)

pf_steps <- c(0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.7, 4.2)
replicates <- 3
noise_sd <- 0.008

rows <- list()
for (i in seq_len(nrow(classes))) {
  cl <- classes[i, ]
  for (r in seq_len(replicates)) {
    theta_true <- vg_theta(pf_to_h(pf_steps), cl$theta_r, cl$theta_s,
                           cl$alpha, cl$n)
    theta_obs <- theta_true + rnorm(length(pf_steps), 0, noise_sd)
    rows[[length(rows) + 1]] <- data.frame(
      sample_id = sprintf("%s_%d", cl$texture, r),
      texture = cl$texture,
      pF = pf_steps,
      theta = round(theta_obs, 4)
    )
  }
}

out <- do.call(rbind, rows)
write.csv(out, "data/example_retention.csv", row.names = FALSE)
cat("Wrote", nrow(out), "rows to data/example_retention.csv\n")
