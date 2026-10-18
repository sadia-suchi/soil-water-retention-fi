# run_analysis.R
# Fits the van Genuchten model to every sample in a retention dataset,
# calculates hydraulic properties, and writes tables and figures.
#
# Usage (from the repository root):
#   Rscript scripts/run_analysis.R                      # example data
#   Rscript scripts/run_analysis.R path/to/your.csv     # your own data
#
# Input CSV columns: sample_id, texture, pF, theta

suppressPackageStartupMessages(library(ggplot2))
source("R/vg_model.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1) {
  input_file <- args[1]
} else {
  input_file <- "data/example_retention.csv"
}

dat <- read.csv(input_file)

required_cols <- c("sample_id", "texture", "pF", "theta")
missing_cols <- setdiff(required_cols, names(dat))
if (length(missing_cols) > 0) {
  stop("input is missing columns: ", paste(missing_cols, collapse = ", "))
}

dir.create("results", showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Fit each sample
# ---------------------------------------------------------------------------
sample_ids <- unique(dat$sample_id)
summary_rows <- list()
curve_rows <- list()

for (sid in sample_ids) {
  s <- dat[dat$sample_id == sid, ]
  fit <- fit_vg(s$pF, s$theta)
  props <- hydraulic_properties(fit)
  ci <- bootstrap_paw(s$pF, s$theta, n_boot = 200)

  summary_rows[[sid]] <- data.frame(
    sample_id = sid,
    texture = s$texture[1],
    theta_r = fit$theta_r,
    theta_s = fit$theta_s,
    alpha = fit$alpha,
    n = fit$n,
    rmse = fit$rmse,
    r2 = fit$r2,
    props,
    paw_ci_lower = ci["lower"],
    paw_ci_upper = ci["upper"]
  )

  # Smooth fitted curve for plotting
  pf_grid <- seq(0, 4.5, by = 0.05)
  curve_rows[[sid]] <- data.frame(
    sample_id = sid,
    texture = s$texture[1],
    pF = pf_grid,
    theta = vg_theta(pf_to_h(pf_grid), fit$theta_r, fit$theta_s,
                     fit$alpha, fit$n)
  )
}

summary_tbl <- do.call(rbind, summary_rows)
curves <- do.call(rbind, curve_rows)
rownames(summary_tbl) <- NULL

write.csv(summary_tbl, "results/fit_summary.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 2. Texture-level means (what you would report in a paper table)
# ---------------------------------------------------------------------------
texture_means <- aggregate(
  cbind(field_capacity, wilting_point, plant_available_water,
        paw_mm_per_m, air_capacity, rmse) ~ texture,
  data = summary_tbl,
  FUN = mean
)
write.csv(texture_means, "results/texture_means.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Does plant-available water differ between textures?
#    One-way ANOVA + Tukey HSD (same approach as in my MSc lab work).
# ---------------------------------------------------------------------------
aov_fit <- aov(paw_mm_per_m ~ texture, data = summary_tbl)
tukey <- TukeyHSD(aov_fit)
sink("results/anova_paw.txt")
cat("One-way ANOVA: plant-available water (mm per m) by texture\n\n")
print(summary(aov_fit))
cat("\nTukey HSD\n\n")
print(tukey)
sink()

# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------
p1 <- ggplot() +
  geom_line(data = curves,
            aes(x = theta, y = pF, group = sample_id, colour = texture),
            linewidth = 0.6, alpha = 0.8) +
  geom_point(data = dat,
             aes(x = theta, y = pF, colour = texture),
             size = 1.6) +
  geom_hline(yintercept = c(2.0, 4.2), linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.02, y = 2.0, label = "field capacity (pF 2.0)",
           hjust = 0, vjust = -0.4, size = 3, colour = "grey30") +
  annotate("text", x = 0.02, y = 4.2, label = "wilting point (pF 4.2)",
           hjust = 0, vjust = -0.4, size = 3, colour = "grey30") +
  labs(x = expression("Volumetric water content (m"^3*" m"^-3*")"),
       y = "pF  (log10 suction, cm)",
       colour = "Texture",
       title = "Soil water retention curves, van Genuchten fits") +
  theme_minimal(base_size = 11)
ggsave("results/retention_curves.png", p1, width = 7, height = 5, dpi = 150)

p2 <- ggplot(summary_tbl,
             aes(x = reorder(texture, paw_mm_per_m), y = paw_mm_per_m)) +
  geom_errorbar(aes(ymin = paw_ci_lower * 1000, ymax = paw_ci_upper * 1000),
                width = 0.15, colour = "grey50",
                position = position_jitter(width = 0.12, seed = 1)) +
  geom_point(size = 2.4, colour = "#2b6c8f",
             position = position_jitter(width = 0.12, seed = 1)) +
  labs(x = NULL, y = "Plant-available water (mm per m soil)",
       title = "Plant-available water by texture",
       subtitle = "Points = samples; bars = 95% bootstrap interval") +
  theme_minimal(base_size = 11)
ggsave("results/paw_by_texture.png", p2, width = 6, height = 4, dpi = 150)

cat("Done. Fitted", length(sample_ids), "samples. Results in results/\n")
print(texture_means, digits = 3)
