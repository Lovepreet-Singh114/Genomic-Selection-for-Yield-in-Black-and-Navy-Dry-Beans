library(tidyverse)
library(metan)
library(ggrepel)

# =============================================================================
# READ DATA
# =============================================================================
blues_env  <- read.csv("blues_environment.csv")
blues_year <- read.csv("blues_yearly.csv")
blues_loc  <- read.csv("blues_location.csv")

rename_for_gge <- function(df) {
  df %>%
    rename(GEN = Name, YIELD = BLUE, ENV = Group) %>%
    select(GEN, ENV, YIELD) %>%
    mutate(GEN = as.factor(GEN), ENV = as.factor(ENV))
}

blues_env  <- rename_for_gge(blues_env)
blues_year <- rename_for_gge(blues_year)
blues_loc  <- rename_for_gge(blues_loc)

# =============================================================================
# FIT GGE MODELS — Yan et al. (2007)
# =============================================================================
gge_env  <- gge(blues_env,  ENV, GEN, YIELD, centering = "environment", scaling = 1, svp = "environment")
gge_year <- gge(blues_year, ENV, GEN, YIELD, centering = "environment", scaling = FALSE, svp = "environment")
gge_loc  <- gge(blues_loc,  ENV, GEN, YIELD, centering = "environment", scaling = 1, svp = "environment")

# =============================================================================
# EXTRACT SCORES
# =============================================================================
extract_scores <- function(gge_obj) {
  obj <- gge_obj$YIELD
  
  env_df <- as.data.frame(obj$coordenv[, 1:2])
  colnames(env_df) <- c("PC1", "PC2")
  env_df$Label <- obj$labelenv
  
  gen_df <- as.data.frame(obj$coordgen[, 1:2])
  colnames(gen_df) <- c("PC1", "PC2")
  gen_df$Label <- obj$labelgen
  
  varexpl <- obj$varexpl[1:2]
  
  list(env = env_df, gen = gen_df, varexpl = varexpl)
}

scores_env  <- extract_scores(gge_env)
scores_year <- extract_scores(gge_year)
scores_loc  <- extract_scores(gge_loc)

# =============================================================================
# MEGA-ENVIRONMENT IDENTIFICATION
# Environments close to AEC line = within mega-environment
# Outliers = angle to AEC > mean + 1 SD
# =============================================================================
assign_mega_env <- function(env_df) {
  aec_x <- mean(env_df$PC1)
  aec_y <- mean(env_df$PC2)
  aec_len <- sqrt(aec_x^2 + aec_y^2)
  
  env_df <- env_df %>%
    mutate(
      vec_length = sqrt(PC1^2 + PC2^2),
      # Cosine of angle between env vector and AEC
      cos_angle = (PC1 * aec_x + PC2 * aec_y) / (vec_length * aec_len),
      # Angle in degrees
      angle_to_AEC = acos(pmin(pmax(cos_angle, -1), 1)) * 180 / pi
    )
  
  # Threshold: mean + 1 SD
  threshold <- mean(env_df$angle_to_AEC) + sd(env_df$angle_to_AEC)
  
  env_df <- env_df %>%
    mutate(
      Status = ifelse(angle_to_AEC <= threshold, "Mega-Environment", "Outlier")
    )
  
  cat("AEC angle threshold:", round(threshold, 2), "degrees\n")
  
  env_df
}

scores_env$env  <- assign_mega_env(scores_env$env)
scores_year$env <- assign_mega_env(scores_year$env)
scores_loc$env  <- assign_mega_env(scores_loc$env)

cat("\n=== Environment Classification ===\n")
cat("\nEnvironment-wise:\n")
print(scores_env$env %>% select(Label, angle_to_AEC, Status) %>% arrange(angle_to_AEC))
cat("\nYear-wise:\n")
print(scores_year$env %>% select(Label, angle_to_AEC, Status) %>% arrange(angle_to_AEC))
cat("\nLocation-wise:\n")
print(scores_loc$env %>% select(Label, angle_to_AEC, Status) %>% arrange(angle_to_AEC))

# =============================================================================
# COMMON STYLING
# =============================================================================
col_mega    <- "#2E8B57"
col_outlier <- "#D62728"
status_colors <- c("Mega-Environment" = col_mega, "Outlier" = col_outlier)

yan_theme <- theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0, size = 15, color = "grey40",
                                 margin = margin(b = 10)),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(face = "bold", size = 13),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(face = "bold", size = 12)
  )

axis_labels <- function(vx) {
  list(
    x = paste0("PC1 (", round(vx[1], 2), "%)"),
    y = paste0("PC2 (", round(vx[2], 2), "%)")
  )
}

subtitle_text <- "Scaling = 0, Centering = 2, SVP = 2"

get_limits <- function(env, gen, pad = 0.15) {
  all_x <- c(env$PC1, gen$PC1)
  all_y <- c(env$PC2, gen$PC2)
  rx <- diff(range(all_x))
  ry <- diff(range(all_y))
  list(
    xlim = c(min(all_x) - pad * rx, max(all_x) + pad * rx),
    ylim = c(min(all_y) - pad * ry, max(all_y) + pad * ry)
  )
}

# =============================================================================
# BIPLOT 1: ENVIRONMENT VECTORS
# =============================================================================
biplot_env_vectors <- function(scores, title_label) {
  env <- scores$env
  gen <- scores$gen
  ax <- axis_labels(scores$varexpl)
  lims <- get_limits(env, gen)
  
  ggplot() +
    geom_point(data = gen, aes(PC1, PC2),
               color = "grey60", size = 1.2, alpha = 0.5) +
    geom_segment(data = env, aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Status),
                 linewidth = 1.2) +
    geom_point(data = env, aes(PC1, PC2, color = Status),
               shape = 18, size = 6) +
    geom_text_repel(data = env, aes(PC1, PC2, label = Label, color = Status),
                    fontface = "bold", size = 5.5, show.legend = FALSE,
                    box.padding = 0.6, point.padding = 0.4,
                    max.overlaps = 20) +
    scale_color_manual(values = status_colors) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(title = paste0("Relationship Among ", title_label, "s"),
         subtitle = subtitle_text,
         x = ax$x, y = ax$y, color = "") +
    yan_theme
}

# =============================================================================
# BIPLOT 2: WHICH-WON-WHERE
# =============================================================================
biplot_which_won <- function(scores, title_label) {
  env <- scores$env
  gen <- scores$gen
  ax <- axis_labels(scores$varexpl)
  lims <- get_limits(env, gen)
  
  hull_idx <- chull(gen$PC1, gen$PC2)
  hull_closed <- c(hull_idx, hull_idx[1])
  hull_df <- gen[hull_closed, ]
  vertex_gen <- gen[hull_idx, ]
  
  ext <- max(abs(c(lims$xlim, lims$ylim))) * 2
  sectors <- data.frame()
  for (i in seq_along(hull_idx)) {
    j <- ifelse(i == length(hull_idx), 1, i + 1)
    mx <- (gen$PC1[hull_idx[i]] + gen$PC1[hull_idx[j]]) / 2
    my <- (gen$PC2[hull_idx[i]] + gen$PC2[hull_idx[j]]) / 2
    len <- sqrt(mx^2 + my^2)
    if (len > 0.01) {
      sectors <- rbind(sectors, data.frame(
        x = 0, y = 0,
        xend = mx / len * ext,
        yend = my / len * ext
      ))
    }
  }
  
  ggplot() +
    geom_segment(data = sectors, aes(x = x, y = y, xend = xend, yend = yend),
                 linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_polygon(data = hull_df, aes(PC1, PC2),
                 fill = NA, color = col_mega, linewidth = 1.2) +
    geom_point(data = gen, aes(PC1, PC2),
               color = "grey60", size = 1.2, alpha = 0.5) +
    geom_text_repel(data = vertex_gen, aes(PC1, PC2, label = Label),
                    color = col_mega, fontface = "bold", size = 4.5,
                    box.padding = 0.4, max.overlaps = 20) +
    geom_point(data = env, aes(PC1, PC2, color = Status),
               shape = 18, size = 6) +
    geom_text_repel(data = env, aes(PC1, PC2, label = Label, color = Status),
                    fontface = "bold", size = 5.5, show.legend = FALSE,
                    box.padding = 0.6, point.padding = 0.4,
                    max.overlaps = 20) +
    scale_color_manual(values = status_colors) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(title = paste0("Which-Won-Where (", title_label, ")"),
         subtitle = subtitle_text,
         x = ax$x, y = ax$y, color = "") +
    yan_theme
}

# =============================================================================
# BIPLOT 3: MEAN VS. STABILITY
# =============================================================================
biplot_mean_stability <- function(scores, title_label) {
  env <- scores$env
  gen <- scores$gen
  ax <- axis_labels(scores$varexpl)
  lims <- get_limits(env, gen)
  
  aec_x <- mean(env$PC1)
  aec_y <- mean(env$PC2)
  aec_len <- sqrt(aec_x^2 + aec_y^2)
  u_aec_x <- aec_x / aec_len
  u_aec_y <- aec_y / aec_len
  u_perp_x <- -u_aec_y
  u_perp_y <- u_aec_x
  
  gen <- gen %>%
    mutate(
      proj_mean = PC1 * u_aec_x + PC2 * u_aec_y,
      proj_stab = PC1 * u_perp_x + PC2 * u_perp_y
    )
  
  top_mean <- gen %>% slice_max(proj_mean, n = 5)
  bot_mean <- gen %>% slice_min(proj_mean, n = 5)
  top_unstable <- gen %>% slice_max(abs(proj_stab), n = 5)
  label_gen <- bind_rows(top_mean, bot_mean, top_unstable) %>%
    distinct(Label, .keep_all = TRUE)
  
  ext <- max(abs(c(lims$xlim, lims$ylim))) * 1.5
  
  ggplot() +
    geom_segment(aes(x = -ext * u_aec_x, y = -ext * u_aec_y,
                     xend = ext * u_aec_x, yend = ext * u_aec_y),
                 color = col_mega, linewidth = 1.0,
                 arrow = arrow(length = unit(0.2, "cm"), ends = "last")) +
    geom_segment(aes(x = aec_x - ext * u_perp_x * 0.4,
                     y = aec_y - ext * u_perp_y * 0.4,
                     xend = aec_x + ext * u_perp_x * 0.4,
                     yend = aec_y + ext * u_perp_y * 0.4),
                 color = col_mega, linewidth = 0.8, linetype = "dashed") +
    geom_point(aes(x = aec_x, y = aec_y),
               shape = 21, size = 4, fill = col_mega, color = "white") +
    geom_point(data = gen, aes(PC1, PC2),
               color = "grey60", size = 1.2, alpha = 0.5) +
    geom_text_repel(data = label_gen, aes(PC1, PC2, label = Label),
                    color = col_mega, fontface = "bold", size = 4,
                    box.padding = 0.3, max.overlaps = 20) +
    geom_point(data = env, aes(PC1, PC2, color = Status),
               shape = 18, size = 6) +
    geom_text_repel(data = env, aes(PC1, PC2, label = Label, color = Status),
                    fontface = "bold", size = 5.5, show.legend = FALSE,
                    box.padding = 0.6, point.padding = 0.4,
                    max.overlaps = 20) +
    scale_color_manual(values = status_colors) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(title = paste0("Mean vs. Stability (", title_label, ")"),
         subtitle = subtitle_text,
         x = ax$x, y = ax$y, color = "") +
    yan_theme
}

# =============================================================================
# BIPLOT 4: REPRESENTATIVENESS VS. DISCRIMINATING
# =============================================================================
biplot_repr_disc <- function(scores, title_label) {
  env <- scores$env
  gen <- scores$gen
  ax <- axis_labels(scores$varexpl)
  lims <- get_limits(env, gen)
  
  aec_x <- mean(env$PC1)
  aec_y <- mean(env$PC2)
  aec_len <- sqrt(aec_x^2 + aec_y^2)
  
  env <- env %>%
    mutate(vec_length = sqrt(PC1^2 + PC2^2))
  
  max_len <- max(env$vec_length) * 1.1
  circle_radii <- seq(max_len / 3, max_len, length.out = 3)
  
  circle_dfs <- map_dfr(circle_radii, function(r) {
    data.frame(theta = seq(0, 2 * pi, length.out = 100), r = r) %>%
      mutate(x = r * cos(theta), y = r * sin(theta), group = as.character(r))
  })
  
  ext <- max_len * 1.5
  
  ggplot() +
    geom_path(data = circle_dfs, aes(x, y, group = group),
              color = "grey80", linewidth = 0.3) +
    geom_segment(aes(x = 0, y = 0,
                     xend = aec_x / aec_len * ext,
                     yend = aec_y / aec_len * ext),
                 color = col_mega, linewidth = 1.0,
                 arrow = arrow(length = unit(0.2, "cm"))) +
    geom_segment(data = env, aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Status),
                 linewidth = 1.0) +
    geom_point(data = env, aes(PC1, PC2, color = Status),
               shape = 18, size = 6) +
    geom_text_repel(data = env, aes(PC1, PC2, label = Label, color = Status),
                    fontface = "bold", size = 5.5, show.legend = FALSE,
                    box.padding = 0.6, point.padding = 0.4,
                    max.overlaps = 20) +
    geom_point(data = gen, aes(PC1, PC2),
               color = "grey60", size = 1.2, alpha = 0.4) +
    scale_color_manual(values = status_colors) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey50", linewidth = 0.4) +
    coord_cartesian(xlim = lims$xlim, ylim = lims$ylim) +
    labs(title = paste0("Representativeness vs. Discriminating (", title_label, ")"),
         subtitle = subtitle_text,
         x = ax$x, y = ax$y, color = "") +
    yan_theme
}

# =============================================================================
# GENERATE AND SAVE ALL
# =============================================================================
generate_and_save <- function(scores, label) {
  plots <- list(
    env_vectors = biplot_env_vectors(scores, label),
    which_won   = biplot_which_won(scores, label),
    mean_stab   = biplot_mean_stability(scores, label),
    repr_disc   = biplot_repr_disc(scores, label)
  )
  
  for (nm in names(plots)) {
    ggsave(
      filename = paste0("GGE_", tolower(label), "_", nm, ".png"),
      plot = plots[[nm]], width = 10, height = 10, dpi = 600
    )
  }
  
  plots
}

plots_env  <- generate_and_save(scores_env,  "Environment")
plots_year <- generate_and_save(scores_year, "Year")
plots_loc  <- generate_and_save(scores_loc,  "Location")

cat("\n=== All plots saved ===\n")
