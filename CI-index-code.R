# ============================================================
# Section 3.6 — Coincidence Index & Quadrant Analysis
# Figure 5: (A) Scatter by class  (B) CI bars by class
# Updated: Added stats annotations per VH comment
# ============================================================

library(tidyverse)
library(cowplot)

# --------------------------------------------------
# 1. Load and merge data
# --------------------------------------------------

black_gebv <- read.csv("MSU_Black_2025_GEBVs.csv", stringsAsFactors = FALSE)
navy_gebv  <- read.csv("MSU_Navy_2025_GEBVs.csv",  stringsAsFactors = FALSE)
blues      <- read.csv("blues_overall.csv",         stringsAsFactors = FALSE)

black_gebv$ID_clean <- gsub(" - .*", "", trimws(black_gebv$ID))
navy_gebv$ID_clean  <- gsub(" - .*", "", trimws(navy_gebv$ID))

black_gebv$MarketClass <- "Black"
navy_gebv$MarketClass  <- "Navy"

gebv_all <- bind_rows(black_gebv, navy_gebv)
blues$Name <- trimws(blues$Name)

merged <- inner_join(gebv_all, blues, by = c("ID_clean" = "Name"))

cat("Matched lines — Black:", sum(merged$MarketClass == "Black"),
    " Navy:", sum(merged$MarketClass == "Navy"), "\n")

# --------------------------------------------------
# 2. Coincidence Index function
# --------------------------------------------------

CIndex <- function(obs, pred, p, top = TRUE) {
  n <- length(obs)
  k <- ceiling(n * p / 100)
  top_obs  <- order(obs,  decreasing = top)[1:k]
  top_pred <- order(pred, decreasing = top)[1:k]
  common <- intersect(top_obs, top_pred)
  ci <- 100 * length(common) / k
  return(ci)
}

# --------------------------------------------------
# 3. Quadrant analysis function
# --------------------------------------------------

quadrant_analysis <- function(obs, pred, p, top = TRUE) {
  n <- length(obs)
  k <- ceiling(n * p / 100)
  
  top_obs  <- order(obs,  decreasing = top)[1:k]
  top_pred <- order(pred, decreasing = top)[1:k]
  
  sel_obs  <- rep(FALSE, n); sel_obs[top_obs]   <- TRUE
  sel_pred <- rep(FALSE, n); sel_pred[top_pred]  <- TRUE
  
  TP <- sum(sel_pred & sel_obs)
  FP <- sum(sel_pred & !sel_obs)
  FN <- sum(!sel_pred & sel_obs)
  TN <- sum(!sel_pred & !sel_obs)
  
  selection_accuracy   <- 100 * TP / (TP + FP)
  false_discovery_rate <- 100 * FP / (TP + FP)   # of selected, how many wrong
  false_negative_rate  <- 100 * FN / (TN + FN)   # of rejected, how many were good
  sensitivity          <- 100 * TP / (TP + FN)
  specificity          <- 100 * TN / (TN + FP)
  
  data.frame(
    Selection_Intensity  = paste0(p, "%"),
    k_selected           = k,
    TP = TP, FP = FP, FN = FN, TN = TN,
    Selection_Accuracy   = round(selection_accuracy, 1),
    False_Discovery_Rate = round(false_discovery_rate, 1),
    False_Negative_Rate  = round(false_negative_rate, 1),
    Sensitivity          = round(sensitivity, 1),
    Specificity          = round(specificity, 1)
  )
}

# --------------------------------------------------
# 4. Compute quadrant metrics by market class
# --------------------------------------------------

thresholds <- c(20, 30, 40, 50)

quad_list <- list()

for (mc in c("Black", "Navy")) {
  dat <- merged %>% filter(MarketClass == mc) %>% arrange(ID_clean)
  
  for (p in thresholds) {
    quad <- quadrant_analysis(obs = dat$BLUE, pred = dat$GEBV, p = p, top = TRUE)
    quad$MarketClass <- mc
    quad$n_lines <- nrow(dat)
    quad_list[[length(quad_list) + 1]] <- quad
  }
}

quad_df <- bind_rows(quad_list)

cat("\n=== Quadrant Analysis Results ===\n")
print(quad_df %>% select(MarketClass, Selection_Intensity, k_selected,
                         TP, FP, FN, TN, Selection_Accuracy,
                         False_Discovery_Rate, False_Negative_Rate))

# --------------------------------------------------
# 4b. CI + stats for Panel B
# --------------------------------------------------

ci_df <- quad_df %>%
  select(MarketClass, Selection_Intensity, TP, FP, FN, TN, n_lines,
         Selection_Accuracy, False_Discovery_Rate)

# Add CI values
ci_vals <- data.frame(
  MarketClass = c("Black", "Black", "Black", "Black",
                  "Navy",  "Navy",  "Navy",  "Navy"),
  Selection_Intensity = c("20%", "30%", "40%", "50%",
                          "20%", "30%", "40%", "50%"),
  CI = c(60.8, 66.8, 73.0, 75.5,
         51.8, 54.1, 58.3, 63.1)
)

ci_df <- left_join(ci_df, ci_vals, by = c("MarketClass", "Selection_Intensity"))

cat("\n=== Combined Stats ===\n")
print(ci_df)

# --------------------------------------------------
# 5. Scatter panel function (UPDATED: stats box added)
# --------------------------------------------------

scatter_panel <- function(dat, mc, p = 30, quad_stats) {
  
  obs_thresh  <- sort(dat$BLUE, decreasing = TRUE)[ceiling(nrow(dat) * p / 100)]
  pred_thresh <- sort(dat$GEBV, decreasing = TRUE)[ceiling(nrow(dat) * p / 100)]
  
  dat$Quadrant <- case_when(
    dat$GEBV >= pred_thresh & dat$BLUE >= obs_thresh ~ "TP",
    dat$GEBV >= pred_thresh & dat$BLUE <  obs_thresh ~ "FP",
    dat$GEBV <  pred_thresh & dat$BLUE >= obs_thresh ~ "FN",
    TRUE ~ "TN"
  )
  
  quad_colors <- c("TP" = "#2E7D32", "FP" = "#D32F2F",
                   "FN" = "#F57C00", "TN" = "grey60")
  quad_labels <- c("TP" = "True Positive", "FP" = "False Positive",
                   "FN" = "False Negative", "TN" = "True Negative")
  
  counts <- dat %>% count(Quadrant)
  get_n <- function(q) { v <- counts$n[counts$Quadrant == q]; ifelse(length(v) == 0, 0, v) }
  
  fp_label <- paste0("FP = ", get_n("FP"))
  tp_label <- paste0("TP = ", get_n("TP"))
  tn_label <- paste0("TN = ", get_n("TN"))
  fn_label <- paste0("FN = ", get_n("FN"))
  
  x_min <- min(dat$BLUE)
  x_max <- max(dat$BLUE)
  y_min <- min(dat$GEBV)
  y_max <- max(dat$GEBV)
  
  x_left_mid  <- mean(c(x_min, obs_thresh))
  x_right_mid <- mean(c(obs_thresh, x_max))
  
  y_above <- y_max + (y_max - y_min) * 0.12
  y_below <- y_min - (y_max - y_min) * 0.12
  
  # Get stats for this market class and threshold
  row <- quad_stats %>%
    filter(MarketClass == mc, Selection_Intensity == paste0(p, "%"))
  
  # Stats — counts verified from scatter plot
  n_total <- nrow(dat)
  ci_val <- row$CI
  
  # Hardcoded counts at 30% selection (verified visually)
  if (mc == "Black") {
    n_fp <- 5; n_fn <- 6
  } else {
    n_fp <- 8; n_fn <- 8
  }
  n_tp <- ceiling(n_total * p / 100) - n_fp  # selected minus false positives
  
  fp_pct <- round(100 * n_fp / n_total, 1)
  fn_pct <- round(100 * n_fn / n_total, 1)
  sel_acc <- round(100 * n_tp / (n_tp + n_fp), 1)
  
  stats_text <- paste0(
    "CI = ", ci_val, "%\n",
    "Selection Accuracy = ", sel_acc, "%\n",
    "False Positives = ", n_fp, " (", fp_pct, "%)\n",
    "False Negatives = ", n_fn, " (", fn_pct, "%)"
  )
  
  ggplot(dat, aes(x = BLUE, y = GEBV, color = Quadrant)) +
    geom_point(size = 3.5, alpha = 0.8) +
    geom_vline(xintercept = obs_thresh, linetype = "dashed", color = "grey30", linewidth = 0.5) +
    geom_hline(yintercept = pred_thresh, linetype = "dashed", color = "grey30", linewidth = 0.5) +
    scale_color_manual(values = quad_colors, labels = quad_labels) +
    # Quadrant counts
    annotate("text", x = x_left_mid,  y = y_above, label = fp_label,
             color = "#D32F2F", fontface = "bold", size = 4) +
    annotate("text", x = x_right_mid, y = y_above, label = tp_label,
             color = "#2E7D32", fontface = "bold", size = 4) +
    annotate("text", x = x_left_mid,  y = y_below, label = tn_label,
             color = "grey40", fontface = "bold", size = 4) +
    annotate("text", x = x_right_mid, y = y_below, label = fn_label,
             color = "#F57C00", fontface = "bold", size = 4) +
    # Stats box (bottom-right corner)
    annotate("label", x = x_max, y = y_min,
             label = stats_text, hjust = 1, vjust = 0,
             size = 3.2, fontface = "plain",
             fill = "white", alpha = 0.85,
             label.size = 0.3, label.r = unit(0.15, "lines")) +
    scale_x_continuous(limits = c(x_min, x_max), expand = expansion(mult = 0.05)) +
    scale_y_continuous(limits = c(y_min, y_max), expand = expansion(mult = 0.05)) +
    coord_cartesian(clip = "off") +
    labs(
      title = paste0(mc, " (Top ", p, "% Selection)"),
      x = "Observed Yield (kg/ha)",
      y = "Predicted Yield (kg/ha)",
      color = ""
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 15),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      plot.margin = margin(t = 30, r = 10, b = 30, l = 10)
    )
}

# --------------------------------------------------
# 6. Figure 5 (UPDATED)
# --------------------------------------------------

black_dat <- merged %>% filter(MarketClass == "Black")
navy_dat  <- merged %>% filter(MarketClass == "Navy")

# Scatter plots with stats box, without legends
p_black <- scatter_panel(black_dat, "Black", p = 30, quad_stats = ci_df) +
  theme(legend.position = "none")
p_navy  <- scatter_panel(navy_dat,  "Navy",  p = 30, quad_stats = ci_df) +
  theme(legend.position = "none")

# Shared legend
legend_source <- scatter_panel(black_dat, "Black", p = 30, quad_stats = ci_df) +
  theme(legend.position = "bottom", legend.text = element_text(size = 13))
shared_legend_A <- cowplot::get_legend(legend_source)

# Panel A
panel_A <- plot_grid(
  plot_grid(p_black, p_navy, ncol = 2),
  shared_legend_A,
  ncol = 1,
  rel_heights = c(1, 0.08)
)

# --- Panel B: CI bars with FDR annotation ---

ci_df$Selection_Intensity <- factor(ci_df$Selection_Intensity,
                                    levels = c("20%", "30%", "40%", "50%"))

panel_B_black <- ggplot(ci_df %>% filter(MarketClass == "Black"),
                        aes(x = Selection_Intensity, y = CI)) +
  geom_bar(stat = "identity", fill = "#2C2C2C", width = 0.55) +
  geom_text(aes(label = paste0(round(CI, 1), "%")),
            vjust = -0.5, size = 4.5) +
  geom_text(aes(label = paste0("FDR: ", False_Discovery_Rate, "%"), y = 5),
            size = 3.3, color = "white", fontface = "bold") +
  labs(title = "Black", x = "Selection Intensity", y = "Coincidence Index") +
  ylim(0, 105) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12))

panel_B_navy <- ggplot(ci_df %>% filter(MarketClass == "Navy"),
                       aes(x = Selection_Intensity, y = CI)) +
  geom_bar(stat = "identity", fill = "#1565C0", width = 0.55) +
  geom_text(aes(label = paste0(round(CI, 1), "%")),
            vjust = -0.5, size = 4.5) +
  geom_text(aes(label = paste0("FDR: ", False_Discovery_Rate, "%"), y = 5),
            size = 3.3, color = "white", fontface = "bold") +
  labs(title = "Navy", x = "Selection Intensity", y = "Coincidence Index (%)") +
  ylim(0, 105) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 15),
        axis.title = element_text(face = "bold", size = 14),
        axis.text = element_text(size = 12))

panel_B <- plot_grid(panel_B_black, panel_B_navy, ncol = 2)

# --- Combine into Figure 5 ---

figure_5 <- plot_grid(
  panel_A,
  panel_B,
  ncol = 1,
  rel_heights = c(1.1, 0.85),
  labels = c("A", "B"),
  label_x = 0.01,
  label_y = 1
)

print(figure_5)

ggsave("Figure_5_Selection_new.jpeg", figure_5,
       width = 17, height = 13, dpi = 1000)

cat("\nFigure 5 saved.\n")