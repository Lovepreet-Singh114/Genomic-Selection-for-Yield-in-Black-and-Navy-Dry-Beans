################################################################################
# Cross-Market-Class Decomposition Figure (Hardcoded Values)
################################################################################

library(tidyverse)
library(patchwork)

# =============================================================================
# PANEL A: ACCURACY DOT PLOT
# =============================================================================
dotplot_data <- data.frame(
  Scenario = c("Navy\n(n=184)", "Black\n(n=299)", 
               "Whole\n(Navy N)", "Whole\n(Black N)", "Whole\n(n=483)"),
  Accuracy = c(0.323, 0.603, 0.601, 0.656, 0.668),
  Color = c("baseline", "baseline", "mixed", "mixed", "mixed")
)

dotplot_data$Scenario <- factor(dotplot_data$Scenario,
                                levels = rev(c("Navy\n(n=184)", "Black\n(n=299)",
                                               "Whole\n(Navy N)", "Whole\n(Black N)", 
                                               "Whole\n(n=483)")))

pA <- ggplot(dotplot_data, aes(x = Accuracy, y = Scenario, color = Color)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey60") +
  geom_point(size = 5) +
  geom_text(aes(label = round(Accuracy, 3)), hjust = -0.3, size = 4, fontface = "bold") +
  scale_color_manual(values = c("baseline" = "#1a3a4a", "mixed" = "#7b2d8e")) +
  labs(title = "(A) rrBLUP Prediction Accuracy",
       x = "Prediction Accuracy (r)", y = "") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  xlim(0.25, 0.8)

# =============================================================================
# PANEL B: BLACK DECOMPOSITION
# =============================================================================
bar_black <- data.frame(
  Step = c("Black\nBaseline", "Shared Genetic\nArchitecture", "Sample Size\nIncrease", "Whole\nPopulation"),
  Value = c(0.603, 0.053, 0.012, 0.668),
  Type = c("cumulative", "gain", "gain", "cumulative"),
  ymin = c(0, 0.603, 0.656, 0),
  ymax = c(0.603, 0.656, 0.668, 0.668)
)

bar_black$Step <- factor(bar_black$Step, 
                         levels = c("Black\nBaseline", "Shared Genetic\nArchitecture", 
                                    "Sample Size\nIncrease", "Whole\nPopulation"))

pB <- ggplot(bar_black, aes(x = Step)) +
  geom_rect(aes(xmin = as.numeric(Step) - 0.35, xmax = as.numeric(Step) + 0.35,
                ymin = ymin, ymax = ymax,
                fill = interaction(Step, Type))) +
  geom_text(data = bar_black %>% filter(Type == "gain"),
            aes(y = (ymin + ymax) / 2, 
                label = paste0("+", round(Value, 3))),
            size = 4, fontface = "bold", color = "white") +
  geom_text(data = bar_black %>% filter(Type == "cumulative"),
            aes(y = ymax + 0.01, label = round(Value, 3)),
            size = 4) +
  scale_fill_manual(values = c(
    "Black\nBaseline.cumulative" = "#1a3a4a",
    "Shared Genetic\nArchitecture.gain" = "#7b2d8e",
    "Sample Size\nIncrease.gain" = "#b8a042",
    "Whole\nPopulation.cumulative" = "#7b2d8e"
  )) +
  labs(title = expression("(B) Black: r +0.065 total gain"),
       x = "", y = "Prediction Accuracy (r)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  ylim(0, 0.72)

# =============================================================================
# PANEL C: NAVY DECOMPOSITION
# =============================================================================
bar_navy <- data.frame(
  Step = c("Navy\nBaseline", "Shared Genetic\nArchitecture", "Sample Size\nIncrease", "Whole\nPopulation"),
  Value = c(0.323, 0.279, 0.067, 0.668),
  Type = c("cumulative", "gain", "gain", "cumulative"),
  ymin = c(0, 0.323, 0.601, 0),
  ymax = c(0.323, 0.601, 0.668, 0.668)
)

bar_navy$Step <- factor(bar_navy$Step, 
                        levels = c("Navy\nBaseline", "Shared Genetic\nArchitecture", 
                                   "Sample Size\nIncrease", "Whole\nPopulation"))

pC <- ggplot(bar_navy, aes(x = Step)) +
  geom_rect(aes(xmin = as.numeric(Step) - 0.35, xmax = as.numeric(Step) + 0.35,
                ymin = ymin, ymax = ymax,
                fill = interaction(Step, Type))) +
  geom_text(data = bar_navy %>% filter(Type == "gain"),
            aes(y = (ymin + ymax) / 2, 
                label = paste0("+", round(Value, 3))),
            size = 4, fontface = "bold", color = "white") +
  geom_text(data = bar_navy %>% filter(Type == "cumulative"),
            aes(y = ymax + 0.01, label = round(Value, 3)),
            size = 4) +
  scale_fill_manual(values = c(
    "Navy\nBaseline.cumulative" = "#1a3a4a",
    "Shared Genetic\nArchitecture.gain" = "#7b2d8e",
    "Sample Size\nIncrease.gain" = "#b8a042",
    "Whole\nPopulation.cumulative" = "#7b2d8e"
  )) +
  labs(title = expression("(C) Navy: r +0.345 total gain"),
       x = "", y = "Prediction Accuracy (r)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none") +
  ylim(0, 0.72)

# =============================================================================
# PANEL D: SHARED GENETIC ARCHITECTURE
# =============================================================================
cor_data <- data.frame(
  Tier = rep(c("All\nMarkers", "Top 10%", "Top 5%", "Top 1%"), each = 2),
  Metric = rep(c("Effect Correlation", "Same-Sign Proportion"), 4),
  Value = c(0.14, 0.53, 0.41, 0.75, 0.54, 0.83, 0.65, 0.90)
)

cor_data$Tier <- factor(cor_data$Tier, 
                        levels = c("All\nMarkers", "Top 10%", "Top 5%", "Top 1%"))

pD <- ggplot(cor_data, aes(x = Tier, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = round(Value, 2)),
            position = position_dodge(width = 0.7), vjust = -0.3, size = 4, fontface = "bold") +
  scale_fill_manual(values = c("Effect Correlation" = "#1a3a4a", 
                               "Same-Sign Proportion" = "#7b2d8e")) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50") +
  labs(title = "(D) Shared Genetic Architecture: Black vs Navy",
       x = "Marker Effect Tier",
       y = "Correlation / Proportion",
       fill = "") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top") +
  ylim(0, 1.05)

# =============================================================================
# COMBINE AND SAVE 
# =============================================================================
p_all <- (pA | pB) / (pC | pD)

ggsave("decomposition_full_figure.jpeg", p_all, width = 16, height = 10, dpi = 300)

cat("\n=== Figure Complete ===\n")