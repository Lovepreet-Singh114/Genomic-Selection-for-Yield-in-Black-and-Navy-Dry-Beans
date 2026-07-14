## ============================================================================ ##
################# Calculate BLUEs (Best Linear Unbiased Estimates)  #############
## ============================================================================ ##

# Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# 1. Load Required Packages
# ------------------------------------------------------------------------------
install.packages('sommer')
library(tidyverse)   # Data manipulation
library(sommer)      # Mixed model fitting (mmer function)

# ------------------------------------------------------------------------------
# 2. Load Data
# ------------------------------------------------------------------------------

# Load cleaned phenotype data
pheno <- read_csv("pheno_clean.csv", show_col_types = FALSE)

# Convert to factors
pheno <- pheno %>%
  mutate(
    Name   = as.factor(Name),
    Year   = as.factor(Year),
    Location = as.factor(Location),
    Experiment_Name = as.factor(Experiment_Name),
    REP    = as.factor(REP),
    IBLK   = as.factor(IBLK)
  )

# Display data summary
cat("  Observations:", nrow(pheno), "\n")
cat("  Genotypes:", length(unique(pheno$Name)), "\n")
cat("  Years:", paste(levels(pheno$Year), collapse = ", "), "\n")
cat("  Locations:", paste(levels(pheno$Location), collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# 3. Calculate Overall BLUEs (Across All Years and Locations)
# ------------------------------------------------------------------------------

model_overall <- mmer(
  fixed  = yield ~ Name,
  random = ~ Experiment_Name + REP:Experiment_Name + IBLK:REP:Experiment_Name,
  rcov   = ~ units,
  data   = pheno,
  verbose = FALSE
)

# Display variance components
print(summary(model_overall)$varcomp)

# Extract BLUEs from fixed effects
beta_df <- model_overall$Beta

# Get intercept (reference genotype)
intercept <- beta_df$Estimate[beta_df$Effect == "(Intercept)"]
ref_geno  <- levels(pheno$Name)[1]

cat("\nReference genotype:", ref_geno, "\n")
cat("Intercept (reference BLUE):", round(intercept, 3), "\n")

# Extract genotype effects
geno_effects <- beta_df %>%
  filter(Effect != "(Intercept)") %>%
  mutate(
    Name = gsub("^Name", "", Effect),
    BLUE = Estimate + intercept
  ) %>%
  dplyr::select(Name, BLUE)

# Add reference genotype
ref_row <- data.frame(Name = ref_geno, BLUE = intercept)
blues_overall <- rbind(ref_row, geno_effects)

# Sort by genotype name
blues_overall <- blues_overall %>%
  arrange(Name)

# Summary statistics
cat("\nOverall BLUEs Summary:\n")
cat("  Number of genotypes:", nrow(blues_overall), "\n")
cat("  Mean yield:", round(mean(blues_overall$BLUE), 2), "\n")
cat("  SD yield:", round(sd(blues_overall$BLUE), 2), "\n")
cat("  Range:", round(min(blues_overall$BLUE), 2), "to", 
    round(max(blues_overall$BLUE), 2), "\n")

# ------------------------------------------------------------------------------
# 4. Calculate Year-Specific BLUEs 
# ------------------------------------------------------------------------------

# Get unique years
years <- sort(unique(as.character(pheno$Year)))
cat("Years:", paste(years, collapse = ", "), "\n\n")

# Initialize list to store results
blues_yearly_list <- list()

# Loop through each year
for (yr in years) {
  
  cat("Processing Year:", yr, "\n")
  
  # Subset data for this year
  pheno_yr <- pheno %>%
    filter(Year == yr) %>%
    droplevels()
  
  cat("  Observations:", nrow(pheno_yr), "\n")
  cat("  Genotypes:", length(unique(pheno_yr$Name)), "\n")
  cat("  Experiments:", length(unique(pheno_yr$Experiment_Name)), "\n")
  
  # Check experimental design structure
  n_expt <- length(unique(pheno_yr$Experiment_Name))
  n_rep  <- length(unique(pheno_yr$REP))
  n_iblk <- length(unique(pheno_yr$IBLK))
  
  # Fit model based on data structure
  tryCatch({
    
    if (n_expt > 1 && n_iblk > 1) {
      # Multiple experiments with incomplete blocks
      model_yr <- mmer(
        fixed  = yield ~ Name,
        random = ~ Experiment_Name + REP:Experiment_Name + IBLK:REP:Experiment_Name,
        rcov   = ~ units,
        data   = pheno_yr,
        verbose = FALSE
      )
    } else if (n_iblk > 1) {
      # Single experiment with incomplete blocks
      model_yr <- mmer(
        fixed  = yield ~ Name,
        random = ~ REP + IBLK:REP,
        rcov   = ~ units,
        data   = pheno_yr,
        verbose = FALSE
      )
    } else {
      # Simple RCBD
      model_yr <- mmer(
        fixed  = yield ~ Name,
        random = ~ REP,
        rcov   = ~ units,
        data   = pheno_yr,
        verbose = FALSE
      )
    }
    
    # Extract BLUEs
    beta_yr <- model_yr$Beta
    
    # Get intercept and reference genotype
    intercept_yr <- beta_yr$Estimate[beta_yr$Effect == "(Intercept)"]
    ref_geno_yr  <- levels(pheno_yr$Name)[1]
    
    # Extract genotype effects
    geno_effects_yr <- beta_yr %>%
      filter(Effect != "(Intercept)") %>%
      mutate(
        Name = gsub("^Name", "", Effect),
        BLUE = Estimate + intercept_yr
      ) %>%
      dplyr::select(Name, BLUE)
    
    # Add reference genotype
    ref_row_yr <- data.frame(Name = ref_geno_yr, BLUE = intercept_yr)
    blues_yr <- rbind(ref_row_yr, geno_effects_yr)
    
    # Add year column
    blues_yr$Year <- yr
    
    # Store in list
    blues_yearly_list[[yr]] <- blues_yr
    
    cat("  BLUEs calculated for", nrow(blues_yr), "genotypes\n")
    cat("  BLUE range:", round(min(blues_yr$BLUE), 2), "to", 
        round(max(blues_yr$BLUE), 2), "\n\n")
    
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n\n")
  })
}

# Combine all year-specific BLUEs
blues_yearly <- do.call(rbind, blues_yearly_list)
rownames(blues_yearly) <- NULL

# Summary
cat("--- Year-Specific BLUEs Summary ---\n")
print(
  blues_yearly %>%
    group_by(Year) %>%
    summarise(
      n_geno = n(),
      mean_BLUE = round(mean(BLUE), 2),
      sd_BLUE = round(sd(BLUE), 2),
      min_BLUE = round(min(BLUE), 2),
      max_BLUE = round(max(BLUE), 2),
      .groups = "drop"
    )
)

# ------------------------------------------------------------------------------
# 5. Create Wide Format for Year-Specific BLUEs
# ------------------------------------------------------------------------------

blues_yearly_wide <- blues_yearly %>%
  pivot_wider(
    names_from = Year,
    values_from = BLUE,
    names_prefix = "BLUE_"
  )

cat("Wide format dimensions:", nrow(blues_yearly_wide), "genotypes x", 
    ncol(blues_yearly_wide), "columns\n")

# ------------------------------------------------------------------------------
# 6. Check Genotype Overlap Across Years
# ------------------------------------------------------------------------------

# Create presence/absence matrix
geno_year_presence <- blues_yearly %>%
  dplyr::select(Name, Year) %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = Year,
    values_from = present,
    values_fill = 0
  )

# Count years per genotype
year_cols <- grep("^20", names(geno_year_presence), value = TRUE)
geno_year_presence$n_years <- rowSums(geno_year_presence[, year_cols])

cat("Genotypes by number of years present:\n")
print(table(geno_year_presence$n_years))

# Genotypes in all years
n_all_years <- sum(geno_year_presence$n_years == length(years))
cat("\nGenotypes present in ALL", length(years), "years:", n_all_years, "\n")

# Genotypes in at least 2 years
n_min_2 <- sum(geno_year_presence$n_years >= 2)
cat("Genotypes present in at least 2 years:", n_min_2, "\n")

# ------------------------------------------------------------------------------
# 7. Correlation Between Years
# ------------------------------------------------------------------------------

# Get genotypes present in all years for correlation
geno_all_years <- geno_year_presence$Name[geno_year_presence$n_years == length(years)]

if (length(geno_all_years) >= 10) {
  
  blues_corr <- blues_yearly_wide %>%
    filter(Name %in% geno_all_years) %>%
    dplyr::select(starts_with("BLUE_"))
  
  cor_matrix <- cor(blues_corr, use = "pairwise.complete.obs")
  
  cat("Correlation matrix (genotypes in all years, n =", length(geno_all_years), "):\n")
  print(round(cor_matrix, 3))
  
} else {
  cat("Insufficient genotypes in all years for correlation analysis\n")
}

# ------------------------------------------------------------------------------
# 8. Save Output Files
# ------------------------------------------------------------------------------

# Save overall BLUEs
write_csv(blues_overall, "blues_overall.csv")


# Save year-specific BLUEs (long format)
write_csv(blues_yearly, "blues_yearly.csv")

# Save year-specific BLUEs (wide format)
write_csv(blues_yearly_wide, "blues_yearly_wide.csv")

# Save genotype-year presence matrix
write_csv(geno_year_presence, "genotype_year_presence.csv")

