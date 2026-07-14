#######################################################################################
###==============         Data Preparation     ===================================== ##
########################################################################################
                                             
### Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# 1. Load Required Packages
# ------------------------------------------------------------------------------

library(tidyverse)   # Data manipulation and visualization
library(data.table)  # Efficient data handling

# ------------------------------------------------------------------------------
# 2. Set Working Directory and File Paths
# ------------------------------------------------------------------------------

# Set working directory (modify as needed)
setwd("~/Parametric_analysis/1. Data Preparation")

# ------------------------------------------------------------------------------
# 3. Load Phenotypic Data
# ------------------------------------------------------------------------------

pheno_raw <- read_csv("raw_phenotypic_data.csv", show_col_types = FALSE)

# checking dimensions and basic info about data
cat("\nPhenotypic data dimensions:", nrow(pheno_raw), "rows x", ncol(pheno_raw), "columns\n")
cat("Column names:", paste(names(pheno_raw), collapse = ", "), "\n")

# ------------------------------------------------------------------------------
# 4. Explore Phenotypic Data Structure
# ------------------------------------------------------------------------------

# Check years
print(table(pheno_raw$Year))

# Check locations
print(table(pheno_raw$Location))

# Check experiments
print(table(pheno_raw$Experiment_Name))

# Check year x location combinations
print(table(pheno_raw$Year, pheno_raw$Location))

# Number of unique genotypes
n_geno_pheno <- length(unique(pheno_raw$Name))
cat("\nNumber of unique genotypes in phenotypic data:", n_geno_pheno, "\n")

# Yield summary
cat("\nYield summary statistics:\n")
print(summary(pheno_raw$yield))

# ------------------------------------------------------------------------------
# 5. Clean Phenotypic Data
# ------------------------------------------------------------------------------

# Convert columns to appropriate types
pheno_clean <- pheno_raw %>%
  mutate(
    Name   = as.factor(Name),
    Year   = as.factor(Year),
    Location = as.factor(Location),
    Experiment_Name = as.factor(Experiment_Name),
    REP    = as.factor(REP),
    IBLK   = as.factor(IBLK),
    yield  = as.numeric(yield)
  )

# Check for missing yield values
n_missing_yield <- sum(is.na(pheno_clean$yield))
cat("Missing yield values:", n_missing_yield, "\n")

# Remove rows with missing yield (if any)
if (n_missing_yield > 0) {
  pheno_clean <- pheno_clean %>% filter(!is.na(yield))
  cat("Removed", n_missing_yield, "rows with missing yield\n")
}

# ------------------------------------------------------------------------------
# 6. Load Genotypic Data
# ------------------------------------------------------------------------------

cat("\n--- Loading Genotypic Data ---\n")

geno_raw <- read_csv("geno_f.csv", show_col_types = FALSE)

# Extract genotype names (first column)
geno_names <- geno_raw$ID

# Convert to matrix (excluding ID column)
geno_matrix <- as.matrix(geno_raw[, -1])
rownames(geno_matrix) <- geno_names

# Display dimensions
n_lines   <- nrow(geno_matrix)
n_markers <- ncol(geno_matrix)

cat("Genotypic data dimensions:", n_lines, "lines x", n_markers, "markers\n")

# ------------------------------------------------------------------------------
# 7. Genotype Quality Control
# ------------------------------------------------------------------------------

# Check marker coding
cat("Marker value distribution:\n")
print(table(as.vector(geno_matrix)))

# 7.1 Check for missing values
missing_per_marker <- colMeans(is.na(geno_matrix))
missing_per_line   <- rowMeans(is.na(geno_matrix))

cat("\nMissing data summary:\n")
cat("  Markers with >10% missing:", sum(missing_per_marker > 0.10), "\n")
cat("  Lines with >10% missing:", sum(missing_per_line > 0.10), "\n")

# 7.2 Calculate Minor Allele Frequency (MAF)
geno_012 <- geno_matrix * 2  # Convert to 0, 1, 2 coding

calc_maf <- function(x) {
  p <- mean(x, na.rm = TRUE) / 2
  return(min(p, 1 - p))
}

maf <- apply(geno_012, 2, calc_maf)

cat("\nMAF distribution:\n")
cat("  Min:", round(min(maf), 4), "\n")
cat("  Max:", round(max(maf), 4), "\n")
cat("  Markers with MAF < 0.05:", sum(maf < 0.05), "\n")

# 7.3 Apply QC filters
maf_threshold <- 0.05
missing_threshold <- 0.10

markers_pass_maf     <- maf >= maf_threshold
markers_pass_missing <- missing_per_marker <= missing_threshold
markers_keep         <- markers_pass_maf & markers_pass_missing

cat("\n--- Marker Filtering Summary ---\n")
cat("  Total markers:", n_markers, "\n")
cat("  Failed MAF filter:", sum(!markers_pass_maf), "\n")
cat("  Failed missing filter:", sum(!markers_pass_missing), "\n")
cat("  Markers retained:", sum(markers_keep), "\n")

# Apply filter
geno_filtered <- geno_matrix[, markers_keep]

cat("\nFiltered genotype matrix:", nrow(geno_filtered), "lines x", ncol(geno_filtered), "markers\n")


# ------------------------------------------------------------------------------
# 8. Save Processed Data
# ------------------------------------------------------------------------------

cat("\n--- Saving Processed Data ---\n")

# Save cleaned phenotype data
write_csv(pheno_clean, "pheno_clean.csv")


# Save filtered genotype matrix
geno_to_save <- as.data.frame(geno_filtered)
geno_to_save <- cbind(Name = rownames(geno_filtered), geno_to_save)
write_csv(geno_to_save, "geno_clean.csv")

