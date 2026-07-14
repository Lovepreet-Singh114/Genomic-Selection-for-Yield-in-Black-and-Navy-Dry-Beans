#!/usr/bin/env Rscript
# ==============================================================================
# Genotypic Data Quality Control, Structure, LD, GRM, and Haplotype Analysis
# Dry Bean Genomic Selection Pipeline
# Market classes: Black and Navy
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(rrBLUP)
  library(SNPRelate)
  library(ggplot2)
  library(reshape2)
})

set.seed(123)

# ------------------------------------------------------------------------------
# 1. Input files
# ------------------------------------------------------------------------------

GENO_FILE <- "geno_f.csv"        # Genotype ID + SNPs (0/1/2)
META_FILE <- "metadata.csv"      # Genotype, MarketClass
MAP_FILE  <- "snp_map.csv"                # SNP, chr, pos

OUT_DIR <- "Genotype_QC"
dir.create(OUT_DIR, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 2. Load genotype data
# ------------------------------------------------------------------------------

geno <- fread(GENO_FILE)
geno_id <- geno[[1]]
G <- as.matrix(geno[, -1])
rownames(G) <- geno_id

meta <- fread(META_FILE)
meta <- meta[MarketClass %in% c("Black", "Navy")]

keep_ids <- intersect(rownames(G), meta$Genotype)
G <- G[keep_ids, ]
meta <- meta[match(keep_ids, meta$Genotype), ]

# ------------------------------------------------------------------------------
# 3. Marker-level QC
# ------------------------------------------------------------------------------

snp_call_rate <- colMeans(!is.na(G))

calc_maf <- function(x) {
  p <- mean(x, na.rm = TRUE) / 2
  min(p, 1 - p)
}

maf <- apply(G, 2, calc_maf)

qc_snps <- which(
  snp_call_rate >= 0.90 &
  maf >= 0.05
)

G_filt <- G[, qc_snps]

# ------------------------------------------------------------------------------
# 4. Individual-level QC
# ------------------------------------------------------------------------------

ind_call_rate <- rowMeans(!is.na(G_filt))
heterozygosity <- rowMeans(G_filt == 1, na.rm = TRUE)

qc_ind <- which(ind_call_rate >= 0.90)

G_filt <- G_filt[qc_ind, ]
meta <- meta[qc_ind, ]

# ------------------------------------------------------------------------------
# 5. Mean imputation
# ------------------------------------------------------------------------------

G_imp <- G_filt
for (j in seq_len(ncol(G_imp))) {
  mu <- mean(G_imp[, j], na.rm = TRUE)
  G_imp[is.na(G_imp[, j]), j] <- mu
}

# ------------------------------------------------------------------------------
# 6. Global PCA (all market classes)
# ------------------------------------------------------------------------------

K_all <- A.mat(G_imp)
eig_all <- eigen(K_all)

PC_all <- data.frame(
  Genotype = rownames(G_imp),
  PC1 = eig_all$vectors[, 1],
  PC2 = eig_all$vectors[, 2]
)

PC_all <- merge(PC_all, meta, by = "Genotype")

ggsave(
  file.path(OUT_DIR, "PCA_All_MarketClasses.png"),
  ggplot(PC_all, aes(PC1, PC2, color = MarketClass)) +
    geom_point(size = 2, alpha = 0.8) +
    theme_classic(),
  width = 6, height = 5, dpi = 300
)

# ------------------------------------------------------------------------------
# 7. Market-class–specific PCA
# ------------------------------------------------------------------------------

for (cls in unique(meta$MarketClass)) {

  idx <- meta$MarketClass == cls
  G_cls <- G_imp[idx, ]

  if (nrow(G_cls) < 10) next

  K_cls <- A.mat(G_cls)
  eig_cls <- eigen(K_cls)

  PC_cls <- data.frame(
    Genotype = rownames(G_cls),
    PC1 = eig_cls$vectors[, 1],
    PC2 = eig_cls$vectors[, 2]
  )

  p <- ggplot(PC_cls, aes(PC1, PC2)) +
    geom_point(size = 2) +
    theme_classic() +
    labs(title = paste("PCA:", cls))

  ggsave(
    file.path(OUT_DIR, paste0("PCA_", cls, ".png")),
    p, width = 5, height = 4, dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 8. GRM diagnostics (within and between market classes)
# ------------------------------------------------------------------------------

K_long <- melt(K_all)
colnames(K_long) <- c("G1", "G2", "Rel")

K_long <- merge(K_long, meta, by.x = "G1", by.y = "Genotype")
K_long <- merge(K_long, meta, by.x = "G2", by.y = "Genotype",
                suffixes = c("_1", "_2"))

K_long$ClassPair <- ifelse(
  K_long$MarketClass_1 == K_long$MarketClass_2,
  K_long$MarketClass_1,
  "Between"
)

ggsave(
  file.path(OUT_DIR, "GRM_By_ClassPair.png"),
  ggplot(K_long, aes(Rel, fill = ClassPair)) +
    geom_density(alpha = 0.5) +
    theme_classic(),
  width = 6, height = 5, dpi = 300
)

# ------------------------------------------------------------------------------
# 9. Market-class–specific LD decay
# ------------------------------------------------------------------------------

map <- fread(MAP_FILE)
map <- map[match(colnames(G_filt), map$SNP)]

for (cls in unique(meta$MarketClass)) {

  idx <- meta$MarketClass == cls
  G_cls <- G_imp[idx, ]

  if (nrow(G_cls) < 10) next

  gds_file <- file.path(OUT_DIR, paste0("geno_", cls, ".gds"))

  snpgdsCreateGeno(
    gds_file,
    genmat = G_cls,
    sample.id = rownames(G_cls),
    snp.id = map$SNP,
    snp.chromosome = map$chr,
    snp.position = map$pos,
    snpfirstdim = FALSE
  )

  genofile <- snpgdsOpen(gds_file)

  ld <- snpgdsLDMat(
    genofile,
    method = "r",
    slide = 500,
    num.thread = 4
  )

  snpgdsClose(genofile)

  ld_df <- data.frame(
    Distance = abs(ld$snp.position[,1] - ld$snp.position[,2]) / 1e6,
    r2 = ld$LD^2
  )

  p <- ggplot(ld_df, aes(Distance, r2)) +
    geom_point(alpha = 0.1) +
    geom_smooth(method = "loess") +
    theme_classic() +
    labs(title = paste("LD decay:", cls),
         x = "Physical distance (Mb)",
         y = expression(r^2))

  ggsave(
    file.path(OUT_DIR, paste0("LD_Decay_", cls, ".png")),
    p, width = 6, height = 5, dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 10. Haplotype-based representation (LD blocks)
# ------------------------------------------------------------------------------

# Simple LD-based haplotype blocks using SNPRelate
geno_gds <- file.path(OUT_DIR, "geno_haplo.gds")

snpgdsCreateGeno(
  geno_gds,
  genmat = G_imp,
  sample.id = rownames(G_imp),
  snp.id = map$SNP,
  snp.chromosome = map$chr,
  snp.position = map$pos,
  snpfirstdim = FALSE
)

genofile <- snpgdsOpen(geno_gds)

hap_blocks <- snpgdsLDpruning(
  genofile,
  ld.threshold = 0.8,
  autosome.only = FALSE
)

snpgdsClose(genofile)

hap_snp_ids <- unlist(hap_blocks, use.names = FALSE)
G_hap <- G_imp[, hap_snp_ids]

# Haplotype-based GRM
K_hap <- A.mat(G_hap)

# ------------------------------------------------------------------------------
# 11. Save outputs for downstream GS comparison
# ------------------------------------------------------------------------------

saveRDS(
  list(
    SNP_Genotypes = G_filt,
    SNP_Imputed = G_imp,
    SNP_GRM = K_all,
    Haplotype_SNPs = hap_snp_ids,
    Haplotype_GRM = K_hap,
    Metadata = meta
  ),
  file = file.path(OUT_DIR, "Genotype_QC_Extended_Objects.rds")
)

cat("Extended genotypic characterization completed successfully.\n")
