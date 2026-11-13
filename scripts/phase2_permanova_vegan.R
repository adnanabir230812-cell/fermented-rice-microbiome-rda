# ==============================================================================
# Phase 2: Two-Way PERMANOVA (adonis2), PERMDISP (betadisper), and ANOSIM in R
# Package Requirement: vegan (v2.6-4)
# ==============================================================================
if (!requireNamespace("vegan", quietly = TRUE)) {
    install.packages("vegan", repos = "https://cloud.r-project.org")
}
library(vegan)

dist_dir <- file.path("data", "distance_matrices")
if (!dir.exists(dist_dir)) {
    dist_dir <- file.path("..", "data", "distance_matrices")
}

read_dm <- function(tsv_path) {
    df <- read.table(tsv_path, sep="\t", header=TRUE, row.names=1, check.names=FALSE)
    as.dist(df)
}

metrics <- list(
    "Weighted_UniFrac" = "weighted_unifrac.tsv",
    "Bray_Curtis" = "bray_curtis.tsv",
    "Unweighted_UniFrac" = "unweighted_unifrac.tsv",
    "Jaccard" = "jaccard.tsv"
)

for (m_name in names(metrics)) {
    f_path <- file.path(dist_dir, metrics[[m_name]])
    dm <- read_dm(f_path)
    samples <- attr(dm, "Labels")
    meta <- data.frame(
        Sample = samples,
        Variety = sub("-(12|16)$", "", samples),
        Time = ifelse(grepl("-16$", samples), "16h", "12h"),
        row.names = samples
    )
    cat("====================================================================\n")
    cat("Metric:", m_name, "\n")
    perm_res <- adonis2(dm ~ Variety + Time, data = meta, permutations = 9999, by = "terms")
    print(perm_res)
    disp_var <- betadisper(dm, meta$Variety)
    cat("PERMDISP (Variety) F =", anova(disp_var)$`F value`[1], "P =", anova(disp_var)$`Pr(>F)`[1], "\n")
}
