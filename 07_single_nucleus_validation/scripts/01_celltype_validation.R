suppressPackageStartupMessages({
    library(limma)
})

# ============================================================
# Files
# ============================================================

count_file <- "raw/GSE138852_counts.csv.gz"
cov_file   <- "metadata/GSE138852_covariates.csv.gz"

bulk_file <- "../04_cohort_specific_dea/results/EC_AD_vs_Control_all_probes.tsv"

dir.create("results", showWarnings=FALSE)
dir.create("tables", showWarnings=FALSE)
dir.create("figures", showWarnings=FALSE)

# ============================================================
# Metadata
# ============================================================

cov <- read.csv(
    gzfile(cov_file),
    check.names=FALSE,
    stringsAsFactors=FALSE
)

# First unnamed column contains nucleus IDs
names(cov)[1] <- "cell_id"

names(cov)[
    names(cov)=="oupSample.batchCond"
] <- "condition"

names(cov)[
    names(cov)=="oupSample.cellType"
] <- "cell_type"

cov$condition[
    cov$condition=="ct"
] <- "Control"

cov$condition[
    cov$condition=="AD"
] <- "AD"

# Standardize cell-type names
cellmap <- c(
    "neuron"="Neuron",
    "astro"="Astrocyte",
    "mg"="Microglia",
    "oligo"="Oligodendrocyte",
    "OPC"="OPC",
    "endo"="Endothelial"
)

cov$cell_type_std <- unname(
    cellmap[cov$cell_type]
)

keep <- !is.na(cov$cell_type_std) &
        cov$condition %in% c("AD","Control")

cov <- cov[keep,,drop=FALSE]

cat("Nuclei retained:",nrow(cov),"\n\n")

cat("Cell-type counts:\n")
print(
    addmargins(
        table(cov$cell_type_std,cov$condition)
    )
)

# ============================================================
# Read counts
# ============================================================

cat("\nReading count matrix...\n")

dat <- read.csv(
    gzfile(count_file),
    check.names=FALSE,
    stringsAsFactors=FALSE
)

gene_id <- dat[[1]]

mat <- as.matrix(dat[,-1,drop=FALSE])
storage.mode(mat) <- "numeric"

rownames(mat) <- gene_id

cat(
    "Count matrix:",
    nrow(mat),"genes x",
    ncol(mat),"nuclei\n"
)

# ============================================================
# Match metadata and matrix
# ============================================================

common <- intersect(
    colnames(mat),
    cov$cell_id
)

cat("Matched nuclei:",length(common),"\n")

mat <- mat[,common,drop=FALSE]

cov <- cov[
    match(common,cov$cell_id),
    ,
    drop=FALSE
]

stopifnot(
    all(colnames(mat)==cov$cell_id)
)

# ============================================================
# Basic filtering
# ============================================================

detected <- rowSums(mat > 0)

mat <- mat[
    detected >= 10,
    ,
    drop=FALSE
]

cat(
    "Genes after detection filtering:",
    nrow(mat),"\n"
)

# ============================================================
# Library-size normalization
#
# logCPM used for validation-level expression comparison.
# ============================================================

lib <- colSums(mat)

lib[lib==0] <- 1

cpm <- t(
    t(mat) / lib * 1e6
)

logcpm <- log2(cpm + 1)

rm(dat,mat,cpm)
gc()

# ============================================================
# Cell-type-specific AD vs Control
# ============================================================

celltypes <- c(
    "Neuron",
    "Astrocyte",
    "Microglia",
    "Oligodendrocyte",
    "OPC",
    "Endothelial"
)

all_results <- list()

for(ct in celltypes){

    cat("\n============================\n")
    cat(ct,"\n")
    cat("============================\n")

    idx <- cov$cell_type_std==ct

    x <- logcpm[,idx,drop=FALSE]
    m <- cov[idx,,drop=FALSE]

    cat(
        "AD =",sum(m$condition=="AD"),
        "Control =",sum(m$condition=="Control"),
        "\n"
    )

    if(sum(m$condition=="AD") < 10 ||
       sum(m$condition=="Control") < 10){

        cat("Insufficient nuclei; skipping.\n")
        next
    }

    # Retain genes detected in >=10% nuclei
    detection <- rowMeans(x > 0)

    x <- x[
        detection >= 0.10,
        ,
        drop=FALSE
    ]

    group <- factor(
        m$condition,
        levels=c("Control","AD")
    )

    design <- model.matrix(~0+group)

    colnames(design) <- c(
        "Control",
        "AD"
    )

    fit <- lmFit(
        x,
        design
    )

    cont <- makeContrasts(
        AD-Control,
        levels=design
    )

    fit <- contrasts.fit(
        fit,
        cont
    )

    fit <- eBayes(
        fit,
        robust=TRUE,
        trend=TRUE
    )

    tt <- topTable(
        fit,
        number=Inf,
        adjust.method="BH",
        sort.by="P"
    )

    tt$Gene <- rownames(tt)
    tt$CellType <- ct

    tt <- tt[
        ,
        c(
            "Gene",
            "CellType",
            "logFC",
            "AveExpr",
            "t",
            "P.Value",
            "adj.P.Val",
            "B"
        )
    ]

    outfile <- paste0(
        "results/GSE138852_",
        ct,
        "_AD_vs_Control.tsv"
    )

    write.table(
        tt,
        outfile,
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    all_results[[ct]] <- tt
}

combined <- do.call(
    rbind,
    all_results
)

write.table(
    combined,
    "results/GSE138852_all_celltypes_AD_vs_Control.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# Summary
# ============================================================

summary_table <- do.call(
    rbind,
    lapply(
        names(all_results),
        function(ct){

            z <- all_results[[ct]]

            data.frame(
                CellType=ct,
                Tested=nrow(z),
                FDR_005=sum(
                    z$adj.P.Val < .05,
                    na.rm=TRUE
                ),
                Up=sum(
                    z$adj.P.Val < .05 &
                    z$logFC > 0,
                    na.rm=TRUE
                ),
                Down=sum(
                    z$adj.P.Val < .05 &
                    z$logFC < 0,
                    na.rm=TRUE
                )
            )
        }
    )
)

write.table(
    summary_table,
    "tables/GSE138852_celltype_DE_summary.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n====================================\n")
cat("SINGLE-NUCLEUS VALIDATION COMPLETE\n")
cat("====================================\n\n")

print(summary_table)

cat("\nResults written to results/\n")
