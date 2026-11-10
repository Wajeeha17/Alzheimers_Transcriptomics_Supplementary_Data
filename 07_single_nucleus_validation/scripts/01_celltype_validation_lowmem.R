suppressPackageStartupMessages({
    library(limma)
})

count_file <- "raw/GSE138852_counts.csv.gz"
cov_file   <- "metadata/GSE138852_covariates.csv.gz"

dir.create("results", showWarnings=FALSE)
dir.create("tables", showWarnings=FALSE)
dir.create("logs", showWarnings=FALSE)

# ------------------------------------------------------------
# METADATA
# ------------------------------------------------------------

cov <- read.csv(
    gzfile(cov_file),
    check.names=FALSE,
    stringsAsFactors=FALSE
)

names(cov)[1] <- "cell_id"
names(cov)[names(cov)=="oupSample.batchCond"] <- "condition"
names(cov)[names(cov)=="oupSample.cellType"] <- "cell_type"

cov$condition[cov$condition=="ct"] <- "Control"

cellmap <- c(
    neuron="Neuron",
    astro="Astrocyte",
    mg="Microglia",
    oligo="Oligodendrocyte",
    OPC="OPC",
    endo="Endothelial"
)

cov$cell_type_std <- unname(
    cellmap[cov$cell_type]
)

cov <- cov[
    !is.na(cov$cell_type_std) &
    cov$condition %in% c("AD","Control"),
    ,
    drop=FALSE
]

cat("\nNuclei retained:",nrow(cov),"\n")

cat("\nCell counts:\n")
print(
    addmargins(
        table(
            cov$cell_type_std,
            cov$condition
        )
    )
)

# ------------------------------------------------------------
# READ COUNTS ONCE
# ------------------------------------------------------------

cat("\nReading count matrix...\n")

dat <- read.csv(
    gzfile(count_file),
    check.names=FALSE,
    stringsAsFactors=FALSE
)

genes <- dat[[1]]

counts <- as.matrix(
    dat[,-1,drop=FALSE]
)

storage.mode(counts) <- "numeric"
rownames(counts) <- genes

rm(dat)
gc()

cat(
    "\nMatrix:",
    nrow(counts),"genes x",
    ncol(counts),"nuclei\n"
)

# ------------------------------------------------------------
# Match cells
# ------------------------------------------------------------

common <- intersect(
    colnames(counts),
    cov$cell_id
)

counts <- counts[
    ,
    common,
    drop=FALSE
]

cov <- cov[
    match(common,cov$cell_id),
    ,
    drop=FALSE
]

stopifnot(
    identical(
        colnames(counts),
        cov$cell_id
    )
)

cat("Matched nuclei:",ncol(counts),"\n")

# ------------------------------------------------------------
# CELL-TYPE ANALYSIS
# Process ONE cell type at a time
# ------------------------------------------------------------

celltypes <- c(
    "Neuron",
    "Astrocyte",
    "Microglia",
    "Oligodendrocyte",
    "OPC",
    "Endothelial"
)

summary_list <- list()

for(ct in celltypes){

    cat("\n====================================\n")
    cat("Processing:",ct,"\n")
    cat("====================================\n")

    idx <- which(
        cov$cell_type_std==ct
    )

    meta <- cov[
        idx,
        ,
        drop=FALSE
    ]

    cat(
        "AD:",sum(meta$condition=="AD"),
        "| Control:",sum(meta$condition=="Control"),
        "| Total:",nrow(meta),"\n"
    )

    if(sum(meta$condition=="AD") < 10 ||
       sum(meta$condition=="Control") < 10){

        cat("Insufficient nuclei. Skipping.\n")
        next
    }

    # Only current cell type enters memory
    x <- counts[
        ,
        idx,
        drop=FALSE
    ]

    # Detection filter: >=10% nuclei
    detected_fraction <- rowMeans(x > 0)

    x <- x[
        detected_fraction >= 0.10,
        ,
        drop=FALSE
    ]

    rm(detected_fraction)
    gc()

    cat(
        "Genes retained:",
        nrow(x),"\n"
    )

    # ---------------------------------------------
    # Normalize only current cell type
    # ---------------------------------------------

    lib <- colSums(x)
    lib[lib <= 0] <- 1

    # modify matrix column-by-column to reduce
    # temporary memory usage
    for(j in seq_len(ncol(x))){
        x[,j] <- log2(
            (x[,j] / lib[j]) * 1e6 + 1
        )
    }

    rm(lib)
    gc()

    # ---------------------------------------------
    # limma
    # ---------------------------------------------

    group <- factor(
        meta$condition,
        levels=c("Control","AD")
    )

    design <- model.matrix(
        ~0 + group
    )

    colnames(design) <- c(
        "Control",
        "AD"
    )

    fit <- lmFit(
        x,
        design
    )

    contrast <- makeContrasts(
        AD-Control,
        levels=design
    )

    fit <- contrasts.fit(
        fit,
        contrast
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

    summary_list[[ct]] <- data.frame(
        CellType=ct,
        AD_n=sum(meta$condition=="AD"),
        Control_n=sum(meta$condition=="Control"),
        Tested=nrow(tt),

        FDR_005=sum(
            tt$adj.P.Val < 0.05,
            na.rm=TRUE
        ),

        Up=sum(
            tt$adj.P.Val < 0.05 &
            tt$logFC > 0,
            na.rm=TRUE
        ),

        Down=sum(
            tt$adj.P.Val < 0.05 &
            tt$logFC < 0,
            na.rm=TRUE
        )
    )

    cat(
        "FDR < 0.05:",
        summary_list[[ct]]$FDR_005,
        "\n"
    )

    rm(
        x,
        fit,
        tt,
        design,
        contrast,
        group,
        meta
    )

    gc()
}

# ------------------------------------------------------------
# SUMMARY
# ------------------------------------------------------------

if(length(summary_list)==0){
    stop("No cell types completed.")
}

summary_table <- do.call(
    rbind,
    summary_list
)

rownames(summary_table) <- NULL

write.table(
    summary_table,
    "tables/GSE138852_celltype_DE_summary.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n====================================\n")
cat("VALIDATION COMPLETED SUCCESSFULLY\n")
cat("====================================\n\n")

print(summary_table)

cat(
    "\nSaved:\n",
    "tables/GSE138852_celltype_DE_summary.tsv\n"
)
