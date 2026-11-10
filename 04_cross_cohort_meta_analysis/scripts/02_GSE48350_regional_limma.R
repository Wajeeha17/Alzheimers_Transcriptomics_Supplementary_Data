suppressPackageStartupMessages(library(limma))

# ============================================================
# Load expression matrix
# ============================================================

expr <- read.delim(
    "data/GSE48350/GSE48350_expression.tsv.gz",
    row.names=1,
    check.names=FALSE
)

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

cat("Expression matrix:",
    nrow(expr),"probes x",
    ncol(expr),"samples\n")

# IMPORTANT:
# GSE48350 processed matrix is used as supplied.
# No additional log2 transformation is applied.

# ============================================================
# Load phenotype
# ============================================================

pheno <- read.delim(
    "metadata/GSE48350_curated_phenotype.tsv",
    stringsAsFactors=FALSE
)

idx <- match(colnames(expr), pheno$sample_id)

if(any(is.na(idx))){
    stop("Expression/phenotype sample mismatch.")
}

pheno <- pheno[idx,]

stopifnot(
    all(colnames(expr) == pheno$sample_id)
)

# ============================================================
# Regions shared with GSE5281
# ============================================================

regions <- c("EC","HIP","SFG")

summary_list <- list()

for(reg in regions){

    cat("\n====================================\n")
    cat("REGION:",reg,"\n")
    cat("====================================\n")

    keep <- pheno$region == reg &
            pheno$status %in% c("Control","AD")

    y <- expr[,keep,drop=FALSE]
    p <- pheno[keep,,drop=FALSE]

    group <- factor(
        p$status,
        levels=c("Control","AD")
    )

    cat("\nSample counts:\n")
    print(table(group))

    # --------------------------------------------------------
    # Design matrix
    # --------------------------------------------------------

    design <- model.matrix(~0 + group)

    colnames(design) <- c(
        "Control",
        "AD"
    )

    # --------------------------------------------------------
    # limma
    # --------------------------------------------------------

    fit <- lmFit(y,design)

    contrast <- makeContrasts(
        AD_vs_Control = AD-Control,
        levels=design
    )

    fit2 <- contrasts.fit(
        fit,
        contrast
    )

    # robust empirical Bayes
    fit2 <- eBayes(
        fit2,
        robust=TRUE
    )

    # --------------------------------------------------------
    # Complete results
    # --------------------------------------------------------

    tt <- topTable(
        fit2,
        coef="AD_vs_Control",
        number=Inf,
        adjust.method="BH",
        sort.by="P"
    )

    tt$probe_id <- rownames(tt)

    tt <- tt[,c(
        "probe_id",
        setdiff(names(tt),"probe_id")
    )]

    # --------------------------------------------------------
    # Standard error for later meta-analysis
    # SE = logFC / moderated t
    # --------------------------------------------------------

    tt$SE <- abs(
        tt$logFC / tt$t
    )

    # handle possible t=0
    tt$SE[!is.finite(tt$SE)] <- NA

    # --------------------------------------------------------
    # DEG classification
    # --------------------------------------------------------

    tt$significance <- "NS"

    tt$significance[
        tt$adj.P.Val < 0.05 &
        tt$logFC >= 0.5
    ] <- "Up"

    tt$significance[
        tt$adj.P.Val < 0.05 &
        tt$logFC <= -0.5
    ] <- "Down"

    # --------------------------------------------------------
    # Save complete table
    # --------------------------------------------------------

    write.table(
        tt,
        paste0(
            "results/GSE48350_",
            reg,
            "_AD_vs_Control_all_probes.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # --------------------------------------------------------
    # Significant DEG table
    # --------------------------------------------------------

    deg <- tt[
        tt$significance != "NS",
        ,
        drop=FALSE
    ]

    write.table(
        deg,
        paste0(
            "results/GSE48350_",
            reg,
            "_AD_vs_Control_DEGs.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # --------------------------------------------------------
    # Counts
    # --------------------------------------------------------

    n_up <- sum(
        tt$significance=="Up",
        na.rm=TRUE
    )

    n_down <- sum(
        tt$significance=="Down",
        na.rm=TRUE
    )

    n_fdr <- sum(
        tt$adj.P.Val < 0.05,
        na.rm=TRUE
    )

    summary_list[[reg]] <- data.frame(
        region=reg,
        samples=ncol(y),
        AD=sum(group=="AD"),
        Control=sum(group=="Control"),
        FDR_005=n_fdr,
        Up=n_up,
        Down=n_down,
        Total_DEG=n_up+n_down
    )

    # --------------------------------------------------------
    # Volcano
    # --------------------------------------------------------

    pdf(
        paste0(
            "figures/GSE48350_",
            reg,
            "_volcano.pdf"
        ),
        width=7,
        height=6
    )

    plot(
        tt$logFC,
        -log10(tt$adj.P.Val),
        pch=20,
        cex=0.45,
        xlab="log2 fold change",
        ylab="-log10(BH-adjusted P-value)",
        main=paste(
            "GSE48350",
            reg,
            "AD vs Control"
        )
    )

    abline(
        v=c(-0.5,0.5),
        lty=2
    )

    abline(
        h=-log10(0.05),
        lty=2
    )

    dev.off()

    # --------------------------------------------------------
    # MA plot
    # --------------------------------------------------------

    pdf(
        paste0(
            "figures/GSE48350_",
            reg,
            "_MA_plot.pdf"
        ),
        width=7,
        height=6
    )

    status <- decideTests(
        fit2,
        p.value=0.05,
        lfc=0.5
    )

    plotMD(
        fit2,
        column=1,
        status=status,
        main=paste(
            "GSE48350",
            reg,
            "AD vs Control"
        )
    )

    abline(
        h=c(-0.5,0.5),
        lty=2
    )

    dev.off()
}

# ============================================================
# Combined summary
# ============================================================

summary_table <- do.call(
    rbind,
    summary_list
)

write.table(
    summary_table,
    "tables/GSE48350_regional_DEG_summary.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n\n====================================\n")
cat("FINAL SUMMARY\n")
cat("====================================\n")

print(
    summary_table,
    row.names=FALSE
)

cat("\nGSE48350 regional DEA completed.\n")
