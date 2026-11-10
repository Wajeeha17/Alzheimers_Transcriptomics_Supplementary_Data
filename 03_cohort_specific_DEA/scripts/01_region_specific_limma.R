suppressPackageStartupMessages(library(limma))

dir.create("results", showWarnings=FALSE)
dir.create("tables", showWarnings=FALSE)
dir.create("figures", showWarnings=FALSE)

# ------------------------------------------------------------
# 1. Load expression matrix
# ------------------------------------------------------------
expr <- read.delim(
    "input/GSE5281_expression.tsv.gz",
    row.names=1,
    check.names=FALSE
)

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

cat("Original matrix:", nrow(expr), "features x",
    ncol(expr), "samples\n")

# GSE5281 values are linear-scale MAS5 intensities.
# Apply log2 transformation only.
expr[expr <= 0] <- NA
expr <- log2(expr)

# Remove probes containing non-finite values after transformation
keep <- apply(expr,1,function(x) all(is.finite(x)))
expr <- expr[keep,,drop=FALSE]

cat("After filtering:",nrow(expr),"features\n")

# ------------------------------------------------------------
# 2. Load phenotype
# ------------------------------------------------------------
pheno <- read.delim(
    "metadata/GSE5281_curated_phenotype.tsv",
    stringsAsFactors=FALSE,
    check.names=FALSE
)

# Match phenotype to expression columns
idx <- match(colnames(expr), pheno$sample_id)

if(any(is.na(idx))){
    stop(
        "Expression samples missing from phenotype: ",
        paste(colnames(expr)[is.na(idx)],collapse=", ")
    )
}

pheno <- pheno[idx,,drop=FALSE]

if(!identical(colnames(expr),pheno$sample_id))
    stop("Sample ordering error")

cat("Matched samples:",nrow(pheno),"\n")

# ------------------------------------------------------------
# 3. Regions
# ------------------------------------------------------------
regions <- c("EC","HIP","MTG","PC","SFG","VCX")

summary_list <- list()

for(region_name in regions){

    cat("\n====================================\n")
    cat("REGION:",region_name,"\n")
    cat("====================================\n")

    sel <- pheno$region == region_name

    x <- expr[,sel,drop=FALSE]
    p <- pheno[sel,,drop=FALSE]

    group <- factor(
        p$status,
        levels=c("Control","AD")
    )

    cat("Samples:",ncol(x),"\n")
    cat("Control:",sum(group=="Control"),"\n")
    cat("AD:",sum(group=="AD"),"\n")

    # Design matrix
    design <- model.matrix(~0 + group)
    colnames(design) <- c("Control","AD")

    # Fit limma model
    fit <- lmFit(x,design)

    contrast <- makeContrasts(
        AD_vs_Control = AD-Control,
        levels=design
    )

    fit2 <- contrasts.fit(fit,contrast)
    fit2 <- eBayes(fit2)

    # Complete ranked results
    res <- topTable(
        fit2,
        coef="AD_vs_Control",
        number=Inf,
        adjust.method="BH",
        sort.by="P"
    )

    res$probe_id <- rownames(res)

    # Put probe ID first
    res <- res[
        ,
        c(
            "probe_id",
            setdiff(names(res),"probe_id")
        )
    ]

    # Classification
    res$significance <- "NS"

    res$significance[
        res$adj.P.Val < 0.05 &
        res$logFC >= 0.5
    ] <- "Up"

    res$significance[
        res$adj.P.Val < 0.05 &
        res$logFC <= -0.5
    ] <- "Down"

    # Save complete results
    write.table(
        res,
        paste0(
            "results/",
            region_name,
            "_AD_vs_Control_all_probes.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    # Significant DEG table
    deg <- res[
        res$significance %in% c("Up","Down"),
        ,
        drop=FALSE
    ]

    write.table(
        deg,
        paste0(
            "results/",
            region_name,
            "_AD_vs_Control_DEGs.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    n_up <- sum(res$significance=="Up")
    n_down <- sum(res$significance=="Down")
    n_fdr <- sum(res$adj.P.Val < 0.05,na.rm=TRUE)

    summary_list[[region_name]] <- data.frame(
        region=region_name,
        samples=ncol(x),
        AD=sum(group=="AD"),
        Control=sum(group=="Control"),
        FDR_005=n_fdr,
        Up=n_up,
        Down=n_down,
        Total_DEG=n_up+n_down
    )

    cat("FDR < 0.05:",n_fdr,"\n")
    cat("Up:",n_up,"\n")
    cat("Down:",n_down,"\n")

    # --------------------------------------------------------
    # Volcano plot
    # --------------------------------------------------------
    pdf(
        paste0(
            "figures/",
            region_name,
            "_volcano.pdf"
        ),
        width=7,
        height=6
    )

    y <- -log10(pmax(res$adj.P.Val,1e-300))

    plot(
        res$logFC,
        y,
        pch=16,
        cex=0.45,
        xlab="log2 fold change (AD - Control)",
        ylab="-log10 adjusted P-value",
        main=paste(region_name,"AD vs Control")
    )

    abline(v=c(-0.5,0.5),lty=2)
    abline(h=-log10(0.05),lty=2)

    dev.off()

    # --------------------------------------------------------
    # MA plot
    # --------------------------------------------------------
    pdf(
        paste0(
            "figures/",
            region_name,
            "_MA_plot.pdf"
        ),
        width=7,
        height=6
    )

    plotMD(
        fit2,
        column=1,
        status=decideTests(
            fit2,
            method="separate",
            adjust.method="BH",
            p.value=0.05,
            lfc=0.5
        ),
        main=paste(region_name,"AD vs Control")
    )

    abline(h=c(-0.5,0.5),lty=2)

    dev.off()

    rm(x,p,fit,fit2,res,deg)
    gc()
}

# ------------------------------------------------------------
# 4. Regional DEG summary
# ------------------------------------------------------------
summary_table <- do.call(rbind,summary_list)
rownames(summary_table) <- NULL

write.table(
    summary_table,
    "tables/GSE5281_regional_DEG_summary.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n====================================\n")
cat("FINAL REGIONAL SUMMARY\n")
cat("====================================\n")

print(summary_table)

cat("\nAnalysis completed successfully.\n")
