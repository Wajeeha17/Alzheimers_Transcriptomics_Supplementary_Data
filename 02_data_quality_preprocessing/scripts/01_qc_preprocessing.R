suppressPackageStartupMessages({
    library(limma)
})

files <- list.files(
    "input/expression_data",
    pattern="expression\\.tsv\\.gz$",
    full.names=TRUE
)

dir.create("results", showWarnings=FALSE)
dir.create("figures", showWarnings=FALSE)
dir.create("tables", showWarnings=FALSE)

qc_summary <- list()

for(f in files){

    dataset <- sub("_platform.*", "", basename(f))

    cat("\n==============================\n")
    cat("DATASET:", dataset, "\n")
    cat("==============================\n")

    x <- read.delim(
        gzfile(f),
        row.names=1,
        check.names=FALSE
    )

    x <- as.matrix(x)
    storage.mode(x) <- "numeric"

    cat("Features:", nrow(x), "\n")
    cat("Samples :", ncol(x), "\n")

    missing_n <- sum(is.na(x))
    nonfinite_n <- sum(!is.finite(x), na.rm=TRUE)

    q <- quantile(
        x,
        probs=c(0,0.25,0.5,0.75,0.99,1),
        na.rm=TRUE
    )

    cat("Quantiles:\n")
    print(q)

    # GEO/limma-style heuristic for determining whether
    # log2 transformation is likely required
    log2_needed <- (q[5] > 100) ||
                   ((q[6] - q[1] > 50) && q[2] > 0)

    if(log2_needed){

        cat("Applying log2 transformation\n")

        x[x <= 0] <- NA
        x <- log2(x)

    } else {

        cat("Expression appears already log-scaled\n")
    }

    # Remove features with no usable measurements
    keep <- rowSums(is.finite(x)) > 0
    x <- x[keep,,drop=FALSE]

    # Median imputation only for residual missing values
    if(anyNA(x)){

        for(i in seq_len(nrow(x))){

            bad <- !is.finite(x[i,])

            if(any(bad)){

                med <- median(x[i,!bad], na.rm=TRUE)

                if(is.finite(med))
                    x[i,bad] <- med
            }
        }
    }

    # Remove zero-variance features
    v <- apply(x,1,var,na.rm=TRUE)

    x <- x[
        is.finite(v) & v > 0,
        ,
        drop=FALSE
    ]

    cat("Features after filtering:", nrow(x), "\n")

    # Save filtered/log-scaled matrix BEFORE optional normalization
    pre_norm_file <- paste0(
        "results/",
        dataset,
        "_filtered_expression.tsv.gz"
    )

    con <- gzfile(pre_norm_file,"wt")

    write.table(
        x,
        con,
        sep="\t",
        quote=FALSE,
        col.names=NA
    )

    close(con)

    # Boxplot before normalization
    pdf(
        paste0("figures/",dataset,"_boxplot_before_normalization.pdf"),
        width=10,
        height=6
    )

    boxplot(
        x,
        outline=FALSE,
        las=2,
        main=paste(dataset,"Before normalization"),
        ylab="Expression"
    )

    dev.off()

    # Quantile normalize within dataset
    xn <- normalizeBetweenArrays(
        x,
        method="quantile"
    )

    # Boxplot after normalization
    pdf(
        paste0("figures/",dataset,"_boxplot_after_normalization.pdf"),
        width=10,
        height=6
    )

    boxplot(
        xn,
        outline=FALSE,
        las=2,
        main=paste(dataset,"After normalization"),
        ylab="Normalized expression"
    )

    dev.off()

    # PCA using most variable features
    vars <- apply(xn,1,var)

    n_top <- min(5000,length(vars))

    top <- order(
        vars,
        decreasing=TRUE
    )[seq_len(n_top)]

    pca <- prcomp(
        t(xn[top,,drop=FALSE]),
        scale.=FALSE
    )

    pca_df <- data.frame(
        sample=rownames(pca$x),
        PC1=pca$x[,1],
        PC2=pca$x[,2]
    )

    write.table(
        pca_df,
        paste0("tables/",dataset,"_PCA_coordinates.tsv"),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    pdf(
        paste0("figures/",dataset,"_PCA.pdf"),
        width=7,
        height=6
    )

    plot(
        pca$x[,1],
        pca$x[,2],
        pch=19,
        xlab="PC1",
        ylab="PC2",
        main=paste(dataset,"PCA")
    )

    text(
        pca$x[,1],
        pca$x[,2],
        labels=rownames(pca$x),
        cex=0.45,
        pos=3
    )

    dev.off()

    # Sample correlation
    cor_mat <- cor(
        xn[top,,drop=FALSE],
        method="pearson",
        use="pairwise.complete.obs"
    )

    write.table(
        cor_mat,
        paste0("tables/",dataset,"_sample_correlations.tsv"),
        sep="\t",
        quote=FALSE,
        col.names=NA
    )

    # Hierarchical clustering
    hc <- hclust(
        as.dist(1-cor_mat),
        method="average"
    )

    pdf(
        paste0("figures/",dataset,"_sample_clustering.pdf"),
        width=12,
        height=7
    )

    plot(
        hc,
        main=paste(dataset,"Sample clustering"),
        xlab="",
        sub="",
        cex=0.5
    )

    dev.off()

    # Save normalized matrix
    outfile <- paste0(
        "results/",
        dataset,
        "_normalized_expression.tsv.gz"
    )

    con <- gzfile(outfile,"wt")

    write.table(
        xn,
        con,
        sep="\t",
        quote=FALSE,
        col.names=NA
    )

    close(con)

    qc_summary[[dataset]] <- data.frame(
        dataset=dataset,
        original_features=nrow(x),
        samples=ncol(x),
        missing_original=missing_n,
        nonfinite_original=nonfinite_n,
        log2_applied=log2_needed,
        normalized_features=nrow(xn)
    )

    rm(x,xn,pca,cor_mat)
    gc()
}

if(length(qc_summary)>0){

    out <- do.call(rbind,qc_summary)

    write.table(
        out,
        "tables/QC_summary.tsv",
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    cat("\n===== QC SUMMARY =====\n")
    print(out)
}
