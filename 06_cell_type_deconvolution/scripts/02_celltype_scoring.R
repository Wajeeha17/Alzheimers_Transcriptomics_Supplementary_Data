# ============================================================
# GSE5281 brain cell-type enrichment scoring
# ============================================================

input_file <- "results/GSE5281_gene_expression_log2.tsv.gz"

if(!file.exists(input_file)){
    stop(
        "Gene-level expression matrix missing: ",
        input_file,
        "\nRun 01_probe_to_gene_expression.R first."
    )
}

x <- read.delim(
    input_file,
    row.names=1,
    check.names=FALSE
)

x <- as.matrix(x)
storage.mode(x) <- "numeric"

cat("Gene expression matrix:",
    nrow(x),"genes x",
    ncol(x),"samples\n\n")

# ============================================================
# Brain cell-type marker panels
#
# These are marker-based enrichment scores, NOT absolute
# cell fractions.
# ============================================================

markers <- list(

    Neuron = c(
        "RBFOX3","SNAP25","SYT1","SYN1","SYP",
        "MAP2","NEFL","NEFM","NEFH","GRIN1",
        "CAMK2A","STMN2","TUBB3","ENO2"
    ),

    Astrocyte = c(
        "GFAP","AQP4","ALDH1L1","SLC1A2",
        "SLC1A3","GLUL","GJA1","S100B",
        "SOX9","CLU","APOE","Aldoc"
    ),

    Microglia = c(
        "AIF1","CX3CR1","P2RY12","TMEM119",
        "CSF1R","TYROBP","TREM2","C1QA",
        "C1QB","C1QC","CTSS","FCER1G",
        "LAPTM5"
    ),

    Oligodendrocyte = c(
        "MBP","PLP1","MOG","MAG","MOBP",
        "CNP","CLDN11","OLIG1","OLIG2",
        "SOX10","ERMN","MYRF"
    ),

    OPC = c(
        "PDGFRA","CSPG4","VCAN","SOX10",
        "OLIG1","OLIG2","NKX2-2","PTPRZ1",
        "BCAN"
    ),

    Endothelial = c(
        "CLDN5","PECAM1","VWF","KDR",
        "FLT1","EMCN","RAMP2","ESAM",
        "CDH5","ENG","ERG"
    )
)

# Correct capitalization
markers$Astrocyte[
    markers$Astrocyte=="Aldoc"
] <- "ALDOC"

# ============================================================
# Match markers to expression matrix
# ============================================================

marker_summary <- data.frame()

matched_markers <- list()

for(ct in names(markers)){

    requested <- unique(markers[[ct]])

    matched <- intersect(
        requested,
        rownames(x)
    )

    matched_markers[[ct]] <- matched

    marker_summary <- rbind(
        marker_summary,
        data.frame(
            CellType=ct,
            Requested=length(requested),
            Matched=length(matched),
            Coverage=length(matched)/length(requested)
        )
    )

    cat(
        ct,":",
        length(matched),"/",
        length(requested),
        "markers matched\n"
    )

    cat(
        paste(matched,collapse=", "),
        "\n\n"
    )
}

write.table(
    marker_summary,
    "tables/GSE5281_celltype_marker_coverage.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# Z-score each gene across samples
#
# This prevents highly expressed genes from dominating the
# marker score.
# ============================================================

gene_z <- t(
    scale(
        t(x),
        center=TRUE,
        scale=TRUE
    )
)

gene_z[!is.finite(gene_z)] <- NA

# ============================================================
# Calculate sample-level enrichment score
# Mean standardized expression of markers
# ============================================================

score_matrix <- matrix(
    NA_real_,
    nrow=ncol(x),
    ncol=length(markers)
)

rownames(score_matrix) <- colnames(x)
colnames(score_matrix) <- names(markers)

for(ct in names(markers)){

    genes <- matched_markers[[ct]]

    if(length(genes)<3){

        warning(
            ct,
            " has fewer than 3 matched markers."
        )

        next
    }

    score_matrix[,ct] <- colMeans(
        gene_z[
            genes,
            ,
            drop=FALSE
        ],
        na.rm=TRUE
    )
}

# ============================================================
# Output
# ============================================================

scores <- data.frame(
    sample_id=rownames(score_matrix),
    score_matrix,
    check.names=FALSE
)

write.table(
    scores,
    "results/GSE5281_celltype_scores.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n====================================\n")
cat("CELL-TYPE SCORING COMPLETED\n")
cat("====================================\n")

cat("\nSamples:",nrow(scores),"\n")

cat("\nScore summary:\n")
print(summary(scores[,-1]))

cat("\nOutput:\n")
cat("results/GSE5281_celltype_scores.tsv\n")
