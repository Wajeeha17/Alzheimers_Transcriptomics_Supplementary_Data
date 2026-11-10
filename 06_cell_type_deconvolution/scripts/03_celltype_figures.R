suppressPackageStartupMessages({
    library(ggplot2)
})

# ============================================================
# INPUT
# ============================================================

score_file <- "results/GSE5281_celltype_scores.tsv"
pheno_file <- "../02_data_acquisition_metadata/raw_metadata/GEO/GSE5281_platform1_metadata.tsv"

if(!file.exists(score_file)){
    stop(
        "\nCell-type score file does not exist:\n",
        score_file,
        "\n\nRun the cell-type scoring/deconvolution step first."
    )
}

dir.create("figures/publication",
           recursive=TRUE,
           showWarnings=FALSE)

scores <- read.delim(
    score_file,
    check.names=FALSE,
    stringsAsFactors=FALSE
)

cat("Score table dimensions:",
    nrow(scores),"x",ncol(scores),"\n")

cat("\nColumns:\n")
print(names(scores))

# ============================================================
# Detect sample-ID column
# ============================================================

possible_ids <- c(
    "sample_id",
    "Sample",
    "sample",
    "GSM",
    "geo_accession"
)

id_col <- possible_ids[
    possible_ids %in% names(scores)
][1]

if(is.na(id_col)){
    stop("Could not identify sample-ID column.")
}

# ============================================================
# Get phenotype from previously curated GSE5281 metadata
# ============================================================

possible_pheno <- c(
    "../04_cohort_specific_dea/metadata/GSE5281_curated_phenotype.tsv",
    "../03_qc_and_preprocessing/metadata/GSE5281_curated_phenotype.tsv"
)

pf <- possible_pheno[file.exists(possible_pheno)][1]

if(is.na(pf)){
    stop("Could not find curated GSE5281 phenotype file.")
}

pheno <- read.delim(
    pf,
    stringsAsFactors=FALSE
)

idx <- match(scores[[id_col]],pheno$sample_id)

if(any(is.na(idx))){
    stop("Some cell-score samples cannot be matched to phenotype.")
}

scores$region <- pheno$region[idx]
scores$status <- pheno$status[idx]

# ============================================================
# Identify cell-type columns
# ============================================================

expected <- c(
    "Neuron",
    "Astrocyte",
    "Microglia",
    "Oligodendrocyte",
    "OPC",
    "Endothelial"
)

celltypes <- expected[
    expected %in% names(scores)
]

if(length(celltypes)<2){
    stop(
        "Expected cell-type columns were not found.\nAvailable columns:\n",
        paste(names(scores),collapse=", ")
    )
}

cat("\nDetected cell types:\n")
print(celltypes)

# ============================================================
# Convert wide -> long without extra packages
# ============================================================

long <- do.call(
    rbind,
    lapply(celltypes,function(ct){

        data.frame(
            sample_id=scores[[id_col]],
            region=scores$region,
            status=scores$status,
            CellType=ct,
            Score=as.numeric(scores[[ct]]),
            stringsAsFactors=FALSE
        )
    })
)

long$region <- factor(
    long$region,
    levels=c("EC","HIP","MTG","PC","SFG","VCX")
)

long$status <- factor(
    long$status,
    levels=c("Control","AD")
)

# ============================================================
# FIGURE 1
# AD vs Control boxplots across regions
# ============================================================

p1 <- ggplot(
    long,
    aes(
        x=status,
        y=Score,
        fill=status
    )
) +
    geom_boxplot(
        width=.65,
        outlier.shape=NA,
        alpha=.75
    ) +
    geom_jitter(
        width=.12,
        size=.65,
        alpha=.45
    ) +
    facet_grid(
        CellType ~ region,
        scales="free_y"
    ) +
    labs(
        title="Cell-type enrichment across Alzheimer’s disease brain regions",
        x=NULL,
        y="Relative cell-type enrichment score"
    ) +
    theme_bw(base_size=10) +
    theme(
        legend.position="none",
        strip.text=element_text(face="bold"),
        axis.text.x=element_text(
            angle=45,
            hjust=1
        ),
        plot.title=element_text(
            face="bold",
            hjust=.5
        )
    )

ggsave(
    "figures/publication/Figure_celltype_boxplots.pdf",
    p1,
    width=14,
    height=12
)

ggsave(
    "figures/publication/Figure_celltype_boxplots.tiff",
    p1,
    width=14,
    height=12,
    dpi=600,
    compression="lzw"
)

# ============================================================
# Calculate regional AD-Control effect sizes
# Hedges' g
# ============================================================

hedges_g <- function(ad,ctrl){

    ad <- ad[is.finite(ad)]
    ctrl <- ctrl[is.finite(ctrl)]

    n1 <- length(ad)
    n0 <- length(ctrl)

    if(n1<3 || n0<3)
        return(NA_real_)

    s_pool <- sqrt(
        ((n1-1)*var(ad) +
         (n0-1)*var(ctrl)) /
        (n1+n0-2)
    )

    if(!is.finite(s_pool) || s_pool==0)
        return(NA_real_)

    d <- (mean(ad)-mean(ctrl))/s_pool

    J <- 1 -
        (3/(4*(n1+n0)-9))

    J*d
}

regions <- c(
    "EC","HIP","MTG",
    "PC","SFG","VCX"
)

stat_list <- list()

k <- 1

for(ct in celltypes){

    for(reg in regions){

        z <- long[
            long$CellType==ct &
            long$region==reg,
        ]

        ad <- z$Score[z$status=="AD"]
        ctrl <- z$Score[z$status=="Control"]

        wt <- tryCatch(
            wilcox.test(
                ad,
                ctrl,
                exact=FALSE
            ),
            error=function(e) NULL
        )

        stat_list[[k]] <- data.frame(
            CellType=ct,
            Region=reg,
            AD_n=sum(is.finite(ad)),
            Control_n=sum(is.finite(ctrl)),
            AD_median=median(ad,na.rm=TRUE),
            Control_median=median(ctrl,na.rm=TRUE),
            Hedges_g=hedges_g(ad,ctrl),
            P_value=if(is.null(wt))
                NA_real_
            else
                wt$p.value
        )

        k <- k+1
    }
}

stats <- do.call(
    rbind,
    stat_list
)

stats$FDR <- p.adjust(
    stats$P_value,
    method="BH"
)

stats$Significance <- ""

stats$Significance[
    stats$FDR < .05
] <- "*"

stats$Significance[
    stats$FDR < .01
] <- "**"

stats$Significance[
    stats$FDR < .001
] <- "***"

write.table(
    stats,
    "tables/GSE5281_celltype_regional_statistics.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# FIGURE 2
# Effect-size heatmap
# ============================================================

stats$Region <- factor(
    stats$Region,
    levels=regions
)

stats$CellType <- factor(
    stats$CellType,
    levels=rev(celltypes)
)

p2 <- ggplot(
    stats,
    aes(
        x=Region,
        y=CellType,
        fill=Hedges_g
    )
) +
    geom_tile(
        colour="white",
        linewidth=.7
    ) +
    geom_text(
        aes(label=Significance),
        size=5,
        fontface="bold"
    ) +
    scale_fill_gradient2(
        midpoint=0,
        name="Hedges' g"
    ) +
    labs(
        title="Regional cell-type alterations in Alzheimer’s disease",
        subtitle="Positive values indicate higher cell-type enrichment scores in AD",
        x="Brain region",
        y=NULL
    ) +
    theme_minimal(base_size=12) +
    theme(
        panel.grid=element_blank(),
        axis.text=element_text(face="bold"),
        plot.title=element_text(
            face="bold",
            hjust=.5
        ),
        plot.subtitle=element_text(
            hjust=.5
        )
    )

ggsave(
    "figures/publication/Figure_celltype_effect_heatmap.pdf",
    p2,
    width=8,
    height=6
)

ggsave(
    "figures/publication/Figure_celltype_effect_heatmap.tiff",
    p2,
    width=8,
    height=6,
    dpi=600,
    compression="lzw"
)

# ============================================================
# FIGURE 3
# Mean cell-type score heatmap by region/status
# ============================================================

means <- aggregate(
    Score ~ region + status + CellType,
    data=long,
    FUN=mean
)

# Z-score within each cell type
means$Z <- ave(
    means$Score,
    means$CellType,
    FUN=function(x){
        if(sd(x)==0)
            rep(0,length(x))
        else
            as.numeric(scale(x))
    }
)

means$Group <- paste(
    means$region,
    means$status,
    sep="_"
)

p3 <- ggplot(
    means,
    aes(
        x=Group,
        y=CellType,
        fill=Z
    )
) +
    geom_tile(
        colour="white",
        linewidth=.5
    ) +
    scale_fill_gradient2(
        midpoint=0,
        name="Z-score"
    ) +
    labs(
        title="Cell-type enrichment profile across brain regions",
        x=NULL,
        y=NULL
    ) +
    theme_minimal(base_size=11) +
    theme(
        panel.grid=element_blank(),
        axis.text.x=element_text(
            angle=45,
            hjust=1
        ),
        axis.text.y=element_text(
            face="bold"
        ),
        plot.title=element_text(
            face="bold",
            hjust=.5
        )
    )

ggsave(
    "figures/publication/Figure_celltype_score_heatmap.pdf",
    p3,
    width=12,
    height=6
)

ggsave(
    "figures/publication/Figure_celltype_score_heatmap.tiff",
    p3,
    width=12,
    height=6,
    dpi=600,
    compression="lzw"
)

cat("\n====================================\n")
cat("CELL-TYPE FIGURES COMPLETED\n")
cat("====================================\n")

cat("\nStatistical results:\n")
print(stats)

cat("\nFigures saved in:\n")
cat("figures/publication/\n")
