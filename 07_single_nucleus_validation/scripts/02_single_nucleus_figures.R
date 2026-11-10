suppressPackageStartupMessages({
    library(ggplot2)
})

dir.create("figures/publication", recursive=TRUE, showWarnings=FALSE)

celltypes <- c(
    "Neuron",
    "Astrocyte",
    "Microglia",
    "Oligodendrocyte",
    "OPC",
    "Endothelial"
)

# ============================================================
# READ ALL CELL-TYPE RESULTS
# ============================================================

res <- list()

for(ct in celltypes){

    f <- paste0(
        "results/GSE138852_",
        ct,
        "_AD_vs_Control.tsv"
    )

    if(!file.exists(f)){
        warning("Missing: ",f)
        next
    }

    z <- read.delim(
        f,
        stringsAsFactors=FALSE
    )

    z$CellType <- ct
    res[[ct]] <- z
}

if(length(res)==0)
    stop("No single-nucleus result files found.")

allres <- do.call(rbind,res)
rownames(allres) <- NULL

allres$minusLog10FDR <- -log10(
    pmax(allres$adj.P.Val,1e-300)
)

# ============================================================
# FIGURE 1 — DEG BURDEN
# ============================================================

summary_file <- "tables/GSE138852_celltype_DE_summary.tsv"

if(!file.exists(summary_file))
    stop("Summary table missing.")

s <- read.delim(
    summary_file,
    stringsAsFactors=FALSE
)

s$CellType <- factor(
    s$CellType,
    levels=celltypes
)

burden <- rbind(
    data.frame(
        CellType=s$CellType,
        Direction="Up",
        Count=s$Up
    ),
    data.frame(
        CellType=s$CellType,
        Direction="Down",
        Count=s$Down
    )
)

p1 <- ggplot(
    burden,
    aes(
        x=CellType,
        y=Count,
        fill=Direction
    )
) +
    geom_col(
        position="dodge",
        width=.72
    ) +
    labs(
        title="Cell-type-specific transcriptional alterations in Alzheimer's disease",
        subtitle="GSE138852 entorhinal cortex single-nucleus transcriptomics",
        x=NULL,
        y="Genes with FDR < 0.05"
    ) +
    theme_classic(base_size=12) +
    theme(
        axis.text.x=element_text(
            angle=35,
            hjust=1
        ),
        plot.title=element_text(
            face="bold"
        ),
        legend.title=element_blank()
    )

ggsave(
    "figures/publication/Figure_SN01_DEG_burden.pdf",
    p1,
    width=8.5,
    height=6
)

ggsave(
    "figures/publication/Figure_SN01_DEG_burden.tiff",
    p1,
    width=8.5,
    height=6,
    dpi=600,
    compression="lzw"
)

# ============================================================
# FIGURE 2 — TOP GENES ACROSS CELL TYPES
#
# Select strongest genes independently within each cell type
# based on absolute logFC among FDR-significant genes.
# ============================================================

top_list <- list()

for(ct in names(res)){

    z <- res[[ct]]

    sig <- z[
        is.finite(z$logFC) &
        is.finite(z$adj.P.Val) &
        z$adj.P.Val < 0.05,
        ,
        drop=FALSE
    ]

    if(nrow(sig)==0)
        next

    sig <- sig[
        order(-abs(sig$logFC)),
        ,
        drop=FALSE
    ]

    top_list[[ct]] <- head(sig,10)
}

if(length(top_list)>0){

    top <- do.call(rbind,top_list)

    genes <- unique(top$Gene)

    heat <- allres[
        allres$Gene %in% genes,
        c(
            "Gene",
            "CellType",
            "logFC",
            "adj.P.Val"
        )
    ]

    # Keep only one row per gene/cell type
    heat <- heat[
        !duplicated(
            paste(
                heat$Gene,
                heat$CellType
            )
        ),
    ]

    heat$Significance <- ""

    heat$Significance[
        heat$adj.P.Val < .05
    ] <- "*"

    heat$Significance[
        heat$adj.P.Val < .01
    ] <- "**"

    heat$Significance[
        heat$adj.P.Val < .001
    ] <- "***"

    heat$CellType <- factor(
        heat$CellType,
        levels=celltypes
    )

    # Gene order based on maximum absolute effect
    gene_order <- aggregate(
        abs(logFC) ~ Gene,
        data=heat,
        FUN=max
    )

    gene_order <- gene_order[
        order(
            gene_order[[2]]
        ),
    ]

    heat$Gene <- factor(
        heat$Gene,
        levels=gene_order$Gene
    )

    p2 <- ggplot(
        heat,
        aes(
            x=CellType,
            y=Gene,
            fill=logFC
        )
    ) +
        geom_tile(
            colour="white",
            linewidth=.35
        ) +
        geom_text(
            aes(label=Significance),
            size=2.7
        ) +
        scale_fill_gradient2(
            midpoint=0,
            name="log2FC"
        ) +
        labs(
            title="Cell-type localization of AD-associated transcriptional changes",
            subtitle="Top cell-type-specific genes in GSE138852",
            x=NULL,
            y=NULL
        ) +
        theme_minimal(base_size=10) +
        theme(
            panel.grid=element_blank(),
            axis.text.x=element_text(
                angle=35,
                hjust=1,
                face="bold"
            ),
            axis.text.y=element_text(
                size=7
            ),
            plot.title=element_text(
                face="bold"
            )
        )

    ggsave(
        "figures/publication/Figure_SN02_celltype_heatmap.pdf",
        p2,
        width=9,
        height=max(8,length(genes)*0.20)
    )

    ggsave(
        "figures/publication/Figure_SN02_celltype_heatmap.tiff",
        p2,
        width=9,
        height=max(8,length(genes)*0.20),
        dpi=600,
        compression="lzw"
    )
}

# ============================================================
# FIGURE 3 — CELL-TYPE EFFECT DISTRIBUTION
#
# This is threshold-independent and therefore complements
# DEG counts.
# ============================================================

plotres <- allres[
    is.finite(allres$logFC),
]

plotres$CellType <- factor(
    plotres$CellType,
    levels=celltypes
)

p3 <- ggplot(
    plotres,
    aes(
        x=CellType,
        y=logFC
    )
) +
    geom_violin(
        trim=TRUE,
        scale="width"
    ) +
    geom_boxplot(
        width=.12,
        outlier.shape=NA
    ) +
    geom_hline(
        yintercept=0,
        linetype=2,
        linewidth=.4
    ) +
    coord_cartesian(
        ylim=quantile(
            plotres$logFC,
            c(.01,.99),
            na.rm=TRUE
        )
    ) +
    labs(
        title="Distribution of AD-associated transcriptional effects by cell type",
        x=NULL,
        y="AD vs Control log2 fold change"
    ) +
    theme_classic(base_size=12) +
    theme(
        axis.text.x=element_text(
            angle=35,
            hjust=1
        ),
        plot.title=element_text(
            face="bold"
        )
    )

ggsave(
    "figures/publication/Figure_SN03_effect_distribution.pdf",
    p3,
    width=8.5,
    height=6
)

ggsave(
    "figures/publication/Figure_SN03_effect_distribution.tiff",
    p3,
    width=8.5,
    height=6,
    dpi=600,
    compression="lzw"
)

# ============================================================
# FIGURE 4 — TOP EFFECTS DOT PLOT
# ============================================================

candidate <- allres[
    is.finite(allres$logFC) &
    is.finite(allres$adj.P.Val) &
    allres$adj.P.Val < .05,
]

if(nrow(candidate)>0){

    # Global top 30 genes by maximum absolute effect
    mx <- aggregate(
        abs(logFC) ~ Gene,
        candidate,
        max
    )

    mx <- mx[
        order(-mx[[2]]),
    ]

    topgenes <- head(
        mx$Gene,
        30
    )

    dot <- allres[
        allres$Gene %in% topgenes,
        c(
            "Gene",
            "CellType",
            "logFC",
            "adj.P.Val"
        )
    ]

    dot$Strength <- -log10(
        pmax(
            dot$adj.P.Val,
            1e-300
        )
    )

    # Cap visual size so extreme P-values do not dominate
    dot$Strength_plot <- pmin(
        dot$Strength,
        20
    )

    dot$CellType <- factor(
        dot$CellType,
        levels=celltypes
    )

    dot$Gene <- factor(
        dot$Gene,
        levels=rev(topgenes)
    )

    p4 <- ggplot(
        dot,
        aes(
            x=CellType,
            y=Gene
        )
    ) +
        geom_point(
            aes(
                size=Strength_plot,
                fill=logFC
            ),
            shape=21
        ) +
        scale_fill_gradient2(
            midpoint=0,
            name="log2FC"
        ) +
        scale_size_continuous(
            name="-log10 FDR\n(capped at 20)",
            range=c(1,7)
        ) +
        labs(
            title="Cell-type specificity of prominent AD-associated genes",
            x=NULL,
            y=NULL
        ) +
        theme_classic(base_size=10) +
        theme(
            axis.text.x=element_text(
                angle=35,
                hjust=1,
                face="bold"
            ),
            plot.title=element_text(
                face="bold"
            )
        )

    ggsave(
        "figures/publication/Figure_SN04_gene_dotplot.pdf",
        p4,
        width=9,
        height=9
    )

    ggsave(
        "figures/publication/Figure_SN04_gene_dotplot.tiff",
        p4,
        width=9,
        height=9,
        dpi=600,
        compression="lzw"
    )
}

cat("\n============================================\n")
cat("SINGLE-NUCLEUS FIGURES COMPLETED\n")
cat("============================================\n")

cat("\nFigures:\n")
print(
    list.files(
        "figures/publication",
        full.names=FALSE
    )
)
