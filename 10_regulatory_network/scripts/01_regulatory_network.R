suppressPackageStartupMessages({
    library(dorothea)
    library(igraph)
    library(ggplot2)
})

regions <- c("EC","HIP","MTG","PC","SFG","VCX")

dir.create("results/regional_networks",recursive=TRUE,showWarnings=FALSE)
dir.create("tables",showWarnings=FALSE)
dir.create("figures",showWarnings=FALSE)

# ============================================================
# DOROTHEA HUMAN REGULONS
# ============================================================

data(dorothea_hs,package="dorothea")

regulons <- as.data.frame(dorothea_hs)

# Standardize DoRothEA column names across package versions
if ("source" %in% names(regulons)) {
    regulons$TF <- regulons$TF
} else if ("tf" %in% names(regulons)) {
    regulons$TF <- regulons$tf
} else {
    stop(
        "Could not identify TF column in DoRothEA. Columns: ",
        paste(names(regulons), collapse=", ")
    )
}

if (!"target" %in% names(regulons))
    stop("DoRothEA target column not found.")

if (!"mor" %in% names(regulons))
    stop("DoRothEA mor column not found.")

# High/intermediate-confidence interactions only
regulons <- regulons[
    regulons$confidence %in% c("A","B","C"),
]

cat("DoRothEA A-C regulatory interactions:",
    nrow(regulons),"\n")

cat("Unique TFs:",
    length(unique(regulons$TF)),"\n")

# ============================================================
# WGCNA MODULE INFORMATION
# ============================================================

module_file <- "input/WGCNA_gene_module_assignments.tsv"

if(file.exists(module_file)){

    modules <- read.delim(
        module_file,
        stringsAsFactors=FALSE
    )

} else {

    modules <- NULL
    warning("WGCNA module table not found.")
}

# ============================================================
# ANALYSE EACH BRAIN REGION
# ============================================================

regionalTF <- list()
allEdges <- list()

for(reg in regions){

    cat("\n====================================\n")
    cat("Region:",reg,"\n")
    cat("====================================\n")

    f <- paste0(
        "input/GSE5281_gene_level_DEA/",
        reg,
        "_gene_level_AD_vs_Control.tsv"
    )

    if(!file.exists(f)){
        warning("Missing ",f)
        next
    }

    de <- read.delim(
        f,
        stringsAsFactors=FALSE
    )

    de <- de[
        !is.na(de$Gene) &
        is.finite(de$logFC) &
        is.finite(de$adj.P.Val),
    ]

    # --------------------------------------------------------
    # Candidate disease-associated targets
    # --------------------------------------------------------

    sig <- de[
        de$adj.P.Val < .05 &
        abs(de$logFC) >= .5,
    ]

    cat("Significant genes:",nrow(sig),"\n")

    if(nrow(sig)<5){
        warning("Too few genes in ",reg)
        next
    }

    # --------------------------------------------------------
    # Regulatory interactions targeting DE genes
    # --------------------------------------------------------

    edges <- regulons[
        regulons$target %in% sig$Gene,
    ]

    if(nrow(edges)==0){
        warning("No regulatory edges for ",reg)
        next
    }

    # Add target expression evidence
    idx <- match(
        edges$target,
        sig$Gene
    )

    edges$Target_logFC <- sig$logFC[idx]
    edges$Target_FDR <- sig$adj.P.Val[idx]

    edges$Region <- reg

    # --------------------------------------------------------
    # Is TF itself differentially expressed?
    # --------------------------------------------------------

    tfidx <- match(
        edges$TF,
        de$Gene
    )

    edges$TF_logFC <- de$logFC[tfidx]
    edges$TF_FDR <- de$adj.P.Val[tfidx]

    # --------------------------------------------------------
    # Regulatory consistency
    #
    # mor = mode of regulation:
    # +1 activation
    # -1 repression
    # --------------------------------------------------------

    edges$Expected_target_direction <-
        edges$TF_logFC * edges$mor

    edges$Direction_consistent <-
        ifelse(
            is.finite(edges$Expected_target_direction),
            sign(edges$Expected_target_direction) ==
            sign(edges$Target_logFC),
            NA
        )

    # --------------------------------------------------------
    # WGCNA support
    # --------------------------------------------------------

    if(!is.null(modules)){

        edges$Target_Module <-
            modules$Module[
                match(
                    edges$target,
                    modules$Gene
                )
            ]

        edges$TF_Module <-
            modules$Module[
                match(
                    edges$TF,
                    modules$Gene
                )
            ]

        edges$Same_Module <-
            !is.na(edges$Target_Module) &
            !is.na(edges$TF_Module) &
            edges$Target_Module ==
            edges$TF_Module

    } else {

        edges$Target_Module <- NA
        edges$TF_Module <- NA
        edges$Same_Module <- NA
    }

    write.table(
        edges,
        paste0(
            "results/regional_networks/",
            reg,
            "_TF_target_network.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    allEdges[[reg]] <- edges

    # --------------------------------------------------------
    # TF SUMMARY
    # --------------------------------------------------------

    tfNames <- unique(edges$TF)

    tfsummary <- do.call(
        rbind,
        lapply(tfNames,function(tf){

            z <- edges[
                edges$TF==tf,
            ]

            data.frame(
                Region=reg,
                TF=tf,
                Targets=nrow(z),

                Up_targets=sum(
                    z$Target_logFC > 0,
                    na.rm=TRUE
                ),

                Down_targets=sum(
                    z$Target_logFC < 0,
                    na.rm=TRUE
                ),

                Mean_abs_target_logFC=
                    mean(
                        abs(z$Target_logFC),
                        na.rm=TRUE
                    ),

                Direction_consistent=
                    sum(
                        z$Direction_consistent %in% TRUE,
                        na.rm=TRUE
                    ),

                Same_WGCNA_module=
                    sum(
                        z$Same_Module %in% TRUE,
                        na.rm=TRUE
                    )
            )
        })
    )

    tfsummary$Consistency_fraction <-
        tfsummary$Direction_consistent /
        tfsummary$Targets

    # Descriptive prioritization score
    tfsummary$Regulatory_score <-
        log1p(tfsummary$Targets) *
        tfsummary$Mean_abs_target_logFC

    tfsummary <- tfsummary[
        order(
            -tfsummary$Regulatory_score
        ),
    ]

    write.table(
        tfsummary,
        paste0(
            "tables/",
            reg,
            "_TF_summary.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    regionalTF[[reg]] <- tfsummary
}

# ============================================================
# COMBINE REGIONS
# ============================================================

if(length(allEdges)==0)
    stop("No regulatory networks were generated.")

allEdges <- do.call(
    rbind,
    allEdges
)

rownames(allEdges) <- NULL

write.table(
    allEdges,
    "results/GSE5281_all_regional_TF_target_edges.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

allTF <- do.call(
    rbind,
    regionalTF
)

rownames(allTF) <- NULL

write.table(
    allTF,
    "tables/GSE5281_all_regional_TF_summary.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# CROSS-REGION TF RECURRENCE
# ============================================================

tfNames <- unique(allTF$TF)

recurrence <- do.call(
    rbind,
    lapply(tfNames,function(tf){

        z <- allTF[
            allTF$TF==tf,
        ]

        data.frame(
            TF=tf,

            Regions=length(
                unique(z$Region)
            ),

            Total_targets=sum(
                z$Targets
            ),

            Mean_regulatory_score=
                mean(
                    z$Regulatory_score
                ),

            Mean_target_effect=
                mean(
                    z$Mean_abs_target_logFC
                ),

            Regions_list=
                paste(
                    sort(unique(z$Region)),
                    collapse=","
                )
        )
    })
)

recurrence <- recurrence[
    order(
        -recurrence$Regions,
        -recurrence$Mean_regulatory_score
    ),
]

write.table(
    recurrence,
    "tables/GSE5281_cross_region_TF_recurrence.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# FIGURE 1 — TOP RECURRENT TFs
# ============================================================

top <- head(
    recurrence,
    25
)

top$TF <- factor(
    top$TF,
    levels=rev(top$TF)
)

p1 <- ggplot(
    top,
    aes(
        x=TF,
        y=Regions
    )
) +
    geom_col() +
    coord_flip() +
    labs(
        title="Recurrent transcriptional regulators across AD brain regions",
        subtitle="DoRothEA A-C confidence regulons",
        x=NULL,
        y="Number of brain regions"
    ) +
    theme_classic(base_size=12)

ggsave(
    "figures/RegulatoryNetwork_TF_recurrence.pdf",
    p1,
    width=8,
    height=7
)

ggsave(
    "figures/RegulatoryNetwork_TF_recurrence.tiff",
    p1,
    width=8,
    height=7,
    dpi=600,
    compression="lzw"
)

# ============================================================
# FIGURE 2 — TF x REGION HEATMAP
# ============================================================

selectedTF <- head(
    recurrence$TF,
    30
)

heat <- allTF[
    allTF$TF %in% selectedTF,
    c(
        "TF",
        "Region",
        "Regulatory_score"
    )
]

heat$TF <- factor(
    heat$TF,
    levels=rev(selectedTF)
)

heat$Region <- factor(
    heat$Region,
    levels=regions
)

p2 <- ggplot(
    heat,
    aes(
        x=Region,
        y=TF,
        fill=Regulatory_score
    )
) +
    geom_tile(
        colour="white",
        linewidth=.3
    ) +
    labs(
        title="Regional distribution of candidate transcriptional regulators",
        subtitle="Descriptive TF regulatory score",
        x=NULL,
        y=NULL,
        fill="Regulatory\nscore"
    ) +
    theme_minimal(base_size=11) +
    theme(
        panel.grid=element_blank()
    )

ggsave(
    "figures/RegulatoryNetwork_TF_region_heatmap.pdf",
    p2,
    width=8,
    height=9
)

ggsave(
    "figures/RegulatoryNetwork_TF_region_heatmap.tiff",
    p2,
    width=8,
    height=9,
    dpi=600,
    compression="lzw"
)

cat("\n============================================\n")
cat("REGULATORY NETWORK RECONSTRUCTION COMPLETED\n")
cat("============================================\n")

cat("\nRegulatory edges:",
    nrow(allEdges),"\n")

cat("Unique TFs:",
    length(unique(allEdges$TF)),"\n")

cat("Unique targets:",
    length(unique(allEdges$target)),"\n")

cat("\nTop recurrent TFs:\n")

print(
    head(
        recurrence,
        15
    )
)
