suppressPackageStartupMessages({
    library(dorothea)
    library(ggplot2)
})

regions <- c("EC","HIP","MTG","PC","SFG","VCX")

dir.create("results/TF_activity",recursive=TRUE,showWarnings=FALSE)

data(dorothea_hs,package="dorothea")

reg <- as.data.frame(dorothea_hs)

# Handle different dorothea package versions
if("source" %in% names(reg)){
    reg$TF <- reg$source
} else if("tf" %in% names(reg)){
    reg$TF <- reg$tf
} else {
    stop("TF column not found.")
}

reg <- reg[
    reg$confidence %in% c("A","B","C"),
]

all_activity <- list()

for(region in regions){

    cat("\nProcessing:",region,"\n")

    f <- paste0(
        "input/GSE5281_gene_level_DEA/",
        region,
        "_gene_level_AD_vs_Control.tsv"
    )

    de <- read.delim(
        f,
        stringsAsFactors=FALSE
    )

    # Full ranked transcriptomic signature
    de <- de[
        !is.na(de$Gene) &
        is.finite(de$t),
        c("Gene","t")
    ]

    de <- de[
        !duplicated(de$Gene),
    ]

    TFs <- unique(reg$TF)

    activity <- lapply(TFs,function(tf){

        r <- reg[
            reg$TF==tf,
        ]

        idx <- match(
            r$target,
            de$Gene
        )

        keep <- !is.na(idx)

        r <- r[keep,]
        idx <- idx[keep]

        # Require adequate regulon coverage
        if(length(idx)<10)
            return(NULL)

        target_t <- de$t[idx]

        # Signed regulon score:
        # activated targets contribute positively;
        # repressed targets contribute inversely.
        signed <- target_t * r$mor

        # Weighted if likelihood is usable,
        # otherwise use signed mean.
        if("likelihood" %in% names(r) &&
           is.numeric(r$likelihood)){

            w <- abs(r$likelihood)

            if(all(is.finite(w)) && sum(w)>0){

                score <- weighted.mean(
                    signed,
                    w=w,
                    na.rm=TRUE
                )

            } else {

                score <- mean(
                    signed,
                    na.rm=TRUE
                )
            }

        } else {

            score <- mean(
                signed,
                na.rm=TRUE
            )
        }

        # Standardized score for comparability
        zscore <-
            mean(signed,na.rm=TRUE) /
            (
                sd(signed,na.rm=TRUE) /
                sqrt(length(signed))
            )

        data.frame(
            Region=region,
            TF=tf,
            Targets_used=length(signed),
            Activity_score=score,
            Activity_Z=zscore
        )
    })

    activity <- activity[
        !vapply(
            activity,
            is.null,
            logical(1)
        )
    ]

    activity <- do.call(
        rbind,
        activity
    )

    activity$P_value <-
        2*pnorm(
            -abs(activity$Activity_Z)
        )

    activity$FDR <-
        p.adjust(
            activity$P_value,
            method="BH"
        )

    activity$Direction <-
        ifelse(
            activity$Activity_score>0,
            "Higher inferred activity in AD",
            "Lower inferred activity in AD"
        )

    activity <- activity[
        order(activity$FDR),
    ]

    write.table(
        activity,
        paste0(
            "results/TF_activity/",
            region,
            "_TF_activity.tsv"
        ),
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    all_activity[[region]] <- activity
}

# ============================================================
# COMBINE REGIONS
# ============================================================

combined <- do.call(
    rbind,
    all_activity
)

rownames(combined) <- NULL

write.table(
    combined,
    "tables/GSE5281_TF_activity_all_regions.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# CROSS-REGION CONSISTENCY
# ============================================================

TFs <- unique(combined$TF)

cross <- do.call(
    rbind,
    lapply(TFs,function(tf){

        z <- combined[
            combined$TF==tf,
        ]

        data.frame(
            TF=tf,

            Regions_tested=nrow(z),

            Positive_regions=sum(
                z$Activity_score>0
            ),

            Negative_regions=sum(
                z$Activity_score<0
            ),

            FDR05_regions=sum(
                z$FDR<.05
            ),

            Mean_activity_Z=
                mean(
                    z$Activity_Z,
                    na.rm=TRUE
                ),

            Mean_abs_activity_Z=
                mean(
                    abs(z$Activity_Z),
                    na.rm=TRUE
                )
        )
    })
)

cross$Direction_consistency <-
    pmax(
        cross$Positive_regions,
        cross$Negative_regions
    ) /
    cross$Regions_tested

cross <- cross[
    order(
        -cross$FDR05_regions,
        -cross$Direction_consistency,
        -cross$Mean_abs_activity_Z
    ),
]

write.table(
    cross,
    "tables/GSE5281_cross_region_TF_activity.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# HEATMAP — strongest recurrent TFs
# ============================================================

topTF <- head(
    cross$TF,
    30
)

plotdat <- combined[
    combined$TF %in% topTF,
]

plotdat$TF <- factor(
    plotdat$TF,
    levels=rev(topTF)
)

plotdat$Region <- factor(
    plotdat$Region,
    levels=regions
)

p <- ggplot(
    plotdat,
    aes(
        x=Region,
        y=TF,
        fill=Activity_Z
    )
) +
    geom_tile(
        colour="white",
        linewidth=.3
    ) +
    scale_fill_gradient2(
        midpoint=0,
        name="Activity Z"
    ) +
    labs(
        title="Inferred transcription-factor activity across AD brain regions",
        subtitle="DoRothEA A-C regulons; AD versus control transcriptomic signatures",
        x=NULL,
        y=NULL
    ) +
    theme_minimal(base_size=11) +
    theme(
        panel.grid=element_blank()
    )

ggsave(
    "figures/TF_activity_region_heatmap.pdf",
    p,
    width=8,
    height=9
)

ggsave(
    "figures/TF_activity_region_heatmap.tiff",
    p,
    width=8,
    height=9,
    dpi=600,
    compression="lzw"
)

cat("\n=========================================\n")
cat("TF ACTIVITY INFERENCE COMPLETED\n")
cat("=========================================\n\n")

print(
    head(
        cross,
        20
    )
)
