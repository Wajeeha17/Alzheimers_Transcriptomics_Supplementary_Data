suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(ggplot2)
})

regions <- c(
    "EC","HIP","MTG",
    "PC","SFG","VCX"
)

dir.create("results/ORA",recursive=TRUE,showWarnings=FALSE)
dir.create("figures/ORA",recursive=TRUE,showWarnings=FALSE)

convert_genes <- function(symbols){

    bitr(
        unique(symbols),
        fromType="SYMBOL",
        toType="ENTREZID",
        OrgDb=org.Hs.eg.db
    )
}

for(reg in regions){

    cat(
        "\n============================\n",
        reg,
        "\n============================\n",
        sep=""
    )

    f <- paste0(
        "results/gene_level_DEA/",
        reg,
        "_gene_level_AD_vs_Control.tsv"
    )

    if(!file.exists(f))
        next

    x <- read.delim(
        f,
        stringsAsFactors=FALSE
    )

    universe_map <- convert_genes(
        x$Gene
    )

    universe <- unique(
        universe_map$ENTREZID
    )

    for(direction in c("Up","Down")){

        symbols <- x$Gene[
            x$Direction==direction
        ]

        cat(
            direction,
            ":",
            length(symbols),
            "genes\n"
        )

        if(length(symbols)<10){
            cat("Too few genes; skipping.\n")
            next
        }

        conv <- convert_genes(
            symbols
        )

        entrez <- unique(
            conv$ENTREZID
        )

        if(length(entrez)<10)
            next

        # ==========================================
        # GO Biological Process
        # ==========================================

        ego <- enrichGO(
            gene=entrez,
            universe=universe,
            OrgDb=org.Hs.eg.db,
            keyType="ENTREZID",
            ont="BP",
            pAdjustMethod="BH",
            pvalueCutoff=.05,
            qvalueCutoff=.20,
            readable=TRUE
        )

        go_table <- as.data.frame(
            ego
        )

        write.table(
            go_table,
            paste0(
                "results/ORA/",
                reg,"_",direction,
                "_GO_BP.tsv"
            ),
            sep="\t",
            quote=FALSE,
            row.names=FALSE
        )

        # ==========================================
        # KEGG
        # ==========================================

        ekegg <- enrichKEGG(
            gene=entrez,
            universe=universe,
            organism="hsa",
            pAdjustMethod="BH",
            pvalueCutoff=.05,
            qvalueCutoff=.20
        )

        kegg_table <- as.data.frame(
            ekegg
        )

        write.table(
            kegg_table,
            paste0(
                "results/ORA/",
                reg,"_",direction,
                "_KEGG.tsv"
            ),
            sep="\t",
            quote=FALSE,
            row.names=FALSE
        )

        # ==========================================
        # FIGURES
        # ==========================================

        if(nrow(go_table)>0){

            p <- dotplot(
                ego,
                showCategory=15
            ) +
                ggtitle(
                    paste(
                        reg,
                        direction,
                        "GO Biological Process"
                    )
                )

            ggsave(
                paste0(
                    "figures/ORA/",
                    reg,"_",direction,
                    "_GO_BP_dotplot.pdf"
                ),
                p,
                width=9,
                height=7
            )

            ggsave(
                paste0(
                    "figures/ORA/",
                    reg,"_",direction,
                    "_GO_BP_dotplot.tiff"
                ),
                p,
                width=9,
                height=7,
                dpi=600,
                compression="lzw"
            )
        }

        if(nrow(kegg_table)>0){

            p <- dotplot(
                ekegg,
                showCategory=15
            ) +
                ggtitle(
                    paste(
                        reg,
                        direction,
                        "KEGG pathways"
                    )
                )

            ggsave(
                paste0(
                    "figures/ORA/",
                    reg,"_",direction,
                    "_KEGG_dotplot.pdf"
                ),
                p,
                width=9,
                height=7
            )
        }


    }
}

cat("\nREGIONAL ORA COMPLETED\n")
