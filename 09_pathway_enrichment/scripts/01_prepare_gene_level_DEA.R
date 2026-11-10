suppressPackageStartupMessages({
    library(hgu133plus2.db)
    library(AnnotationDbi)
})

regions <- c(
    "EC","HIP","MTG",
    "PC","SFG","VCX"
)

dir.create(
    "results/gene_level_DEA",
    recursive=TRUE,
    showWarnings=FALSE
)

# Representative probe mapping previously generated
mapfile <- "../07_cell_type_deconvolution/tables/GSE5281_selected_probe_per_gene.tsv"

if(!file.exists(mapfile))
    stop("Representative probe mapping not found.")

mapping <- read.delim(
    mapfile,
    stringsAsFactors=FALSE
)

cat(
    "Representative genes:",
    nrow(mapping),"\n"
)

for(reg in regions){

    f <- paste0(
        "input/GSE5281_DEA/",
        reg,
        "_AD_vs_Control_all_probes.tsv"
    )

    if(!file.exists(f)){
        warning("Missing: ",f)
        next
    }

    x <- read.delim(
        f,
        stringsAsFactors=FALSE,
        check.names=FALSE
    )

    # Detect probe ID
    if("probe_id" %in% names(x)){
        x$PROBEID <- x$probe_id
    } else if("ProbeID" %in% names(x)){
        x$PROBEID <- x$ProbeID
    } else {
        x$PROBEID <- rownames(x)
    }

    idx <- match(
        mapping$PROBEID,
        x$PROBEID
    )

    keep <- !is.na(idx)

    gene <- x[
        idx[keep],
        ,
        drop=FALSE
    ]

    gene$Gene <- mapping$SYMBOL[keep]

    gene <- gene[
        !duplicated(gene$Gene),
        ,
        drop=FALSE
    ]

    gene$Direction <- "NS"

    gene$Direction[
        gene$adj.P.Val < .05 &
        gene$logFC >= .5
    ] <- "Up"

    gene$Direction[
        gene$adj.P.Val < .05 &
        gene$logFC <= -.5
    ] <- "Down"

    outfile <- paste0(
        "results/gene_level_DEA/",
        reg,
        "_gene_level_AD_vs_Control.tsv"
    )

    write.table(
        gene,
        outfile,
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    cat(
        reg,
        ":",
        nrow(gene),"genes |",
        sum(gene$Direction=="Up"),"up |",
        sum(gene$Direction=="Down"),"down\n"
    )
}
