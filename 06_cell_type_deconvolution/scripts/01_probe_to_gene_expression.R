suppressPackageStartupMessages({
    library(hgu133plus2.db)
    library(AnnotationDbi)
})

f <- "../02_data_acquisition_metadata/expression_data/GSE5281_platform1_expression.tsv.gz"

x <- read.delim(
    f,
    row.names=1,
    check.names=FALSE
)

x <- as.matrix(x)
storage.mode(x) <- "numeric"

# GSE5281 is linear-scale processed expression
x[x <= 0] <- NA
x <- log2(x)

keep <- apply(x,1,function(z) all(is.finite(z)))
x <- x[keep,,drop=FALSE]

cat("Input probes:",nrow(x),"\n")

ann <- AnnotationDbi::select(
    hgu133plus2.db,
    keys=rownames(x),
    keytype="PROBEID",
    columns=c("SYMBOL")
)

ann <- ann[
    !is.na(ann$SYMBOL) &
    ann$SYMBOL != "",
]

ann <- ann[
    !duplicated(ann$PROBEID),
]

common <- intersect(
    rownames(x),
    ann$PROBEID
)

x <- x[common,,drop=FALSE]

ann <- ann[
    match(common,ann$PROBEID),
]

# Variance calculated without disease labels
probe_variance <- apply(
    x,
    1,
    var,
    na.rm=TRUE
)

mapping <- data.frame(
    PROBEID=rownames(x),
    SYMBOL=ann$SYMBOL,
    variance=probe_variance,
    stringsAsFactors=FALSE
)

# Highest-variance probe retained for each gene
mapping <- mapping[
    order(
        mapping$SYMBOL,
        -mapping$variance
    ),
]

mapping <- mapping[
    !duplicated(mapping$SYMBOL),
]

gene_expr <- x[
    mapping$PROBEID,
    ,
    drop=FALSE
]

rownames(gene_expr) <- mapping$SYMBOL

write.table(
    mapping,
    "tables/GSE5281_selected_probe_per_gene.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

con <- gzfile(
    "results/GSE5281_gene_expression_log2.tsv.gz",
    "wt"
)

write.table(
    gene_expr,
    con,
    sep="\t",
    quote=FALSE,
    col.names=NA
)

close(con)

cat("Mapped genes:",nrow(gene_expr),"\n")
cat("Samples:",ncol(gene_expr),"\n")
