suppressPackageStartupMessages({
    library(WGCNA)
})

options(stringsAsFactors=FALSE)
allowWGCNAThreads(nThreads=2)

dir.create("results",showWarnings=FALSE)
dir.create("tables",showWarnings=FALSE)
dir.create("figures",showWarnings=FALSE)

# ============================================================
# EXPRESSION
# ============================================================

expr <- read.delim(
    "input/GSE5281_gene_expression.tsv.gz",
    row.names=1,
    check.names=FALSE
)

expr <- as.matrix(expr)
storage.mode(expr) <- "numeric"

cat("Initial matrix:",
    nrow(expr),"genes x",
    ncol(expr),"samples\n")

# WGCNA expects samples x genes
datExpr <- t(expr)

rm(expr)
gc()

# ============================================================
# PHENOTYPE
# ============================================================

pheno <- read.delim(
    "input/GSE5281_phenotype.tsv",
    stringsAsFactors=FALSE
)

idx <- match(
    rownames(datExpr),
    pheno$sample_id
)

if(any(is.na(idx)))
    stop("Expression/phenotype sample mismatch.")

pheno <- pheno[idx,]

stopifnot(
    all(
        rownames(datExpr)==
        pheno$sample_id
    )
)

# ============================================================
# QUALITY CONTROL
# ============================================================

gsg <- goodSamplesGenes(
    datExpr,
    verbose=3
)

if(!gsg$allOK){

    datExpr <- datExpr[
        gsg$goodSamples,
        gsg$goodGenes
    ]

    pheno <- pheno[
        gsg$goodSamples,
    ]
}

cat(
    "After QC:",
    nrow(datExpr),"samples x",
    ncol(datExpr),"genes\n"
)

# ============================================================
# TOP 5000 VARIABLE GENES
# ============================================================

geneVariance <- apply(
    datExpr,
    2,
    var,
    na.rm=TRUE
)

geneVariance <- sort(
    geneVariance,
    decreasing=TRUE
)

nKeep <- min(
    5000,
    length(geneVariance)
)

keepGenes <- names(
    geneVariance[
        seq_len(nKeep)
    ]
)

datExpr <- datExpr[
    ,
    keepGenes,
    drop=FALSE
]

write.table(
    data.frame(
        Gene=keepGenes,
        Variance=geneVariance[keepGenes]
    ),
    "tables/WGCNA_selected_variable_genes.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat(
    "Genes entering WGCNA:",
    ncol(datExpr),"\n"
)

# ============================================================
# SAMPLE CLUSTERING
# ============================================================

pdf(
    "figures/WGCNA_sample_clustering.pdf",
    width=12,
    height=7
)

sampleTree <- hclust(
    dist(datExpr),
    method="average"
)

plot(
    sampleTree,
    main="GSE5281 sample clustering",
    xlab="",
    sub="",
    cex=.6
)

dev.off()

# ============================================================
# SOFT-THRESHOLD SELECTION
# ============================================================

powers <- c(
    1:10,
    seq(12,20,2)
)

sft <- pickSoftThreshold(
    datExpr,
    powerVector=powers,
    networkType="signed",
    corFnc="bicor",
    verbose=5
)

fit <- sft$fitIndices

write.table(
    fit,
    "tables/WGCNA_soft_threshold_statistics.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

pdf(
    "figures/WGCNA_soft_threshold.pdf",
    width=10,
    height=5
)

par(mfrow=c(1,2))

plot(
    fit[,1],
    -sign(fit[,3])*fit[,2],
    xlab="Soft threshold",
    ylab="Signed R²",
    type="n",
    main="Scale-free topology"
)

text(
    fit[,1],
    -sign(fit[,3])*fit[,2],
    labels=fit[,1],
    cex=.8
)

abline(
    h=.80,
    lty=2
)

plot(
    fit[,1],
    fit[,5],
    xlab="Soft threshold",
    ylab="Mean connectivity",
    type="n",
    main="Mean connectivity"
)

text(
    fit[,1],
    fit[,5],
    labels=fit[,1],
    cex=.8
)

dev.off()

# ============================================================
# AUTOMATIC POWER SELECTION
# ============================================================

signedR2 <- -sign(fit[,3])*fit[,2]

candidate <- fit[,1][
    which(signedR2 >= .80)
]

if(length(candidate)>0){

    softPower <- min(candidate)

} else {

    softPower <- fit[,1][
        which.max(signedR2)
    ]
}

cat(
    "Selected soft power:",
    softPower,"\n"
)

writeLines(
    paste(
        "Selected soft power:",
        softPower
    ),
    "results/WGCNA_selected_power.txt"
)

# ============================================================
# NETWORK CONSTRUCTION
# ============================================================

net <- blockwiseModules(
    datExpr,

    power=softPower,

    networkType="signed",

    TOMType="signed",

    corType="bicor",

    minModuleSize=30,

    mergeCutHeight=.25,

    numericLabels=TRUE,

    pamRespectsDendro=FALSE,

    maxBlockSize=5000,

    saveTOMs=FALSE,

    verbose=3
)

moduleColors <- labels2colors(
    net$colors
)

# ============================================================
# MODULE ASSIGNMENTS
# ============================================================

module_table <- data.frame(
    Gene=colnames(datExpr),
    Module=moduleColors
)

write.table(
    module_table,
    "tables/WGCNA_gene_module_assignments.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

module_counts <- as.data.frame(
    table(moduleColors)
)

names(module_counts) <- c(
    "Module",
    "Genes"
)

write.table(
    module_counts,
    "tables/WGCNA_module_sizes.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# DENDROGRAM
# ============================================================

pdf(
    "figures/WGCNA_gene_dendrogram_modules.pdf",
    width=14,
    height=7
)

plotDendroAndColors(
    net$dendrograms[[1]],
    moduleColors[
        net$blockGenes[[1]]
    ],
    "Module",
    dendroLabels=FALSE,
    hang=.03,
    addGuide=TRUE,
    guideHang=.05
)

dev.off()

# ============================================================
# MODULE EIGENGENES
# ============================================================

MEs <- orderMEs(
    net$MEs
)

write.table(
    data.frame(
        sample_id=rownames(datExpr),
        MEs,
        check.names=FALSE
    ),
    "results/WGCNA_module_eigengenes.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# TRAITS
# AD + REGIONS
# ============================================================

traits <- data.frame(
    AD=as.numeric(
        pheno$status=="AD"
    )
)

regions <- c(
    "EC","HIP","MTG",
    "PC","SFG","VCX"
)

for(r in regions){

    traits[[r]] <- as.numeric(
        pheno$region==r
    )
}

rownames(traits) <- pheno$sample_id

# ============================================================
# MODULE-TRAIT ASSOCIATIONS
# ============================================================

moduleTraitCor <- cor(
    MEs,
    traits,
    use="p"
)

moduleTraitP <- corPvalueStudent(
    moduleTraitCor,
    nrow(datExpr)
)

write.table(
    moduleTraitCor,
    "tables/WGCNA_module_trait_correlations.tsv",
    sep="\t",
    quote=FALSE,
    col.names=NA
)

write.table(
    moduleTraitP,
    "tables/WGCNA_module_trait_Pvalues.tsv",
    sep="\t",
    quote=FALSE,
    col.names=NA
)

# BH correction across all module-trait tests
FDR <- matrix(
    p.adjust(
        as.vector(moduleTraitP),
        method="BH"
    ),
    nrow=nrow(moduleTraitP),
    ncol=ncol(moduleTraitP),
    dimnames=dimnames(moduleTraitP)
)

write.table(
    FDR,
    "tables/WGCNA_module_trait_FDR.tsv",
    sep="\t",
    quote=FALSE,
    col.names=NA
)

# ============================================================
# MODULE-TRAIT HEATMAP
# ============================================================

textMatrix <- paste0(
    sprintf("%.2f",moduleTraitCor),
    "\nFDR=",
    sprintf("%.2g",FDR)
)

dim(textMatrix) <- dim(
    moduleTraitCor
)

pdf(
    "figures/WGCNA_module_trait_heatmap.pdf",
    width=10,
    height=max(
        6,
        nrow(moduleTraitCor)*.45
    )
)

labeledHeatmap(
    Matrix=moduleTraitCor,

    xLabels=names(traits),

    yLabels=names(MEs),

    ySymbols=names(MEs),

    colorLabels=FALSE,

    colors=blueWhiteRed(50),

    textMatrix=textMatrix,

    setStdMargins=FALSE,

    cex.text=.65,

    zlim=c(-1,1),

    main="Module–trait relationships"
)

dev.off()

# ============================================================
# MODULE MEMBERSHIP
# ============================================================

MM <- cor(
    datExpr,
    MEs,
    use="p"
)

MMp <- corPvalueStudent(
    MM,
    nrow(datExpr)
)

write.table(
    MM,
    "results/WGCNA_module_membership.tsv",
    sep="\t",
    quote=FALSE,
    col.names=NA
)

# ============================================================
# GENE SIGNIFICANCE FOR AD
# ============================================================

GS_AD <- cor(
    datExpr,
    traits$AD,
    use="p"
)

GS_AD_p <- corPvalueStudent(
    GS_AD,
    nrow(datExpr)
)

gs <- data.frame(
    Gene=rownames(GS_AD),
    GS_AD=as.numeric(GS_AD),
    P_value=as.numeric(GS_AD_p)
)

gs$FDR <- p.adjust(
    gs$P_value,
    method="BH"
)

write.table(
    gs,
    "results/WGCNA_gene_significance_AD.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# HUB CANDIDATES
#
# Rank genes within each module by module membership.
# This is network centrality, not automatically biological
# importance.
# ============================================================

hubList <- list()

for(mod in setdiff(
    unique(moduleColors),
    "grey"
)){

    genes <- module_table$Gene[
        module_table$Module==mod
    ]

    meName <- paste0(
        "ME",
        mod
    )

    if(!meName %in% colnames(MM))
        next

    tab <- data.frame(
        Gene=genes,
        Module=mod,
        ModuleMembership=MM[
            genes,
            meName
        ],
        GS_AD=gs$GS_AD[
            match(
                genes,
                gs$Gene
            )
        ]
    )

    tab <- tab[
        order(
            -abs(tab$ModuleMembership)
        ),
    ]

    hubList[[mod]] <- head(
        tab,
        50
    )
}

hubs <- do.call(
    rbind,
    hubList
)

write.table(
    hubs,
    "tables/WGCNA_top_module_hub_candidates.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# SAVE OBJECT
# ============================================================

save(
    net,
    MEs,
    moduleColors,
    traits,
    datExpr,
    pheno,
    softPower,
    file="results/GSE5281_WGCNA_objects.RData"
)

cat("\n========================================\n")
cat("WGCNA COMPLETED\n")
cat("========================================\n")

cat(
    "\nSamples:",
    nrow(datExpr)
)

cat(
    "\nGenes:",
    ncol(datExpr)
)

cat(
    "\nSoft power:",
    softPower
)

cat(
    "\nModules:",
    length(
        setdiff(
            unique(moduleColors),
            "grey"
        )
    ),
    "\n\n"
)

print(
    module_counts[
        order(
            -module_counts$Genes
        ),
    ]
)
