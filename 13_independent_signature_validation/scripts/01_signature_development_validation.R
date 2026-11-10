suppressPackageStartupMessages({
    library(glmnet)
    library(pROC)
    library(ggplot2)
})

set.seed(20260922)

dir.create("results",showWarnings=FALSE)
dir.create("tables",showWarnings=FALSE)
dir.create("figures",showWarnings=FALSE)
dir.create("models",showWarnings=FALSE)

matched_regions <- c("EC","HIP","SFG")

# ============================================================
# HELPER: READ EXPRESSION MATRIX
# ============================================================

read_expr <- function(file){

    x <- read.delim(
        gzfile(file),
        check.names=FALSE,
        stringsAsFactors=FALSE
    )

    ids <- x[[1]]

    x <- x[,-1,drop=FALSE]

    x <- as.matrix(x)

    storage.mode(x) <- "numeric"

    rownames(x) <- ids

    x
}

# ============================================================
# 1. CANDIDATE GENES
# ============================================================

cat("1. Reading multilayer candidates...\n")

priority <- read.delim(
    "input/multilayer_prioritization.tsv",
    stringsAsFactors=FALSE
)

priority <- priority[
    !is.na(priority$Gene) &
    priority$Gene!="",
]

# Require substantial multilayer support.
# Cap candidate pool to prevent p >> n feature-selection instability.
candidate <- priority[
    priority$Evidence_layers_positive >= 4 &
    priority$DEG_regions >= 2,
]

candidate <- candidate[
    order(
        -candidate$Evidence_layers_positive,
        -candidate$Integrated_score
    ),
]

candidate <- head(
    candidate,
    100
)

candidate_genes <- unique(
    candidate$Gene
)

cat("Candidate genes:",
    length(candidate_genes),"\n")

write.table(
    candidate,
    "tables/signature_candidate_pool.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# 2. GSE5281 DISCOVERY EXPRESSION
# ============================================================

cat("\n2. Reading GSE5281...\n")

expr1 <- read_expr(
    "input/GSE5281_gene_expression.tsv.gz"
)

ph1 <- read.delim(
    "input/GSE5281_phenotype.tsv",
    stringsAsFactors=FALSE
)

ph1 <- ph1[
    ph1$region %in% matched_regions &
    ph1$status %in% c("AD","Control"),
]

common_samples <- intersect(
    ph1$sample_id,
    colnames(expr1)
)

ph1 <- ph1[
    match(
        common_samples,
        ph1$sample_id
    ),
]

expr1 <- expr1[
    ,
    common_samples,
    drop=FALSE
]

cat("Discovery samples:",
    ncol(expr1),"\n")

cat("AD:",
    sum(ph1$status=="AD"),
    " Control:",
    sum(ph1$status=="Control"),
    "\n")

# Candidate genes present in discovery matrix
genes1 <- intersect(
    candidate_genes,
    rownames(expr1)
)

cat("Candidate genes available in GSE5281:",
    length(genes1),"\n")

if(length(genes1)<5)
    stop("Too few candidate genes found in GSE5281.")

# samples x genes
X1 <- t(
    expr1[
        genes1,
        ,
        drop=FALSE
    ]
)

y1 <- ifelse(
    ph1$status=="AD",
    1,
    0
)

# ============================================================
# 3. REMOVE ZERO-VARIANCE FEATURES
# ============================================================

v <- apply(
    X1,
    2,
    var,
    na.rm=TRUE
)

keep <- is.finite(v) & v>0

X1 <- X1[,keep,drop=FALSE]

cat("Features after variance filtering:",
    ncol(X1),"\n")

# ============================================================
# 4. REGION COVARIATES
#
# Genes are penalized.
# Region covariates are included but unpenalized.
# ============================================================

region_matrix <- model.matrix(
    ~ region,
    data=ph1
)[,-1,drop=FALSE]

colnames(region_matrix) <- paste0(
    "REGION_",
    colnames(region_matrix)
)

Xmodel <- cbind(
    X1,
    region_matrix
)

penalty <- c(
    rep(1,ncol(X1)),
    rep(0,ncol(region_matrix))
)

# ============================================================
# 5. REPEATED INTERNAL CROSS-VALIDATION
# ============================================================

cat("\n3. Repeated elastic-net cross-validation...\n")

nrep <- 30
alpha_value <- 0.5

selected_list <- vector(
    "list",
    nrep
)

cv_auc <- rep(
    NA_real_,
    nrep
)

lambda_values <- rep(
    NA_real_,
    nrep
)

for(r in seq_len(nrep)){

    set.seed(
        20260922 + r
    )

    # Stratified 5-fold assignment
    foldid <- integer(
        length(y1)
    )

    for(cl in c(0,1)){

        ii <- which(y1==cl)

        foldid[ii] <- sample(
            rep(
                1:5,
                length.out=length(ii)
            )
        )
    }

    fit <- cv.glmnet(
        x=Xmodel,
        y=y1,
        family="binomial",
        alpha=alpha_value,
        foldid=foldid,
        type.measure="auc",
        standardize=TRUE,
        penalty.factor=penalty
    )

    lambda_values[r] <-
        fit$lambda.1se

    cf <- coef(
        fit,
        s="lambda.1se"
    )

    nz <- rownames(cf)[
        as.vector(cf)!=0
    ]

    # Remove intercept and region covariates
    nz <- setdiff(
        nz,
        c(
            "(Intercept)",
            colnames(region_matrix)
        )
    )

    selected_list[[r]] <- nz

    cv_auc[r] <- max(
        fit$cvm,
        na.rm=TRUE
    )

    cat(
        "Repeat",r,
        " selected",length(nz),
        "genes\n"
    )
}

selection_frequency <- sort(
    table(
        unlist(selected_list)
    ) / nrep,
    decreasing=TRUE
)

selection_table <- data.frame(
    Gene=names(selection_frequency),
    Selection_frequency=
        as.numeric(selection_frequency)
)

write.table(
    selection_table,
    "tables/elastic_net_selection_frequency.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# Stable genes: selected in >=50% repeated CV runs
stable_genes <- selection_table$Gene[
    selection_table$Selection_frequency >= .50
]

# If excessively large, retain the 20 most stable.
if(length(stable_genes)>20){
    stable_genes <- head(
        selection_table$Gene,
        20
    )
}

# If too few survive, take five most recurrent.
if(length(stable_genes)<2){

    stable_genes <- head(
        selection_table$Gene,
        min(
            5,
            nrow(selection_table)
        )
    )
}

if(length(stable_genes)<2)
    stop("Stable multigene signature could not be identified.")

cat("\nStable signature genes:",
    length(stable_genes),"\n")

print(stable_genes)

# ============================================================
# 6. FINAL FROZEN DISCOVERY MODEL
# ============================================================

Xsig <- X1[
    ,
    stable_genes,
    drop=FALSE
]

# Standardization parameters derived ONLY from discovery data
mu <- apply(
    Xsig,
    2,
    mean
)

sdv <- apply(
    Xsig,
    2,
    sd
)

Xsig_z <- sweep(
    Xsig,
    2,
    mu,
    "-"
)

Xsig_z <- sweep(
    Xsig_z,
    2,
    sdv,
    "/"
)

final_data <- data.frame(
    status=y1,
    Xsig_z,
    region=ph1$region,
    check.names=FALSE
)

gene_terms <- paste(
    paste0("`",stable_genes,"`"),
    collapse=" + "
)

formula_final <- as.formula(
    paste(
        "status ~",
        gene_terms,
        "+ region"
    )
)

final_model <- glm(
    formula_final,
    data=final_data,
    family=binomial()
)

saveRDS(
    final_model,
    "models/GSE5281_frozen_signature_model.rds"
)

saveRDS(
    list(
        genes=stable_genes,
        mean=mu,
        sd=sdv,
        matched_regions=matched_regions
    ),
    "models/signature_preprocessing.rds"
)

coef_table <- data.frame(
    Term=names(coef(final_model)),
    Coefficient=as.numeric(
        coef(final_model)
    )
)

write.table(
    coef_table,
    "tables/frozen_signature_coefficients.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# 7. DISCOVERY APPARENT PERFORMANCE
#
# Descriptive only — NOT external validation.
# ============================================================

pred1 <- predict(
    final_model,
    type="response"
)

roc1 <- roc(
    y1,
    pred1,
    quiet=TRUE
)

ci1 <- ci.auc(
    roc1
)

cat(
    "\nDiscovery apparent AUC:",
    as.numeric(auc(roc1)),
    "\n"
)

# ============================================================
# 8. PREPARE GSE48350 USING FROZEN GPL570 MAPPING
# ============================================================

cat("\n4. Preparing independent GSE48350 validation...\n")

mapping <- read.delim(
    "input/frozen_probe_gene_mapping.tsv",
    stringsAsFactors=FALSE,
    check.names=FALSE
)

cat(
    "Mapping columns:",
    paste(names(mapping),collapse=", "),
    "\n"
)

# Automatically identify columns
gene_candidates <- c(
    "SYMBOL","Gene","gene","symbol"
)

probe_candidates <- c(
    "PROBEID","ProbeID","probe_id","probe"
)

gene_col <- gene_candidates[
    gene_candidates %in% names(mapping)
][1]

probe_col <- probe_candidates[
    probe_candidates %in% names(mapping)
][1]

if(
    is.na(gene_col) ||
    is.na(probe_col)
){
    stop(
        "Could not identify probe/gene columns in frozen mapping."
    )
}

expr2_probe <- read_expr(
    "input/GSE48350_expression.tsv.gz"
)

ph2 <- read.delim(
    "input/GSE48350_phenotype.tsv",
    stringsAsFactors=FALSE
)

ph2 <- ph2[
    ph2$region %in% matched_regions &
    ph2$status %in% c("AD","Control"),
]

common2 <- intersect(
    ph2$sample_id,
    colnames(expr2_probe)
)

ph2 <- ph2[
    match(
        common2,
        ph2$sample_id
    ),
]

expr2_probe <- expr2_probe[
    ,
    common2,
    drop=FALSE
]

cat(
    "Validation samples:",
    ncol(expr2_probe),
    "\n"
)

cat(
    "AD:",
    sum(ph2$status=="AD"),
    " Control:",
    sum(ph2$status=="Control"),
    "\n"
)

# Frozen GSE5281-selected probes only
m <- mapping[
    mapping[[gene_col]] %in% stable_genes,
]

m <- m[
    !duplicated(
        m[[gene_col]]
    ),
]

idx <- match(
    m[[probe_col]],
    rownames(expr2_probe)
)

m <- m[
    !is.na(idx),
]

idx <- idx[
    !is.na(idx)
]

cat(
    "Signature genes recovered in validation:",
    nrow(m),
    "/",
    length(stable_genes),
    "\n"
)

if(nrow(m)!=length(stable_genes)){

    missingGenes <- setdiff(
        stable_genes,
        m[[gene_col]]
    )

    cat(
        "Missing:",
        paste(
            missingGenes,
            collapse=", "
        ),
        "\n"
    )

    stop(
        "Not all frozen signature genes were recoverable in GSE48350."
    )
}

expr2 <- expr2_probe[
    idx,
    ,
    drop=FALSE
]

rownames(expr2) <- m[[gene_col]]

expr2 <- expr2[
    stable_genes,
    ,
    drop=FALSE
]

X2 <- t(expr2)

# ============================================================
# IMPORTANT SCALE HANDLING
#
# GSE5281 and GSE48350 processed matrices have different
# numerical distributions. We therefore calculate within-
# cohort gene standardization for score transportability.
#
# This preserves the frozen genes and coefficients but means
# validation evaluates a standardized relative-expression
# signature rather than absolute cross-platform expression.
# ============================================================

X2_z <- scale(X2)

X2_z <- as.matrix(
    X2_z
)

validation_data <- data.frame(
    X2_z,
    region=ph2$region,
    check.names=FALSE
)

y2 <- ifelse(
    ph2$status=="AD",
    1,
    0
)

# ============================================================
# 9. EXTERNAL VALIDATION
# ============================================================

pred2 <- predict(
    final_model,
    newdata=validation_data,
    type="response"
)

roc2 <- roc(
    y2,
    pred2,
    quiet=TRUE
)

auc2 <- as.numeric(
    auc(roc2)
)

ci2 <- ci.auc(
    roc2
)

cat(
    "\nExternal validation AUC:",
    auc2,
    "\n"
)

cat(
    "95% CI:",
    ci2[1],
    "-",
    ci2[3],
    "\n"
)

# ============================================================
# REGION-SPECIFIC EXTERNAL VALIDATION
# ============================================================

performance <- data.frame(
    Cohort="GSE48350",
    Region="Combined",
    N=length(y2),
    AD=sum(y2==1),
    Control=sum(y2==0),
    AUC=auc2,
    CI_low=as.numeric(ci2[1]),
    CI_high=as.numeric(ci2[3])
)

for(reg in matched_regions){

    ii <- which(
        ph2$region==reg
    )

    if(
        length(unique(y2[ii]))<2
    )
        next

    rr <- roc(
        y2[ii],
        pred2[ii],
        quiet=TRUE
    )

    cc <- ci.auc(
        rr
    )

    performance <- rbind(
        performance,
        data.frame(
            Cohort="GSE48350",
            Region=reg,
            N=length(ii),
            AD=sum(y2[ii]==1),
            Control=sum(y2[ii]==0),
            AUC=as.numeric(
                auc(rr)
            ),
            CI_low=as.numeric(
                cc[1]
            ),
            CI_high=as.numeric(
                cc[3]
            )
        )
    )
}

write.table(
    performance,
    "tables/GSE48350_external_validation_performance.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# 10. PREDICTIONS
# ============================================================

pred_table <- data.frame(
    sample_id=ph2$sample_id,
    region=ph2$region,
    status=ph2$status,
    Signature_score=
        as.numeric(pred2)
)

write.table(
    pred_table,
    "results/GSE48350_signature_predictions.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

# ============================================================
# 11. ROC FIGURE
# ============================================================

pdf(
    "figures/Signature_external_validation_ROC.pdf",
    width=6.5,
    height=6.5
)

plot(
    roc2,
    legacy.axes=TRUE,
    print.auc=TRUE,
    main="Independent validation of transcriptomic signature"
)

abline(
    a=0,
    b=1,
    lty=2
)

dev.off()

tiff(
    "figures/Signature_external_validation_ROC.tiff",
    width=6.5,
    height=6.5,
    units="in",
    res=600,
    compression="lzw"
)

plot(
    roc2,
    legacy.axes=TRUE,
    print.auc=TRUE,
    main="Independent validation of transcriptomic signature"
)

abline(
    a=0,
    b=1,
    lty=2
)

dev.off()

# ============================================================
# 12. SIGNATURE SCORE DISTRIBUTION
# ============================================================

pd <- pred_table

p <- ggplot(
    pd,
    aes(
        x=status,
        y=Signature_score
    )
) +
    geom_boxplot(
        outlier.shape=NA
    ) +
    geom_jitter(
        width=.15,
        alpha=.55,
        size=1.5
    ) +
    facet_wrap(
        ~region,
        scales="free_y"
    ) +
    labs(
        title="Independent validation of the transcriptomic signature",
        subtitle="Frozen GSE5281-derived gene signature evaluated in GSE48350",
        x=NULL,
        y="Signature score"
    ) +
    theme_classic(base_size=11)

ggsave(
    "figures/Signature_score_by_region.pdf",
    p,
    width=9,
    height=4.5
)

ggsave(
    "figures/Signature_score_by_region.tiff",
    p,
    width=9,
    height=4.5,
    dpi=600,
    compression="lzw"
)

# ============================================================
# 13. SUMMARY
# ============================================================

cv_summary <- data.frame(
    Metric=c(
        "Repeated_CV_runs",
        "Median_internal_CV_AUC",
        "IQR_internal_CV_AUC",
        "Signature_genes"
    ),
    Value=c(
        nrep,
        median(
            cv_auc,
            na.rm=TRUE
        ),
        IQR(
            cv_auc,
            na.rm=TRUE
        ),
        length(
            stable_genes
        )
    )
)

write.table(
    cv_summary,
    "tables/GSE5281_internal_CV_summary.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n============================================\n")
cat("SIGNATURE DEVELOPMENT/VALIDATION COMPLETED\n")
cat("============================================\n")

cat(
    "\nSignature genes:\n",
    paste(
        stable_genes,
        collapse=", "
    ),
    "\n"
)

cat(
    "\nMedian repeated-CV AUC:",
    median(
        cv_auc,
        na.rm=TRUE
    ),
    "\n"
)

cat(
    "External GSE48350 AUC:",
    auc2,
    "\n"
)

cat("\nExternal performance by region:\n")

print(
    performance,
    row.names=FALSE
)
