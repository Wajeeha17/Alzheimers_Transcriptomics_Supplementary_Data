m <- read.delim(
    "metadata/GSE48350_raw_metadata.tsv",
    check.names=FALSE,
    quote="",
    stringsAsFactors=FALSE
)

sample_id <- m$geo_accession
title <- m$title
source <- m$source_name_ch1

# ---------------------------------------------------------
# Diagnosis
# ---------------------------------------------------------
status <- ifelse(
    grepl("_AD$|_AD\\b", source, ignore.case=TRUE) |
    grepl("_AD$|_AD\\b", title, ignore.case=TRUE),
    "AD",
    ifelse(
        grepl("indiv", title, ignore.case=TRUE),
        "Control",
        NA_character_
    )
)

# ---------------------------------------------------------
# Brain region
# ---------------------------------------------------------
txt <- paste(title, source)

region <- ifelse(
    grepl("entorhinal", txt, ignore.case=TRUE), "EC",
    ifelse(
        grepl("hippocamp", txt, ignore.case=TRUE), "HIP",
        ifelse(
            grepl("superior.?frontal", txt, ignore.case=TRUE), "SFG",
            ifelse(
                grepl("post.?central|postcentral", txt, ignore.case=TRUE),
                "PCG",
                NA_character_
            )
        )
    )
)

# ---------------------------------------------------------
# Sex
# ---------------------------------------------------------
sex <- ifelse(
    grepl("female", txt, ignore.case=TRUE), "Female",
    ifelse(
        grepl("male", txt, ignore.case=TRUE), "Male",
        NA_character_
    )
)

# ---------------------------------------------------------
# Age from title
# ---------------------------------------------------------
age <- suppressWarnings(
    as.numeric(
        sub(
            ".*_([0-9]+)yrs_.*",
            "\\1",
            title
        )
    )
)

age[!grepl("[0-9]+yrs",title)] <- NA

# ---------------------------------------------------------
# Individual ID
# ---------------------------------------------------------
individual_id <- ifelse(
    grepl("indiv",title,ignore.case=TRUE),
    sub(
        ".*_(indiv[^_]+)$",
        "\\1",
        title,
        ignore.case=TRUE
    ),
    NA_character_
)

pheno <- data.frame(
    sample_id=sample_id,
    title=title,
    source=source,
    status=status,
    region=region,
    sex=sex,
    age=age,
    individual_id=individual_id,
    stringsAsFactors=FALSE
)

write.table(
    pheno,
    "metadata/GSE48350_curated_phenotype.tsv",
    sep="\t",
    quote=FALSE,
    row.names=FALSE
)

cat("\n===== DIAGNOSIS =====\n")
print(table(pheno$status,useNA="ifany"))

cat("\n===== REGION × DIAGNOSIS =====\n")
print(addmargins(table(
    pheno$region,
    pheno$status,
    useNA="ifany"
)))

cat("\n===== SEX × DIAGNOSIS =====\n")
print(table(pheno$sex,pheno$status,useNA="ifany"))

cat("\nAge summary:\n")
print(tapply(pheno$age,pheno$status,summary))

cat("\nUnclassified diagnosis:",
    sum(is.na(pheno$status)),"\n")

cat("Unclassified region:",
    sum(is.na(pheno$region)),"\n")
