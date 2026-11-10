suppressPackageStartupMessages(library(GEOquery))

ids <- readLines("raw_metadata/candidate_GSE.txt")
ids <- trimws(ids)
ids <- unique(ids[grepl("^GSE[0-9]+$", ids)])

dir.create("geo_cache", showWarnings=FALSE)
dir.create("raw_metadata/GEO", recursive=TRUE, showWarnings=FALSE)
dir.create("expression_data", showWarnings=FALSE)

summary_list <- list()

for (id in ids) {

    cat("\n====================================\n")
    cat("Retrieving:", id, "\n")
    cat("====================================\n")

    tryCatch({

        gset <- getGEO(
            id,
            GSEMatrix=TRUE,
            getGPL=FALSE,
            destdir="geo_cache"
        )

        for (i in seq_along(gset)) {

            eset <- gset[[i]]
            pd   <- pData(eset)
            expr <- exprs(eset)

            platform <- annotation(eset)

            meta_file <- sprintf(
                "raw_metadata/GEO/%s_platform%d_metadata.tsv",
                id, i
            )

            expr_file <- sprintf(
                "expression_data/%s_platform%d_expression.tsv.gz",
                id, i
            )

            write.table(
                pd,
                meta_file,
                sep="\t",
                quote=FALSE,
                row.names=FALSE
            )

            con <- gzfile(expr_file, "wt")

            write.table(
                expr,
                con,
                sep="\t",
                quote=FALSE,
                col.names=NA
            )

            close(con)

            summary_list[[length(summary_list)+1]] <- data.frame(
                dataset_id=id,
                platform_index=i,
                platform=platform,
                samples=ncol(expr),
                features=nrow(expr),
                metadata_columns=ncol(pd),
                stringsAsFactors=FALSE
            )

            cat("Platform:", platform, "\n")
            cat("Samples:", ncol(expr), "\n")
            cat("Features:", nrow(expr), "\n")
            cat("Saved metadata:", meta_file, "\n")
            cat("Saved expression:", expr_file, "\n")
        }

    }, error=function(e) {

        cat("FAILED:", id, "\n")
        cat("Reason:", conditionMessage(e), "\n")

    })
}

if (length(summary_list) > 0) {

    summary_df <- do.call(rbind, summary_list)

    write.table(
        summary_df,
        "tables/GEO_retrieval_summary.tsv",
        sep="\t",
        quote=FALSE,
        row.names=FALSE
    )

    cat("\n=== RETRIEVAL SUMMARY ===\n")
    print(summary_df)

} else {

    cat("\nNo datasets successfully retrieved.\n")

}
