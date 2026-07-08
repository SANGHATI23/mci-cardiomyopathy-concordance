
if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}

needed <- c("DESeq2", "apeglm", "readr", "dplyr", "tibble")

for (pkg in needed) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        if (pkg %in% c("DESeq2", "apeglm")) {
            BiocManager::install(pkg, ask=FALSE, update=FALSE)
        } else {
            install.packages(pkg, repos="https://cloud.r-project.org")
        }
    }
}
