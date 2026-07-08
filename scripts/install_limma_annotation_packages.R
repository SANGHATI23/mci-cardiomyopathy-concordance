
if (!requireNamespace("BiocManager", quietly=TRUE)) {
    install.packages("BiocManager", repos="https://cloud.r-project.org")
}

needed_bioc <- c("limma", "AnnotationDbi", "org.Hs.eg.db")
needed_cran <- c("readr", "dplyr", "tibble")

for (pkg in needed_bioc) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        BiocManager::install(pkg, ask=FALSE, update=FALSE)
    }
}

for (pkg in needed_cran) {
    if (!requireNamespace(pkg, quietly=TRUE)) {
        install.packages(pkg, repos="https://cloud.r-project.org")
    }
}
