
# Optional Manuscript Text: scRNA-seq Contextualization Module

Bulk heart transcriptomic cohorts are affected by cellular composition, tissue-region sampling, disease stage, and genotype heterogeneity. Therefore, unstable bulk MCI values should not automatically be interpreted as absence of biological relevance. To support future resource expansion, we created a single-cell RNA-seq contextualization module that maps MCI genes onto public human cardiac single-cell or single-nucleus RNA-seq datasets.

The module is designed to annotate each MCI gene by cardiac cell-type expression, including cardiomyocyte, fibroblast, endothelial, smooth-muscle, and immune-cell expression where those labels are available. For each gene, the planned output reports mean expression by cell type, dominant cell type, cell-type specificity score, MCI tier, GTEx confidence status, and a row-level interpretation.

This layer is intended to help distinguish truly unstable bulk transcriptomic evidence from signals that may be diluted by cell-type composition. For example, a sarcomeric gene with low bulk MCI but strong cardiomyocyte-specific expression may require cardiomyocyte-resolved analysis, proteomics, or genotype-stratified validation rather than being deprioritized solely on the basis of bulk RNA-seq instability.
