suppressMessages(library("rARPACK"))
suppressMessages(library("logisticPCA"))
suppressMessages(library("optparse"))

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input (.tsv) of binary matrix"),
  make_option(c("-o", "--outdir"), type="character", help="Output directory")
)

set.seed(42)

parser <- OptionParser(usage="%prog [options] file", option_list=option_list)
args <- parse_args(parser)
opt <- args$options
file <- args$args

X.df <- as.data.frame(read.table(args$input, sep='\t', header=T, row.names=1))
X <- as.matrix(sapply(X.df, as.numeric))

# Determine m via logisticSVD
logsvd_model = logisticSVD(X, k = 2)
logpca_cv = cv.lpca(X, ks = 2, ms = 1:10)
write.table(logpca_cv, file.path(args$outdir, "logistic_pca_cv_nll.tsv") ,sep='\t')

# Run logistic PCA models (we  use logPCA, not convexLogisticPCA)
logpca_model = logisticPCA(X, k = 2, m = which.min(logpca_cv))
clogpca_model = convexLogisticPCA(X, k = 2, m = which.min(logpca_cv))

# Write out Results for PCA Loadings & PCA Components
rownames(logpca_model$U) <- colnames(X)
colnames(logpca_model$U) <- c('PC1','PC2')

write.table(logpca_model$U, file.path(args$outdir, "logistic_pca_loadings.tsv"), sep='\t')
write.table(logpca_model$PCs, file.path(args$outdir, "logistic_pca.tsv"), sep='\t')
