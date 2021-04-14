suppressMessages(library("edgeR"))
suppressMessages(library("limma"))

RunDiffExprAnalysisLimma <- function(counts.df, var.df, covariates.df=NULL, genes.df=NULL, use.sva=TRUE, n.sv=NULL, return_mod=FALSE) {
    # Computes differential expression using a combination of SmartSVA and voom-limma
    #
    # Args:
    #   counts.df: data.frame with read counts (#genes x #samples)
    #   var: variable of interest
    #   genes.df: data.frame mapping gene IDs to gene names
    #   n.sv: number of surrogates for SVA. Automatically determined if NULL.
    #
    # Returns:
    #   voom-limma output augmented with Storey q-values
    #   SVA surrogate variables
    stopifnot(dim(unique(var.df))[1]>1 && dim(var.df)[2]==1)
    var <- colnames(var.df)[1]

    design <- cbind(1, var.df)

    # apply edgeR normalization (TMM) to counts
    dge <- DGEList(counts=counts.df)
    dge <- calcNormFactors(dge)

    # define model
    mod <- model.matrix(as.formula(paste0('~', var)), var.df)

    if (use.sva) {
        cat("  * running SmartSVA\n")
        v <- voom(dge, design=mod)  # run on transformed counts

        Y.r <- t(resid(lm(as.formula(paste0('t(v$E) ~ ', var)), data=var.df)))

        if (is.null(n.sv)){
          n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
        }

        cat(paste0("  * SVs: ", n.sv, "\n"))
        sv.obj <- smartsva.cpp(as.matrix(v$E), mod, mod0=NULL, n.sv=n.sv, alpha=1, B=200, VERBOSE=TRUE)

        # update model to include SVs
        mod <- model.matrix(as.formula(paste0('~', var, '+sv.obj$sv')), var.df)  # intercept, var, SVs
    }

    if (!is.null(covariates.df)) {
        covar.var.df <- cbind(var.df, covariates.df)

        if (use.sva){
          mod <- model.matrix(as.formula(paste0('~', paste(paste(colnames(covar.var.df),collapse = '+'),"+sv.obj$sv"))), covar.var.df)
        } else{
          mod <- model.matrix(as.formula(paste0('~', paste(colnames(covar.var.df),collapse = '+'))), covar.var.df)
        }
    }

    cat(paste0("  * model matrix dimensions: ", dim(mod)[2], "\n"))

    if (return_mod) {
      return(mod)
    }

    # run limma
    v <- voom(dge, design=mod)
    fit <- lmFit(v, design=mod)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=ncol(design), n=Inf, sort.by="none")

    cat(paste0("Differentially expressed genes at 0.05 FDR: ", sum(res[, 'adj.P.Val']<=0.05), "\n"))
    if (!is.null(genes.df)) {
        res[, 'gene_name'] <- genes.df[row.names(res), 'Description']
        res <- res[, colnames(res)[c(8,1:7)]]
    }

    if (use.sva) {
        sva.df <- data.frame(sv.obj$sv)
        rownames(sva.df) <- colnames(counts.df)
        return(list("res"=res, "C"=sva.df))
    } else {
        return(res)
    }
}
