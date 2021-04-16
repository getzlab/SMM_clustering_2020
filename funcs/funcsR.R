suppressMessages(library(fgsea))
suppressMessages(library('ggplot2'))
suppressMessages(library('ggpubr'))
suppressMessages(library('ggdendro'))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(xCell))

fixID <- function(s){
    s <- as.character(s)
    x <- as.integer(strsplit(s, "_")[[1]][2])
    x <- x+1
    return(paste("C",as.character(x),sep=''))
}

fixID2 <- function(s){
    x <- as.integer(s)
    x <- x+1
    return(paste("C",as.character(x),sep=''))
}

#----------------------------------
# fGSEA funcs
#----------------------------------
runGSEA <- function(df, id, how, GMTS, filt.pval=NULL){
    df.filt <- df[(df$id==id),]

    if(!is.null(filt.pval)){
        df.filt <- df.filt[df.filt$adj.P.Val < filt.pval,]
    }

    if (how=='logfc'){
        df.filt$rank <- df.filt$logFC
    } else{
        df.filt$rank <- -log(df.filt$adj.P.Val) * df.filt$logFC
    }

    R <- df.filt$rank
    names(R) <- df.filt$gene_name
    e.df <- fgsea(GMTS, R, 10000, minSize=3)
    e.df$id <- id

    return(e.df[order(e.df$padj),])
}

runAllGSEA <- function(df, GMTS, how="logfc", filt.pval=NULL, seed=NULL){
    if(!is.null(seed)){
      set.seed(seed)
    }

    e.df <- list()
    c=1
    for (i in unique(df$id)){
        e.df[[c]] = runGSEA(de.df, i, how, GMTS, filt.pval=filt.pval)
        c = c + 1
    }

    e.df <- do.call("rbind", e.df)
    e.df$leadingEdge <- sapply(e.df$leadingEdge, paste, collapse=",")
    return(e.df)
}

#----------------------------------
# Plotting Funcs
#----------------------------------
plotTPM <- function(
    df,
    y,
    my_comparisons,
    gene.df,
    x='consensus_nmf',
    w=8,
    h=4,
    label.y=10,
    label.x=1,
    ylim.min=0,
    ylim.max=NA,
    label.pos="right",
    palette="jco",
    title=NULL
    ){
    id <- rownames(gene.df[gene.df$Description==y,])[1]

    if(is.null(title)){
        title <- y
    }

    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggboxplot(
        df,
        x = x,
        y = id,
        color = x,
        palette = palette,
        add='jitter',
        title=title,
        ylab='logTPM+1',
        xlab='',
        font.label = list(size = 6, color = "black")
    ) +
    ylim(ylim.min,ylim.max) +
    stat_compare_means(comparisons = my_comparisons, method="wilcox.test")+
    stat_compare_means(label.y = label.y, label.x=label.x)+
    theme(legend.position=label.pos)

    return(p)
}

plotSig <- function(df, y, my_comparisons, x='consensus_nmf', w=8, h=4, label.y=10, label.x=1, ylim.min=0, ylim.max=NA, label.pos="right"){
    options(repr.plot.width=w, repr.plot.height=h)

    ggboxplot(
        df,
        x = x,
        y = y,
        color = x,
        palette = "jco",
        add='jitter',
        title=y,
        ylab='Mean logTPM+1',
        xlab='',
        font.label = list(size = 6, color = "black"),
        label.pos=label.pos
    ) +
    ylim(ylim.min,ylim.max) +
    stat_compare_means(comparisons = my_comparisons, method="wilcox.test")+
    stat_compare_means(label.y = label.y, label.x=label.x)+
    theme(legend.position="right")

}

plotVolcano <- function(de.df, w=10, h=12, xlim=NA, ylim=NA, ...){
    options(repr.plot.width=w, repr.plot.height=h)

    de.df$gene_name <- as.character(de.df$gene_name)

    keyvals <- ifelse(
    de.df$adj.P.Val > 0.1, 'grey',
      ifelse(abs(de.df$logFC) < 0.5, 'grey',
          ifelse(
              de.df$logFC < 0, 'royalblue',
              'red3'
          )
        )

    )

    names(keyvals)[keyvals == 'grey'] <- 'NS'
    names(keyvals)[keyvals == 'red3'] <- 'Up'
    names(keyvals)[keyvals == 'royalblue'] <- 'Down'

    EnhancedVolcano(de.df,
        lab = de.df$gene_name,
        x = 'logFC',
        y = 'adj.P.Val',
        pCutoff=0.1,
        FCcutoff=0.5,
        xlim=c(-xlim,xlim),
        ylim=c(0,ylim),
        col=c('grey', 'grey', 'grey', 'red3'),
        ylab = bquote(~-Log[10]~ 'Adj. P-val'),
        colCustom = keyvals,
        ...
    )
}

plotGSEA <- function(e_df, pval.thresh=0.1, filter=NULL, palette='RdBu', h=13, w=15, s_color='black', fix_id=T){
    if (fix_id){
        e_df$id <- sapply(e_df$id,fixID)
    }
    e_df$sig <- e_df$padj<pval.thresh
    e_df$logpval <- -log10(e_df$padj)
    e_df <- e_df[e_df$pathway %in% e_df[e_df$sig,]$pathway,]

    if(!is.null(filter)){
        e_df <- dplyr::filter(e_df, grepl(filter, pathway))
    }

    ### Order axis by dendrogram
    # Load data
    X <- e_df[,c('pathway','id','NES')]
    X <- reshape(X[,c('pathway','id','NES')], timevar='id', idvar='pathway', direction='wide',)
    rownames(X) <- X$pathway
    X$pathway <- NULL

    X[is.na(X)] <- 0

    # Build the dendrogram
    dend <- as.dendrogram(hclust(d = dist(x = X)))
    dendro.plot <- ggdendrogram(dend,rotate = TRUE)

    # Use dendrogram order to order colomn
    order <- order.dendrogram(dend) # dendrogram order
    e_df$pathway <- factor(x = e_df$pathway, levels = unique(e_df$pathway)[order], ordered = TRUE)

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        e_df,
        x="id",
        y="pathway",
        fill = "NES",
        size="logpval",
        color=ifelse(e_df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(e_df$NES)) * c(-1, 1))+
        labs(x="", y="", fill="Enrichment", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=0))

    return(p)
}

plotGSEA_v2 <- function(e_df, pval.thresh=0.1, filter=NULL, palette='RdBu', h=13, w=15, s_color='black', ncol=NA){
    e_df$sig <- e_df$padj<pval.thresh
    e_df$logpval <- -log10(e_df$padj)
    e_df <- e_df[e_df$pathway %in% e_df[e_df$sig,]$pathway,]

    if(!is.null(filter)){
        e_df <- dplyr::filter(e_df, grepl(filter, pathway))
    }

    ### Order axis by dendrogram
    # Load data
    X <- e_df[,c('pathway','id','NES')]
    X <- reshape(X[,c('pathway','id','NES')], timevar='id', idvar='pathway', direction='wide',)
    rownames(X) <- X$pathway
    X$pathway <- NULL

    X[is.na(X)] <- 0

    # Build the dendrogram
    dend <- as.dendrogram(hclust(d = dist(x = X)))
    dendro.plot <- ggdendrogram(dend,rotate = TRUE)

    # Use dendrogram order to order colomn
    order <- order.dendrogram(dend) # dendrogram order
    e_df$pathway <- factor(x = e_df$pathway, levels = unique(e_df$pathway)[order], ordered = TRUE)

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        e_df,
        x="id",
        y="pathway",
        fill = "NES",
        size="logpval",
        color=ifelse(e_df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(e_df$NES)) * c(-1, 1))+
        labs(x="", y="", fill="Enrichment", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=0))+
        facet_grid(grouping ~ ., scales = "free", space = "free")

    return(p)
}
