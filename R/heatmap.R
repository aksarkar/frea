requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('gdata')
requireNamespace('grid')
requireNamespace('gridExtra')
requireNamespace('gtable')
requireNamespace('reshape2')

#' @include scales.R theme_nature.R
NULL

data(roadmap_sample_info)

scale_heatmap <- function(...) {
    scale_fill_gradient(low='#fee8c8', high='#e34a33', ...)
}

heatmap <- function(plot) {
    plot + geom_raster() + coord_fixed() + theme_nature
}

sparse_heatmap <- function(plot) {
    heatmap(plot) + theme(panel.grid.major=element_line(color='gray90'))
}

epigenome_by_tissue <- function(keep=NULL) {
    if (is.null(keep)) {
        keep <- roadmap_sample_info$EID
    }
    eids <- subset(roadmap_sample_info, EID %in% keep)
    eids$EID <- factor(eids$EID, levels=eid_ordering)
    (heatmap(ggplot(eids, aes(y=EID, x=rep(1), fill=EID))) +
     scale_y_discrete(limits=rev(eid_ordering)) +
     fill_by_eid +
     theme(axis.ticks=element_blank(),
           axis.title=element_blank(),
           axis.text=element_blank(),
           axis.line=element_blank()))
}

cluster_by_tissue <- function(cluster_density, keep=NULL) {
    if (is.null(keep)) {
        keep <- row.names(cluster_density)
    }
    (heatmap(ggplot(data.frame(x=keep, y=rep(1)), aes(x, y, fill=factor(x)))) +
     ## scale_x_discrete(limits=row.names(cluster_density), drop=TRUE) +
     scale_y_discrete() +
     scale_fill_manual(values=color_by_cluster_top(cluster_density=cluster_density)) +
     theme(axis.ticks=element_blank(),
           axis.title=element_blank(),
           axis.text=element_blank(),
           axis.line=element_blank()))
}    

density_by_cluster <- function(cluster_density, keep=NULL) {
    long_form <- melt(cluster_density)
    if (!is.null(keep)) {
        long_form <- subset(long_form, Var1 %in% keep)
    }
    long_form$Var1 <- factor(long_form$Var1, levels=row.names(cluster_density))
    long_form$Var2 <- factor(long_form$Var2, levels=rev(eid_ordering))
    (heatmap(ggplot(long_form, aes(x=Var1, y=Var2, fill=value))) +
     scale_x_discrete(name='Enhancer module') +
     scale_y_discrete(name='Reference epigenome') +
     scale_fill_gradient(name='Cluster weight', low='white', high='black') +
     theme(axis.text=element_blank()))
}

master_regulator_by_cofactor <- function(cofactors) {
    cofactors$V3 <- droplevels(cofactors$V3)
    cofactors$V4 <- droplevels(cofactors$V4)
    match_matrix <- acast(cofactors, formula=V3 ~ V4, value.var='V5', fill=0, fun.aggregate=max)
    row_order <- hclust(dist(match_matrix))$order
    cofactors$V3 <- gdata::reorder.factor(cofactors$V3, new.order=row_order)
    col_order <- apply(match_matrix, 2, function (x) {row_order[which.max(x)]})
    cofactors$V4 <- gdata::reorder.factor(cofactors$V4, new.levels=col_order)
    (sparse_heatmap(ggplot(cofactors, aes(x=V4, y=V3, fill=V5))) +
     scale_heatmap() +
     scale_x_discrete(name='Cofactor', limits=levels(cofactors$V4)) +
     scale_y_discrete(name='Master regulator', limits=rev(levels(cofactors$V3))) +
     theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)))
}

master_regulator_by_phenotype <- function(enrichments) {
    enrichment_matrix <- acast(enrichments, formula=V1 ~ V3, value.var='V4', fill=0, fun.aggregate=max)
    row_order <- hclust(dist(enrichment_matrix))$order
    enrichments$V1 <- with(enrichments, factor(V1, levels=row.names(enrichment_matrix)[order(row_order)]))
    col_order <- apply(enrichment_matrix, 2, function (x) {row_order[which.max(x)]})
    enrichments$V3 <- with(enrichments, factor(V3, levels=levels(V3)[order(col_order)]))
    (sparse_heatmap(ggplot(enrichments, aes(x=V3, y=V1, fill=log10(V4)))) +
     scale_heatmap() +
     scale_y_discrete(limits=rev(levels(enrichments$V1))) +
     labs(x='Master regulator', y='Phenotype') +
     theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)))
}

plot_motif_enrichments <- function(filename) {
    data(honeybadger_cluster_density)
    enrichments <- read.delim(gzfile(filename), header=FALSE, sep=' ')
    enrichments$V1 <- toupper(sub('^[^-]*-', '', enrichments$V1))
    enrichments$V2 <- factor(enrichments$V2, levels=row.names(honeybadger_cluster_density))
    Cairo(type='pdf', file=sub('.txt.gz$', '.pdf', filename), width=120, height=50, units='mm')
    print(master_regulator_by_phenotype(enrichments))
    dev.off()
}

plot_motif_cofactors <- function(filename) {
    data(honeybadger_cluster_density)
    cofactors <- read.delim(gzfile(filename), header=FALSE, sep=' ')
    cofactors$V1 <- toupper(sub('^[^-]*-', '', cofactors$V1))
    cofactors$V2 <- factor(cofactors$V2, levels=row.names(honeybadger_cluster_density))
    d_ply(cofactors, .(V1),
          function (x) {
              Cairo(type='pdf', file=sprintf('%s/%s-cofactors.pdf', dirname(filename), x$V1[1]),
                    width=120, height=60, units='mm')
              print(master_regulator_by_cofactor(x))
              dev.off()
          })}

