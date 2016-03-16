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
