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
    long_form <- merge(long_form, roadmap_sample_info, by.x='Var2', by.y='EID')
    long_form <- ddply(long_form, .(group_name, Var1),
                       function(x) {
                           rep <- x[which.max(x$value),];
                           data.frame(tissue=tail(strsplit(as.character(rep$group_name), ' & ')[[1]], 1),
                                      cluster=rep$Var1, value=rep$value)
                       })
    long_form$cluster <- factor(long_form$cluster, levels=row.names(cluster_density))
    (heatmap(ggplot(long_form, aes(x=cluster, y=tissue, fill=value))) +
     scale_x_discrete(name='Enhancer module') +
     scale_y_discrete(name='Tissue', limits=rev(tissue_ordering), drop=TRUE) +
     scale_fill_gradient(name='Cluster weight', low='white', high='black') +
     theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
           legend.position='right'))
}
