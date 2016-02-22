requireNamespace('Cairo')
requireNamespace('ggplot2')

#' @include scales.R theme_nature.R
NULL

enrichment_by_cluster <- function(enrichments, cluster_density) {
    (ggplot(enrichments, aes(x=V3, y=(V5 - V6)/ sqrt(V7), fill=V3)) +
     geom_bar(stat='identity') +
     scale_x_discrete(name='Enhancer module', drop=FALSE) +
     scale_y_continuous(name='Enrichment z-score', breaks=c(0, 5, 10)) +
     scale_fill_manual(values=color_by_cluster_top(cluster_density)) +
     facet_wrap(~ V1, ncol=1, scales='free') +
     theme_nature +
     theme(axis.text.x=element_blank()))
}

plot_enhancer_enrichments <- function(filename, cluster_density) {
    enrichments <- read.delim(filename, header=FALSE, sep=' ')
    enrichments$V1 <- toupper(sub('^[^-]*-', '', enrichments$V1))
    enrichments$V3 <- factor(enrichments$V3, levels=row.names(cluster_density))
    Cairo(type='pdf', file=sub('.in$', '.pdf', filename), width=190, height=100, units='mm')
    print(enrichment_by_cluster(enrichments, cluster_density))
    dev.off()
}
