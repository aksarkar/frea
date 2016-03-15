requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('grid')
requireNamespace('gtable')

#' @include heatmap.R scales.R theme_nature.R
#' @importFrom grid unit
NULL

enrichment_by_cluster <- function(enrichments, cluster_density) {
    (heatmap(ggplot(enrichments, aes(x=cluster, y=pheno, fill=(V5 - V6)/ sqrt(V7)))) +
     scale_heatmap(name='z-score') +
     labs(x='Enhancer module', y='Phenotype') +
     theme_nature +
     theme(axis.text.x=element_blank(),
           legend.position='right'))
}

plot_enhancer_enrichments <- function(filename, cluster_density) {
    enrichments <- read.delim(filename, header=FALSE, sep=' ')
    enrichments$pheno <- toupper(sub('^[^-]*-', '', enrichments$V1))
    enrichments$cluster <- factor(enrichments$V3, levels=row.names(cluster_density))

    my_density <- (density_by_cluster(cluster_density, keep=unique(enrichments$cluster)) +
                   theme(legend.position='right',
                         axis.title.x=element_blank(),
                         axis.text.y=element_text(margin=margin(2))))
    my_gtable <- gtable:::rbind.gtable(ggplotGrob(my_density),
                                       ggplotGrob(enrichment_by_cluster(enrichments, cluster_density)),
                                       size='last')
    ## Add tissue colors
    cluster_by_tissue_grob <- ggplotGrob(cluster_by_tissue(cluster_density, keep=unique(enrichments$cluster)))
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(4, 'mm'), 2)
    my_gtable <- gtable::gtable_add_grob(my_gtable, cluster_by_tissue_grob, t=3, l=4)
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(4, 'mm'), 7)
    my_gtable <- gtable::gtable_add_grob(my_gtable, cluster_by_tissue_grob, t=8, l=4)
    epigenome_by_tissue_grob <- ggplotGrob(epigenome_by_tissue())
    my_gtable <- gtable::gtable_add_cols(my_gtable, unit(2, 'mm'), 2)
    my_gtable <- gtable::gtable_add_grob(my_gtable, epigenome_by_tissue_grob$grobs[[4]], 4, 3)

    ## Add tissue legend
    my_gtable <- gtable::gtable_add_cols(my_gtable, unit(20, 'mm'), 0)
    my_legend <- roadmap_tissue_legend(list(legend.position='left'), list(direction='vertical'))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, 4, 1)

    Cairo(type='pdf', file=sub('.in$', '.pdf', filename), width=210, height=270, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}
