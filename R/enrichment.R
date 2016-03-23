requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('grid')
requireNamespace('gtable')

#' @include heatmap.R scales.R theme_nature.R
#' @importFrom grid unit
NULL

enrichment_by_cluster <- function(enrichments, cluster_density, scale_name) {
    (heatmap(ggplot(enrichments, aes(x=cluster, y=pheno, fill=val))) +
     scale_heatmap(name=scale_name) +
     scale_y_discrete(levels=rev(levels(enrichments$pheno))) +
     labs(x='Enhancer module', y='Phenotype') +
     theme_nature +
     theme(axis.text.x=element_blank(),
           legend.position='right'))
}

parse_enhancer_enrichments <- function(filename) {
    enrichments <- read.delim(filename, header=FALSE, sep=' ')
    enrichments$pheno <- factor(toupper(sub('^[^-]*-', '', enrichments$V1)))
    subset(enrichments, V6 > 0 & V7 > 0)
}

plot_enhancer_enrichments_by_eid <- function(filename) {
    data(roadmap_sample_info)
    enrichments <- parse_enhancer_enrichments(filename)
    enrichments$eid <- factor(enrichments$V3, levels=eid_ordering)
    enrichments <- transform(enrichments, val=(V5 - V6)/ sqrt(V7))

    my_ggplot <- (heatmap(ggplot(enrichments, aes(x=eid, y=pheno, fill=val))) +
                  scale_heatmap(name='z-score') +
                  scale_y_discrete(limits=rev(levels(enrichments$pheno))) +
                  labs(x='Reference epigenome', y='Phenotype') +
                  theme(axis.text.x=element_blank(),
                        legend.position='bottom'))
    my_gtable <- ggplotGrob(my_ggplot)

    ## Add tissue colors
    tissue_grob <- ggplotGrob(epigenome_by_tissue() + coord_flip())$grobs[[4]]
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(2, 'mm'), 0)
    my_gtable <- gtable::gtable_add_grob(my_gtable, tissue_grob, t=1, l=4, r=4)

    ## Add tissue legend
    my_legend <- roadmap_tissue_legend(list(legend.position='bottom'), list(direction='horizontal', nrow=2))
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(10, 'mm'))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, t=dim(my_gtable)[1], l=1, r=dim(my_gtable)[2])

    Cairo(type='pdf', file=sub('.in$', '.pdf', filename), width=210, height=50, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}

plot_enhancer_enrichments <- function(filename, cluster_density, plot_log_fold=FALSE, flip=FALSE) {
    enrichments <- parse_enhancer_enrichments(filename)
    enrichments$cluster <- factor(enrichments$V3, levels=row.names(cluster_density))

    if (plot_log_fold) {
        enrichments <- transform(enrichments, val=log10(V5 / V6))
        scale_name <- "Log fold enrichment"
    }
    else {
        enrichments <- transform(enrichments, val=(V5 - V6)/ sqrt(V7))
        scale_name <- "z-score"
    }
    my_density <- (density_by_cluster(cluster_density, keep=unique(enrichments$cluster)) +
                   theme(legend.position='right',
                         axis.title.x=element_blank(),
                         axis.text.y=element_text(margin=margin(2))))
    my_gtable <- gtable:::rbind.gtable(ggplotGrob(my_density),
                                       ggplotGrob(enrichment_by_cluster(enrichments, cluster_density, scale_name)),
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

    if (flip) {
        w <- 270
        h <- 210
    }
    else {
        w <- 210
        h <- 270
    }

    Cairo(type='pdf', file=sub('.in$', '.pdf', filename), width=w, height=h, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}
