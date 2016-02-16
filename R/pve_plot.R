requireNamespace('Cairo')
requireNamespace('ggplot2')

#' @include scales.R theme_nature.R

pve_by_annotation <- function (h2g, aes) {
    (ggplot(h2g, aes) +
     geom_bar(stat='identity') +
     scale_x_discrete(drop=FALSE) +
     scale_y_continuous(name='Average PVE per SNP', labels=comma) +
     facet_wrap(~ V1, ncol=1, scales='free') +
     theme_nature)
}

pve_by_cluster <- function(h2g) {
    h2g$V2 <- factor(h2g$V2, levels=honeybadger_ordering)
    (pve_by_annotation(h2g, aes(x=V2, y=V3/V6, fill=V2)) +
     scale_fill_manual(values=color_by_cluster_majority(honeybadger_cluster_density)) +
     labs(x='Enhancer module'))
}

pve_by_cluster_tissue_group <- function (h2g) {
    h2g <- transform(h2g, tissue=sub('_', ' ', V2))
    h2g$tissue <- factor(h2g$tissue, levels=tissue_ordering)
    (pve_by_annotation(h2g, aes(x=tissue, y=V3/V6, fill=tissue)) +
     scale_fill_manual(values=color_by_tissue) +
     labs(x='Tissue group') +
     theme(axis.text.x=element_text(angle=45,hjust=1)))
}
