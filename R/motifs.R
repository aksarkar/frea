requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('gdata')
requireNamespace('gtable')
requireNamespace('reshape2')

#' @include heatmap.R scales.R theme_nature.R

master_regulator_by_phenotype <- function(enrichments) {
    (heatmap(ggplot(enrichments, aes(x=tf, y=pheno, fill=log10(V5)))) +
     scale_heatmap(name='Log odds ratio', limits=c(0, max(log10(enrichments$V5))), breaks=seq(0, 2)) +
     scale_y_discrete(limits=rev(levels(enrichments$pheno))) +
     labs(x='Master regulator', y='Phenotype') +
     theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
           legend.position='bottom'))
}

parse_motif_enrichments <- function(enrichments_file) {
    data(honeybadger_cluster_density)
    enrichments <- read.delim(gzfile(enrichments_file), header=FALSE, sep=' ')
    enrichments$pheno <- toupper(sub('^[^-]*-', '', enrichments$V1))
    enrichments$tf <- sub('_.*$', '', enrichments$V4)
    enrichments
}

plot_motif_enrichments <- function(enrichments_file) {
    enrichments <- parse_motif_enrichments(enrichments_file)
    enrichment_matrix <- acast(enrichments, formula=pheno ~ tf, value.var='V5', fill=0, fun.aggregate=max)
    row_order <- hclust(dist(enrichment_matrix))$order
    enrichments$pheno <- with(enrichments, factor(pheno, levels=row.names(enrichment_matrix)[order(row_order)]))
    col_order <- apply(enrichment_matrix, 2, function (x) {row_order[which.max(x)]})
    enrichments$tf <- factor(enrichments$tf, levels=names(col_order[order(col_order)]))

    Cairo(type='pdf', file=sub('.txt.gz$', '-regulators.pdf', enrichments_file), width=190, height=50, units='mm')
    print(master_regulator_by_phenotype(enrichments))
    dev.off()
}

plot_master_regulator_counts <- function(enrichments_file, cluster_density) {
    enrichments <- parse_motif_enrichments(enrichments_file)
    enrichments$cluster <- factor(enrichments$V2, levels=rownames(cluster_density))
    master_regulators_by_cluster <- ddply(enrichments, .(cluster),
                                          function(x) {data.frame(count=length(unique(x$tf)))})
    my_plot <- (ggplot(master_regulators_by_cluster,
                       aes(x=cluster, y=count, fill=cluster)) +
                geom_bar(stat='identity') +
                labs(y='Count of master regulators') +
                scale_fill_manual(values=color_by_cluster_top(cluster_density)) +
                theme_nature +
                theme(axis.text.x=element_blank(),
                      axis.title.x=element_blank(),
                      plot.margin=unit(rep(2, 4), 'mm')))
    my_grobs <- lapply(list(my_plot, density_by_cluster(cluster_density, keep=unique(enrichments$cluster))),
                       ggplotGrob)
    ## Add placeholder for top legend so rbind doesn't complain
    my_grobs[[1]] <- gtable::gtable_add_cols(my_grobs[[1]], unit(1, 'lines'))
    my_gtable <- do.call(gtable:::rbind.gtable, c(my_grobs, size='last'))

    ## Fix up widths/heights
    my_gtable$widths <- grid:::unit.list(my_gtable$widths)
    my_gtable$heights <- grid:::unit.list(my_gtable$heights)
    widths <- list(unit(length(unique(enrichments$cluster)), 'null'))
    heights <- lapply(list(10, nrow(roadmap_tissue_info)), unit, 'null')
    my_gtable$widths[my_gtable$layout$l[grepl("panel", my_gtable$layout$name)]] <- widths
    my_gtable$heights[my_gtable$layout$t[grepl("panel", my_gtable$layout$name)]] <- heights

    ## Add tissue legend
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(1, 'lines'))
    my_legend <- roadmap_tissue_legend(list(legend.position='bottom'), list(direction='horizontal', nrow=2))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, t=-1, l=1, r=-1)

    Cairo(type='pdf', file=sub('.txt.gz$', '-counts.pdf', enrichments_file), width=190, height=120, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}

