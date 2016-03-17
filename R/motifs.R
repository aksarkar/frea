requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('gdata')
requireNamespace('gtable')
requireNamespace('reshape2')

#' @include heatmap.R scales.R theme_nature.R

master_regulator_by_cofactor <- function(cofactors) {
    (sparse_heatmap(ggplot(cofactors, aes(x=tf, y=co_tf, fill=V5))) +
     scale_fill_gradient(name='Motif score', low='white', high='black') +
     scale_x_discrete(drop=FALSE) +
     labs(x='Master regulator', y='Cofactor') +
     theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)))
}

master_regulator_by_phenotype <- function(enrichments) {
    (sparse_heatmap(ggplot(enrichments, aes(x=tf, y=pheno, fill=log10(V5)))) +
     scale_heatmap(name='Log odds ratio', limits=c(0, max(log10(enrichments$V5))), breaks=seq(0, 2)) +
     scale_y_discrete(limits=rev(levels(enrichments$pheno))) +
     labs(x='Master regulator', y='Phenotype') +
     theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1),
           legend.position='bottom'))
}

plot_motif_enrichments <- function(enrichments_file, cofactors_file) {
    data(honeybadger_cluster_density)

    enrichments <- read.delim(gzfile(enrichments_file), header=FALSE, sep=' ')
    enrichments$pheno <- toupper(sub('^[^-]*-', '', enrichments$V1))
    enrichments$tf <- sub('_.*$', '', enrichments$V4)
    enrichment_matrix <- acast(enrichments, formula=pheno ~ tf, value.var='V5', fill=0, fun.aggregate=max)
    row_order <- hclust(dist(enrichment_matrix))$order
    enrichments$pheno <- with(enrichments, factor(pheno, levels=row.names(enrichment_matrix)[order(row_order)]))
    col_order <- apply(enrichment_matrix, 2, function (x) {row_order[which.max(x)]})
    enrichments$tf <- factor(enrichments$tf, levels=names(col_order[order(col_order)]))

    cofactors <- read.delim(gzfile(cofactors_file), header=FALSE, sep=' ')
    cofactors$pheno <- toupper(sub('^[^-]*-', '', cofactors$V1))
    cofactors$tf <- factor(sub('_.*$', '', cofactors$V3), levels=levels(enrichments$tf))
    cofactors$co_tf <- sub('_.*$', '', cofactors$V4)
    match_matrix <- acast(cofactors, formula=tf ~ co_tf, value.var='V5', fill=0, fun.aggregate=max)
    col_order <- apply(match_matrix, 2, function (x) {row_order[which.max(x)]})
    cofactors$co_tf <- factor(cofactors$co_tf, levels=names(col_order[order(col_order)]))

    my_grobs <- lapply(list(master_regulator_by_phenotype(enrichments),
                            master_regulator_by_cofactor(cofactors)),
                       ggplotGrob)
    my_gtable <- do.call(gtable:::rbind.gtable, c(my_grobs, size='last'))

    Cairo(type='pdf', file=sub('.txt.gz$', '-regulators.pdf', enrichments_file), width=190, height=700, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}
