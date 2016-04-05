requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('gdata')
requireNamespace('grid')
requireNamespace('gridExtra')
requireNamespace('gtable')
requireNamespace('reshape2')

#' @importFrom grid unit
#' @include heatmap.R scales.R theme_nature.R

expression_by_tf <- function(expression) {
    expression_by_tf <- ddply(expression, .(eid, tf), summarize, rpkm=mean(rpkm))
    my_expression <- ddply(expression_by_tf, .(tf), summarize, eid=eid, rpkm=rpkm / max(rpkm))
    (heatmap(ggplot(my_expression, aes(x=eid, y=tf, fill=rpkm))) +
     scale_fill_gradient(name='Relative expression', low='white', high='black') +
     labs(x='Reference epigenome') +
     theme_nature +
     theme(legend.position='bottom',
           legend.key.height=unit(2.5, 'mm'),
           legend.key.width=unit(5, 'mm'),
           axis.text.x=element_text(angle=90, hjust=1, vjust=0.5),
           axis.text.y=element_blank(),
           axis.title.y=element_blank()))
}

pheno_by_tf <- function(expression) {
    my_expression <- ddply(expression, .(pheno, tf), summarize, odds=log10(unique(odds)))
    (heatmap(ggplot(my_expression, aes(x=pheno, y=tf, fill=odds))) +
     labs(x='Phenotype', y='Transcription factor') +
     scale_heatmap(name='Log odds ratio') +
     theme(legend.position='bottom',
           legend.key.height=unit(2.5, 'mm'),
           legend.key.width=unit(5, 'mm'),
           axis.text.x=element_text(angle=90, hjust=1, vjust=0.5)))
}

plot_tf_expression <- function(enrichment_file, ensembl_file, sample_file, rpkm_file) {
    data(honeybadger_cluster_density)
    # Parse the data
    constitutive <- row.names(honeybadger_cluster_density)[1:3]
    enrichments <- subset(read.delim(gzfile(enrichment_file), header=FALSE, sep=' '), V2 %in% constitutive)
    enrichments$tf <- sub("_.*$", "", enrichments$V4)
    enrichments$pheno <- toupper(sub("[^-]*-", "", enrichments$V1))
    enrichments <- ddply(enrichments, .(pheno, tf), function (x){data.frame(odds=max(x$V5))})

    samples <- read.delim(sample_file, header=FALSE)
    genes <- read.delim(ensembl_file, header=FALSE)[c('V1', 'V7')]
    rpkm_by_eid <- read.delim(gzfile(rpkm_file), header=FALSE, skip=1)
    colnames(rpkm_by_eid) <- c('id', as.character(samples$V1))

    enrichment_by_ensembl_id <- merge(enrichments, genes, by.x='tf', by.y='V7')
    enriched_tf_expr <- melt(merge(enrichment_by_ensembl_id, rpkm_by_eid, by.x='V1', by.y='id'),
                             id.vars=c('pheno', 'tf', 'odds', 'V1'),
                             variable.name='eid',
                             value.name='rpkm')
    enriched_tf_expr$eid <- factor(enriched_tf_expr$eid, levels=eid_ordering)
    enriched_tf_expr <- subset(enriched_tf_expr, !is.na(eid))
    best <- ddply(enriched_tf_expr, .(tf), function (x) {x$eid[which.max(x$rpkm)]})
    enriched_tf_expr$tf <- factor(enriched_tf_expr$tf, levels=rev(with(best, tf[order(V1)])))

    # Bind the heatmap grobs
    best_tissue_by_tf <- (heatmap(ggplot(ddply(enriched_tf_expr, .(tf),
                                               function (x) {
                                                   as.character(x$eid[which.max(x$rpkm)])
                                               }),
                                         aes(x=rep(1), y=tf, fill=V1)) +
                                  fill_by_eid +
                                  theme(axis.ticks=element_blank(),
                                        axis.title=element_blank(),
                                        axis.text=element_blank(),
                                        axis.line=element_blank())))
    my_grobs <- lapply(list(pheno_by_tf(enriched_tf_expr),
                            expression_by_tf(enriched_tf_expr),
                            epigenome_by_tissue(unique(enriched_tf_expr$eid)),
                            best_tissue_by_tf),
                       ggplotGrob)
    my_gtable <- gtable:::cbind_gtable(my_grobs[[1]], my_grobs[[2]], size='last')

    # Fix up the widths/heights since these get clobbered by bind
    my_gtable$widths <- grid:::unit.list(my_gtable$widths)
    my_gtable$heights <- grid:::unit.list(my_gtable$heights)
    widths <- lapply(list(enriched_tf_expr$pheno, enriched_tf_expr$eid), #
                     function(x) unit(length(unique(x)), 'null'))
    height <- unit(length(unique(enriched_tf_expr$tf)), 'null')
    my_gtable$widths[my_gtable$layout$l[grepl("panel", my_gtable$layout$name)]] <- widths
    my_gtable$heights[my_gtable$layout$t[grepl("panel", my_gtable$layout$name)]] <- list(height)

    # Add the tissue colors grob
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(3, 'mm'), 1)
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_grobs[[3]]$grobs[[4]], 2, 9, name='tissues')
    my_gtable <- gtable::gtable_add_cols(my_gtable, unit(3, 'mm'), 6)
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_grobs[[4]]$grobs[[4]], t=4, l=7, name='best')
    
    ## Add the tissue legend
    my_legend <- roadmap_tissue_legend(theme_args=list(legend.position='bottom'),
                                       guide_args=list(direction='horizontal', nrow=2))
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(1, 'lines'))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, t=-1, l=1, r=-1)

    Cairo(type='pdf', file=sub('.txt.gz$', '-expression.pdf', enrichment_file), width=190, height=120, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}
