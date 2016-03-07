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
    (heatmap(ggplot(eids, aes(x=EID, y=rep(1), fill=EID))) +
     scale_y_discrete() +
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
    (ggplot(data.frame(x=keep, y=rep(1)), aes(x, y, fill=factor(x))) +
     scale_x_discrete(limits=row.names(cluster_density)) +
     scale_y_discrete() +
     color_by_cluster_top(cluster_density=cluster_density) +
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
     scale_fill_gradient(low='white', high='black') +
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

expression_by_tf <- function(expression) {
    expression_by_tf <- ddply(expression, .(eid, tf), summarize, rpkm=mean(rpkm))
    my_expression <- ddply(expression_by_tf, .(tf), summarize, eid=eid, rpkm=rpkm / max(rpkm))
    (heatmap(ggplot(my_expression, aes(x=eid, y=tf, fill=rpkm))) +
     scale_fill_gradient(name='Relative expression', low='white', high='black') +
     labs(x='Reference epigenome') +
     theme_nature +
     theme(legend.position='bottom',
           legend.key.height=grid::unit(2.5, 'mm'),
           legend.key.width=grid::unit(5, 'mm'),
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
           legend.key.height=grid::unit(2.5, 'mm'),
           legend.key.width=grid::unit(5, 'mm'),
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
    my_grobs <- lapply(list(pheno_by_tf(enriched_tf_expr),
                            expression_by_tf(enriched_tf_expr),
                            epigenome_by_tissue(unique(enriched_tf_expr$eid))),
                       ggplotGrob)
    my_gtable <- gtable:::cbind_gtable(my_grobs[[1]], my_grobs[[2]], size='last')

    # Fix up the widths/heights since these get clobbered by bind
    my_gtable$widths <- grid:::unit.list(my_gtable$widths)
    my_gtable$heights <- grid:::unit.list(my_gtable$heights)
    widths <- lapply(list(enriched_tf_expr$pheno, enriched_tf_expr$eid), #
                     function(x) grid::unit(length(unique(x)), 'null'))
    height <- grid::unit(length(unique(enriched_tf_expr$tf)), 'null')
    my_gtable$widths[my_gtable$layout$l[grepl("panel", my_gtable$layout$name)]] <- widths
    my_gtable$heights[my_gtable$layout$t[grepl("panel", my_gtable$layout$name)]] <- list(height)

    # Add the tissue colors grob
    my_gtable <- gtable::gtable_add_rows(my_gtable, grid::unit(8, 'mm'), 1)
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_grobs[[3]]$grobs[[4]], 2, 9, name='tissues')

    Cairo(type='pdf', file=sub('.txt.gz$', '.pdf', enrichment_file), width=190, height=100, units='mm')
    grid::grid.draw(my_gtable)
    dev.off()
}
