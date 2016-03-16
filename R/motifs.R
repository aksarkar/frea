requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('gdata')
requireNamespace('gtable')
requireNamespace('reshape2')

#' @include heatmap.R scales.R theme_nature.R

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

