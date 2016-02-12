requireNamespace('Cairo')
requireNamespace('ggplot2')
requireNamespace('grid')
requireNamespace('gtable')

#' @include scales.R theme_nature.R
NULL

data(roadmap_sample_info)

scale_heatmap <- function(...) {
    scale_fill_gradient(low='#fee8c8', high='#e34a33', ...)
}

heatmap <- function(plot) {
    plot + geom_raster() + coord_fixed() + theme_nature
}

epigenome_by_tissue <- 
    (heatmap(ggplot(roadmap_sample_info, aes(x=EID, y=rep(1), fill=EID))) +
     scale_x_discrete(limits=eid_ordering) +
     scale_y_discrete() +
     fill_by_eid +
     theme(axis.ticks=element_blank(),
           axis.title=element_blank(),
           axis.text=element_blank(),
           axis.line=element_blank()))

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

enrichment_by_cluster <- function(enrichments, cluster_density) {
    enrichments$V2 <- factor(enrichments$V2, levels=row.names(cluster_density))
    enrichment_by_cluster_panel <-
        (heatmap(ggplot(enrichments, aes(x=V1, y=V2, fill=V3))) +
         scale_heatmap() +
         scale_x_discrete(name='Phenotype') +
         scale_y_discrete(name='Enhancer module') +
         coord_fixed() +
         theme(axis.text.x=element_blank()))
    density_by_cluster_panel <- density_by_cluster(cluster_density, enrichments$V2) + coord_flip()
    grobs <- lapply(list(density_by_cluster_panel, enrichment_by_cluster_panel), ggplotGrob)
    do.call(cbind, c(grobs, size='last'))
}

preview <- function() {
    data(honeybadger2_impute_p2_cluster_density)
    enrichments <- read.delim('/broad/compbio/aksarkar/projects/gwas/wtccc1/EC21/results/matched/hb2-imputed-features/sig', header=F, sep=' ')
    CairoPDF(file='/broad/hptmp/aksarkar/test.pdf', width=11, height=8.5)
    grid.draw(enrichment_by_cluster(enrichments, honeybadger2_impute_p2_cluster_density))
    dev.off()
}
