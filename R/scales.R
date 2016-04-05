requireNamespace('ggplot2')
requireNamespace('plyr')
requireNamespace('reshape2')

data(roadmap_sample_info)

color_by_eid <- with(roadmap_sample_info, scale_color_manual(values=setNames(as.character(color), EID)))

fill_by_eid <- with(roadmap_sample_info, scale_fill_manual(values=setNames(as.character(color), EID)))

#' Color scheme by cluster majority
#'
#' Threshold cluster densities to get subset of reference epigenomes, then take
#' majority color (tissue group)
color_by_cluster_majority <- function(cluster_density, threshold=0.25) {
    long_form_density <- subset(melt(cluster_density), value > threshold)
    annotated <- merge(long_form_density, roadmap_sample_info, by.x='Var2', by.y='EID')
    majority_color <- ddply(annotated, .(Var1), function (x) {names(which.max(table(x$color)))})
    with(majority_color, setNames(as.character(V1), Var1))
}

color_by_cluster_top <- function(cluster_density) {
    constitutive <- subset(ddply(subset(melt(cluster_density), value > .25), .(Var1), nrow), V1 > .5 * dim(cluster_density)[2])$Var1
    top_cell_by_cluster <- ddply(melt(cluster_density), .(Var1),
                                 function (x) {
                                     if (x$Var1[1] %in% constitutive) {
                                         data.frame(Var1=x$Var1[1], Var2='E129')
                                     }
                                     else {
                                         x[which.max(x$value),]
                                     }})
    annotated <- merge(top_cell_by_cluster, roadmap_sample_info, by.x='Var2', by.y='EID')
    with(annotated, setNames(as.character(color), Var1))
}

eid_ordering <- with(roadmap_sample_info, EID[order(position)])

roadmap_tissue_info <-
    ddply(roadmap_sample_info, .(group_name),
          function (x) {
              data.frame(tissue=tail(strsplit(as.character(x[1,]$group_name), ' & ')[[1]], 1),
                         position=min(x$position),
                         color=names(which.max(table(x$color))))
          })

tissue_ordering <- with(roadmap_tissue_info, tissue[order(position)])

color_by_tissue <- with(roadmap_tissue_info, setNames(as.character(color), tissue))

roadmap_tissue_legend <- function(theme_args, guide_args) {
    roadmap_tissue_info$tissue <- factor(roadmap_tissue_info$tissue, levels=tissue_ordering)
    dummy_plot <- (ggplot(roadmap_tissue_info, aes(x=tissue, fill=tissue)) +
                   geom_bar() +
                   scale_fill_manual(name='Tissue', values=color_by_tissue,
                                     guide=do.call(guide_legend, guide_args)) +
                   theme_nature +
                   do.call(theme, theme_args))
    my_gtable <- ggplotGrob(dummy_plot)
    my_gtable$grobs[grepl('guide-box', my_gtable$layout$name)][[1]]
}

phenotype_ordering <- c('CD', 'RA', 'T1D', 'AD', 'BIP', 'SCZ', 'CAD', 'T2D')

phenotype_legend <- function(theme_args, guide_args) {
    pheno <- factor(phenotype_ordering, ordered=TRUE)
    dummy_plot <- (ggplot(data.frame(x=pheno, y=1), aes(x, y, fill=pheno)) +
                   geom_bar(stat='identity') +
                   scale_fill_brewer(name='Phenotype', palette='Dark2',
                                     guide=do.call(guide_legend, guide_args)) +
                   theme_nature +
                   do.call(theme, theme_args))
    my_gtable <- ggplotGrob(dummy_plot)
    my_gtable$grobs[grepl('guide-box', my_gtable$layout$name)][[1]]
}

