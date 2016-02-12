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
    annotated <- merge(long_form_density, roadmap_sample_info, by.x="Var2", by.y="EID")
    majority_color <- ddply(annotated, .(Var1), function (x) {names(which.max(table(x$color)))})
    with(majority_color, scale_color_manual(values=setNames(as.character(V1), Var1)))
}

color_by_cluster_top <- function(cluster_density) {
    top_cell_by_cluster <- ddply(melt(honeybadger_cluster_density), .(Var1), function (x) {x[which.max(x$value),]})
    annotated <- merge(top_cell_by_cluster, roadmap_sample_info, by.x="Var2", by.y="EID")
    with(annotated, scale_color_manual(values=setNames(as.character(color), Var1)))
}

eid_ordering <- with(roadmap_sample_info, EID[order(position)])

roadmap_tissue_info <-
    ddply(roadmap_sample_info, .(group_name),
          function (x) {
              data.frame(tissue=tail(strsplit(as.character(x[1,]$group_name), " & ")[[1]], 1),
                         position=min(x$position),
                         color=names(which.max(table(x$color))))
          })

tissue_ordering <- with(roadmap_tissue_info, tissue[order(position)])

color_by_tissue <- with(roadmap_tissue_info, setNames(as.character(color), tissue))
