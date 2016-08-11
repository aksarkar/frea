requireNamespace('gtable')
requireNamespace('grid')
requireNamespace('ggplot2')

#' @include scales.R theme_nature.R
#' @importFrom grid unit

rrplot_elbow <- function(X) {
    inflections <- ddply(X, .(study, phenotype, feature, eid),
                         function(x) {
                             hull <- x[order(chull(x$total, x$y)),]
                             ddy <- diff(hull$y, differences=2)
                             hull[which(diff(sign(ddy)) != 0)[1] + 1,]$total
                         })
    rank_cutoff <- ddply(inflections, .(study, phenotype), summarize, xintercept=max(V1, na.rm=TRUE))
}

rrplot_highlights <- function() {
    data(roadmap_sample_info)
    t_cell <- subset(roadmap_sample_info, group_name == 'Blood & T-cell')$EID
    b_cell <- subset(roadmap_sample_info, group_name == 'HSC & B-cell')$EID
    brain <- subset(roadmap_sample_info, group_name == 'Brain')$EID
    psych <- unlist(lapply(list(t_cell, brain), as.character))
    list('CD' = t_cell,
         'RA' = t_cell,
         'T1D' = t_cell,
         'AD' = b_cell,
         'BIP' = psych,
         'SCZ' = psych,
         'CAD' = 'E075',
         'T2D' = c('E087', 'E109'))
}

rrplot <- function(X, total_cutoff, axis_labels, rank_cutoff, highlight=FALSE) {
    X <- subset(X, total <= total_cutoff)
    highlights <- rrplot_highlights()
    if (highlight) {
        X <- ddply(X, .(phenotype),
                   function(x) {
                       x$alpha = with(x, ifelse(eid %in% unlist(highlights[[x$phenotype[1]]]), 1, 0.9))
                       x
                   })
    }
    else {
        X$alpha = rep(1)
    }
    (ggplot(X, aes(x=total, y=y, color=factor(eid))) +
     geom_line(aes(alpha=alpha), size=I(.35 / ggplot2:::.pt)) +
     geom_hline(yintercept=0, color='black', size=I(.5 / ggplot2:::.pt)) +
     geom_vline(data=rank_cutoff, aes(xintercept=xintercept), color='red', size=I(.5 / ggplot2:::.pt)) +
     scale_x_continuous(labels=comma, limits=c(0, total_cutoff),
                        breaks=seq(0, total_cutoff, total_cutoff / 4),
                        expand=c(0, 0)) +
     scale_y_continuous(expand=c(0, 0)) +
     facet_wrap(~ phenotype, ncol=4, scales='free') +
     color_by_eid +
     axis_labels +
     theme_nature +
     theme(legend.position='none',
           panel.margin.x=unit(7, 'mm'),
           panel.margin.y=unit(4, 'mm'),
           plot.margin=unit(c(0, 7, 2, 0), 'mm')))
}

plot_rrplot <- function(counts_file, xlab, cutoff) {
    X <- read.csv(gzfile(counts_file), header=FALSE)
    colnames(X) <- c('total', 'phenotype', 'eid', 'feature', 'count', 'expected')
    cumulative_deviation <- ddply(X, .(phenotype, feature, eid), transform,
                                  study=phenotype,
                                  phenotype=toupper(sub("^.*-", "", phenotype)),
                                  y=(count - expected) / max(count))
    cumulative_deviation$phenotype <- factor(cumulative_deviation$phenotype, levels=phenotype_ordering)
    rank_cutoff <- rrplot_elbow(cumulative_deviation)
    write.table(x=rank_cutoff, file=sub('.txt.gz$', '.ranks', counts_file),
                quote=FALSE, row.names=FALSE, col.names=FALSE)

    my_gtable <- ggplotGrob(rrplot(cumulative_deviation, cutoff,
                                   labs(x=xlab, y='Cumulative deviation'),
                                   rank_cutoff))

    ## Add the legend
    my_gtable <- gtable::gtable_add_rows(my_gtable, unit(1, 'lines'))
    my_legend <- roadmap_tissue_legend(theme_args=list(legend.position='bottom'),
                                       guide_args=list(direction='horizontal', nrow=2))
    my_gtable <- gtable::gtable_add_grob(my_gtable, my_legend, name='legend',
                                         t=dim(my_gtable)[1], l=1, r=dim(my_gtable)[2])

    Cairo(file=sub('.txt.gz$', ".pdf", counts_file), type='pdf', width=190,
          height=80, units='mm', family='Helvetica')
    grid::grid.draw(my_gtable)
    dev.off()
}
