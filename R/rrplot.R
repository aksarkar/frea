requireNamespace('directlabels')
requireNamespace('gtable')
requireNamespace('grid')
requireNamespace('ggplot2')

#' @include directlabels.R scales.R theme_nature.R
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

rrplot <- function(X, total_cutoff, axis_labels, rank_cutoff) {
    X <- subset(X, total <= total_cutoff)
    direct_labels <- geom_dl(ggplot2::aes(label=eid),
                             method=list(cex=7/16,
                                         last.points,
                                         calc.boxes,
                                         dl.trans(x = x + .1),
                                         top(10),
                                         'bumpdown'))
    (ggplot(X, aes(x=total, y=y, color=factor(eid))) +
     geom_line(size=I(.35 / ggplot2:::.pt)) +
     geom_hline(yintercept=0, color='black', size=I(.5 / ggplot2:::.pt)) +
     geom_vline(data=rank_cutoff, aes(xintercept=xintercept), color='red', size=I(.5 / ggplot2:::.pt)) +
     direct_labels +
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
    rank_cutoff <- rrplot_elbow(cumulative_deviation)
    write.table(x=rank_cutoff, file=sub('.txt.gz$', '.ranks', counts_file),
                quote=FALSE, row.names=FALSE, col.names=FALSE)

    my_gtable <- ggplotGrob(rrplot(cumulative_deviation, cutoff,
                                   labs(x=xlab, y='Cumulative deviation'),
                                   rank_cutoff))

    ## Don't clip direct labels
    my_gtable$layout$clip[grepl('panel', my_gtable$layout$name)] <- 'off'

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
