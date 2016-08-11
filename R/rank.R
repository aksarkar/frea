requireNamespace('ggplot2')
requireNamespace('dplyr')

#' @include theme_nature.R

plot_rank_correlation <- function(correlation_file) {
    correlations <- read.table(correlation_file)
    correlations$V2 <- factor(correlations$V2,
                              levels=c('WTCCC', 'NARAC1', 'NARAC2', 'EIRA',
                                       'CANADA', 'BRASS', 'Overall'))
    p <- (ggplot(correlations, aes(x=V1, y=V3, color=factor(V2))) +
          geom_line(size=I(.25)) +
          geom_hline(yintercept=0, color="black", size=I(.1), linetype='dashed') +
          scale_x_continuous(limits=c(0, 1e5), expand=c(0, 0), name='Top n SNPs (hold-out absolute z-score)') +
          scale_y_continuous(limits=c(0, 1), name='Hold-out Pearson correlation') +
          scale_color_brewer(palette="Dark2", name="Cohort") +
          theme_nature +
          theme(legend.position="right"))

    Cairo(file=sub('.txt$', '.pdf', correlation_file), type='pdf', height=50, width=89, unit='mm')
    print(p)
    dev.off()
}

qvalue_at_threshold <- function(marginal_file, ranks_file) {
    marginal <- read.delim(gzfile(marginal_file), header=FALSE)
    ranks <- read.table(ranks_file, sep=' ')
    (marginal %>%
     group_by(V1) %>%
     mutate(p = 10 ** -V3) %>%
     mutate(q = qvalue::qvalue(p)$qvalues) %>%
     select(V1, q) %>%
     arrange(q) %>%
     slice(ranks[ranks$V1 == V1[1],]$V3))
}
