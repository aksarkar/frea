requireNamespace('ggplot2')

#' @include theme_nature.R

plot_rank_correlation <- function(correlation_file) {
    correlations <- read.table(correlation_file)
    p <- (ggplot(correlations, aes(x=V1, y=V3, color=factor(V2))) +
          labs(xlab='Top n SNPs (full meta-analysis)', ylab='Hold-out Pearson correlation') +
          geom_line(size=I(.25)) +
          geom_hline(yintercept=0, color="black", size=I(.1), shape='dashed') +
          scale_x_continuous(limits=c(0, 1e5)) +
          scale_color_brewer(palette="Dark2", name="Cohort") +
          theme_nature +
          theme(legend.position="right"))

    Cairo(file=sub('.txt$', '.pdf', args[1]), type='pdf', height=50, width=89, unit='mm')
    print(p)
    dev.off()
}
