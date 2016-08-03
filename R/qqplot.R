#' QQ plot
#'
#' @param data Data frame. Required columns are observed and expected
#' @param theoretical Theoretical distribution
#' @param subsample Only draw every nth point along the curve
qqplot <- function(data, theoretical, subsample=1000) {
    ## Expected QQ line under the null
    quantiles <- c(.25, .75)
    y <- quantile(data$observed, quantiles)
    x <- theoretical(quantiles)
    m <- diff(y) / diff(x)
    b <- data$observed[1L] - m * data$expected[1L]
    if (subsample > 1) {
        subsampled <- data[seq(1, nrow(data), subsample),]
    }
    else {
        subsampled <- data
    }
    (ggplot(subsampled, aes(x=expected, y=observed)) +
     geom_path() +
     geom_abline(slope=m, intercept=b, linetype='dashed') +
     theme_nature)
}

plot_qqplot <- function(summary_file, ...) {
    marginal_stats <- read.delim(gzfile(summary_file), header=FALSE)
    marginal_stats$observed <- marginal_stats$V5
    marginal_stats <- marginal_stats[order(marginal_stats$observed),]
    marginal_stats$expected <- with(marginal_stats, -log10(rev(ppoints(observed))))
    Cairo(type='pdf', file=sub('.bed.gz', '.pdf', summary_file), width=80, height=80, units='mm')
    grid::grid.draw(qqplot(marginal_stats, theoretical=punif, ...))
    dev.off()
}
