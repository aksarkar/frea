#' QQ plot
#'
#' @param data Data frame. Required columns are observed and expected
#' @param theoretical Theoretical distribution
#' @param subsample Only draw every nth point along the curve
qqplot <- function(data, theoretical, subsample=1000) {
  # Expected QQ line under the null
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
   geom_abline(slope=m, yintercept=b) +
   theme_nature)
}
