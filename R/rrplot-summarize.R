#' Compute cumulative deviations
#' 
#' @param X Data frame. Required columns are phenotype, feature, celltype,
#' count, expected
dev <- function(X) {
  ddply(X, .(phenotype, feature, celltype), transform,
        y=(count - expected) / max(count))
}

#' Compute cumulative fraction recovered.
#' 
#' This is the basic idea behind IDR; however, this implementation is much
#' simpler than their definition of correspondence curves
#'
#' @param X Data frame. Required columns are phenotype, feature, celltype,
#' count, expected
cumulative.fraction <- function(X) {
  ddply(X, .(phenotype, feature, celltype), transform,
        y=count / max(count))
}

#' Look for peaks
#'
#' This is one way to find the elbow of the curve, i.e. the cutoff beyond which
#' we stop seeing enrichment
#'
#' @param series Sequence of scores
#' @param span Window size
#' @param ties.method 
peaks <- function(series, span=3, ties.method='first') {
  stopifnot((span <- as.integer(span)) %% 2 == 1)
  z <- embed(series, span)
  s <- span %/% 2
  v <- max.col(z, ties.method=ties.method) == 1 + s
  pad <- rep(FALSE, s)
  result <- c(pad, v, pad)
  result
}
