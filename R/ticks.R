#' Draw ticks at increasing p-value thresholds
logp.ticks <- function(Y) {
  T <- rev(table(cut(Y$V5, labels=seq(1, 10), breaks=seq(0, 10), right=FALSE)))
  Z <- data.frame(thresh=names(T), rank=cumsum(T), row.names=NULL)
  geom_vline(aes(xintercept=rank), data=Z, linetype='dashed', size=I(.35 / ggplot2:::.pt))
}

#' Draw ticks at increasing chi-square thresholds
#'
#' These correspond to p-value thresholds in powers of 10
chisq.ticks <- function(Y) {
  T <- rev(table(cut(Y$V5, labels=seq(1, 10), breaks=qchisq(1 - 10 ** -seq(0, 10), 1), right=FALSE)))
  Z <- data.frame(thresh=names(T), rank=cumsum(T), row.names=NULL)
  ## write.table(Z, file='ticks', col.names=FALSE, row.names=FALSE, quote=FALSEn)
  geom_vline(aes(xintercept=rank), data=Z, linetype='dashed', size=I(.35 / ggplot2:::.pt))
}

#' Draw ticks at pre-computed thresholds
id.ticks <- function(Y) {
  colnames(Y) <- c('thresh', 'rank')
  geom_vline(aes(xintercept=rank), data=Y, linetype='dashed', size=I(.35 / ggplot2:::.pt))
}

