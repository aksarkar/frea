requireNamespace('ggplot2')

rrplot.annotate <- function(X, cutoff) {
  Y <- X[X$total == cutoff, ]
  write.table(rev(Y[order(Y$y),]$celltype), sub('.in.gz$', '.txt', args[1]),
              col.names=FALSE, row.names=FALSE, quote=FALSE)
}

rrplot <- function(X, zero, cutoff, scales.labels, direct.labels, margin=7) {
  stopifnot(is.data.frame(X))
  stopifnot(is.numeric(zero))
  stopifnot(is.numeric(cutoff))
  stopifnot(cutoff > 0)
  X <- X[X$total <= cutoff, ]
  exponent <- ceiling(log10(max(abs(X$y))))
  (ggplot(X, aes(x=total, y=y, color=factor(celltype))) +
   geom_line(size=I(.35 / ggplot2:::.pt)) +
   geom_hline(yintercept=zero, color='black', size=I(.5 / ggplot2:::.pt)) +
   direct.labels +
   annotate('text', x=-Inf, y=Inf, size=(5 / ggplot2:::.pt), hjust=1, vjust=0, label=sprintf('(10^{%d})', floor(log10(max(abs(X$y))))), parse=TRUE) +
   scale_x_continuous(labels=comma, limits=c(0, cutoff),                      breaks=seq(0, cutoff, cutoff / 4),                      expand=c(0, 0)) +
   scale_y_continuous(labels=function(x) {sub('e.*', '', sprintf('%.1e', x))},                           expand=c(0, 0)) +
   scale_color_roadmap +
   scales.labels +
   theme_nature +
   theme(legend.position='none',
         axis.text=element_text(size=5),
         plot.margin=unit(c(7 / ggplot2:::.pt, margin, 0, 0), 'mm')))
}

rrplot.dev <- function(X, xlab, cutoff) {
  direct.labels <- geom_dl(aes(label=celltype),
                           method=list(cex=7/16,
                             last.points,
                             calc.boxes,
                             dl.trans(x = x + .1),
                             top(10),
                             'bumpdown'))
  scales.labels <- labs(x=xlab, y='Cumulative deviation')
  rrplot(dev(X), 0, as.numeric(cutoff), scales.labels, direct.labels)
}

rrplot.snps <- function(X, cutoff=100000) {
  rrplot.dev(X, 'SNP rank by p-value', cutoff)
}

rrplot.corr <- function(X) {
  scales.labels <- labs(x='WTCCC1 rank by p-value', y='Cumulative fraction recovered')
  (rrplot(cumulative.fraction(X), 0, signif(max(X$total), 2), scales.labels, NULL) +
   geom_abline(slope=1/max(X$total), yintercept=0, linetype='dotted', size=I(.35/ggplot2:::.pt), color='black') +
   guides(color=guide_legend('Threshold')) +
   theme(legend.position='right'))
}

rrplot.loci <- function(X, cutoff=10000) {
  rrplot.dev(X[seq(1, nrow(X), 10),], 'Independent loci by tag p-value', cutoff)
}

rrplot.example <- function(X) {
  scales.labels <- labs(x='SNP rank by p-value', y='Cumulative deviation')
  rrplot(dev(X), 0, 6.2e6, scales.labels, NULL)
}

rrplot.clusters <- function(X) {
  direct.labels <- geom_dl(aes(label=celltype),
                           method=list(cex=5/16,
                             last.points,
                             calc.boxes,
                             dl.trans(x = x + .1),
                             top(10),
                             'bumpdown'))
  Y <- dev(X)
  dY <- ddply(Y, .(phenotype, feature, celltype), transform, dy=c(0, diff(y))^2)
  Z <- cast(dY, celltype ~ total, value='dy')
  h <- hclust(dist(Z))
  groups <- cutree(h, k=8)
  Y$celltype <- groups[Y$celltype]
  W <- ddply(Y, .(total, celltype),
             function(A) {data.frame(y=quantile(A$y, .5),
                                     ymin=quantile(A$y, .25),
                                     ymax=quantile(A$y, .75))})
  (rrplot(W, 0, 100000, labs(x='SNPs by p-value', y='Cumulative deviation'), direct.labels) +
   geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=factor(celltype)), color=NA, alpha=.15) +
   scale_color_brewer(palette='Dark2') +
   scale_fill_brewer(palette='Dark2'))
}

rrplot.device <- function(X, aspect, type) {
  panelsize <- 60 - margin
  h <- panelsize * length(table(X$phenotype))
  w <- (aspect * panelsize + margin) * length(table(X$feature))
  Cairo(file=sub('.in.gz$', sprintf('.%s', type), args[1]), type=type, width=w, height=h, units='mm', dpi='auto', family='Helvetica')
}

rrplot.draw <- function(bin.file, assoc.file, aspect, type, ticks.fn, rrplot.fn, ...) {
  X <- read.csv(gzfile(bin.file), header=FALSE)
  Y <- read.table(gzfile(assoc.file), header=FALSE)
  colnames(X) <- c('total', 'phenotype', 'celltype', 'feature', 'count', 'expected')
  rrplot.device(X, as.numeric(aspect), type)
  P <- do.call(rrplot.fn, list(X, ...)) + do.call(ticks.fn, list(Y))
  t <- ggplot_gtable(ggplot_build(P))
  t$layout$clip[t$layout$name == 'panel'] <- 'off'
  grid.draw(t)
  print(P)
  dev.off()
}
