#' Only show the top labels
#' @param n how many labels to show
top <- function(n) {
  function(d, ...) {
    head(d[order(d$y, decreasing=TRUE),], n=n)
  }
}

#' Only show a set of labels
#' @param show vector of labels to include
filter <- function(show) {
  function(d, ...) {
    d[d$celltype %in% show,]
  }
}

#' Sequentially bump labels down if they overlap, starting from the top
#' @param d bounding boxes of labels
bumpdown <- function(d, ...) {
  if (nrow(d) > 1) {
    d <- calc.boxes(d)[order(d$y, decreasing=TRUE),]
    d$h <- d$h * 1.1
    d <- calc.borders(d)
    '%between%' <- function(v,lims)lims[1]<v&v<lims[2]
    obox <- function(x,y){
      tocheck <- with(x,c(left,(right-left)/2+left,right))
      tocheck %between% with(y,c(left,right))
    }
    for(i in 2:nrow(d)){
      dif <- d$top[i]-d$bottom[i-1]
      overlap <- c(obox(d[i,],d[i-1,]),obox(d[i-1,],d[i,]))
      if(dif > 0 && any(overlap)){
        d$bottom[i] <- d$bottom[i] - dif
        d$top[i] <- d$top[i] - dif
        d$y[i] <- d$y[i] - dif
      }
    }
  }
  d
}

