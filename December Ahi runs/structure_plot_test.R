structure.plot.test <- function (q.mat, pop.col = 3, prob.col = 4, sort.probs = F, 
          label.pops = TRUE, col = c("red", "green4", "blue", "magenta", "orange", "grey10"), cex = 0.5, ...) 
{
  prob.cols <- prob.col:ncol(q.mat)
  if (sort.probs) {
    q.mat <- do.call(rbind, by(q.mat, list(q.mat[, pop.col]), 
                               function(x) {
                                 prob.median <- sapply(prob.cols, function(i) median(x[, 
                                                                                       i]))
                                 pop.order <- prob.cols[order(prob.median, decreasing = TRUE)]
                                 order.list <- lapply(pop.order, function(i) x[, 
                                                                               i])
                                 order.list$decreasing = TRUE
                                 x[do.call(order, order.list), ]
                               }, simplify = FALSE))
  }
  q.mat <- q.mat[order(q.mat[, pop.col]), ]
  q.mat[, prob.cols] <- prop.table(as.matrix(q.mat[, prob.cols]), 
                                   1)
  assign.mat <- t(q.mat[, prob.cols])
  if (is.null(col)) 
    col <- rainbow(length(prob.cols))
  bp <- barplot(assign.mat, axisnames = FALSE, col = col, ...)
  tx <- tapply(bp, q.mat[, pop.col], min) + 0.1
  tx <- c(tx - tx[1], max(bp) + tx[1])
  axis(1, at = tx, labels = FALSE)
  if (label.pops) {
    lbl.x <- sapply(1:(length(tx) - 1), function(i) tx[i] + 
                      (tx[i + 1] - tx[i])/2)
    names(lbl.x) <- names(tx)[1:(length(tx) - 1)]
    mtext(names(lbl.x), side = 1, at = lbl.x, cex = cex)
  }
  x <- invisible(list(q.mat = q.mat, bar.centers = bp, pop.ticks = tx))
  abline(v=x$pop.ticks, col = "white", lwd =1.5)
}