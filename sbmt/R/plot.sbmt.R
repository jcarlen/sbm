#' generic function to plot block-to-block edges (i.e. sum of edges of all nodes in the block) over time (reverse direction is show in )
#' 
#' @param x output of sbmt of class sbmt
#' @param mai sets plot margin width
#' @param show_blocks is list of length two that sets which block-to-block edges to show (defaults to all)
#' @param show_reverse if true will show the reverse direction as a dashed gray line (only useful for directed, but may not want even then)
#' @param ... ignored
#' @rdname plot.sbmt
#' @export
plot.sbmt <- function(x, mai = rep(.5,4), show_blocks = list(1:nrow(x$EdgeMatrix[[1]]), 1:nrow(x$EdgeMatrix[[1]])), show_reverse = TRUE, ...) {
  plot_data_out = do.call("cbind", x$EdgeMatrix)
  plot_data_in = do.call("rbind", x$EdgeMatrix)
  G = nrow(x$EdgeMatrix[[1]])
  n_slices = length(x$EdgeMatrix)
  par(mfrow = c(length(show_blocks[[1]]), length(show_blocks[[2]])))
  par(mai = mai)
  for (i in show_blocks[[1]]) {
    for (j in show_blocks[[2]]) {
        plot(plot_data_out[i,seq(j, G*n_slices, by = G)],  type = "l", xlab = "", ylab = "", main = paste0(i,"->",j))
        if (i != j & show_reverse) {
          points(plot_data_out[j, seq(i,G*n_slices, by = G)], type = "l", lty = 2, col = "grey")
        }
    }
  }
}