#' @title Plot Motif Match Positions
#'
#' @description
#' Plot the distribution of motif matches across sequences
#'
#' @param matches Output from \link{getPwmMatches}
#' @param binwidth Width of bins to use when plotting
#' @param abs logical(1) Plot absolute distances from centre
#' @param type Plot match density, the CDF or a binned heatmap
#' @param ... Passed to individual geom* functions
#'
#'
#' @importFrom S4Vectors splitAsList
#' @import ggplot2
#' @export
plotMatchPos <- function(
    matches, binwidth = 10, abs = FALSE, type = c("density", "cdf", "heatmap"),
    ...
){
  .checkInputBM(matches)
  type <- match.arg(type)
  if (!is(matches, "DataFrame")) {
    if (is.null(names(matches))) names(matches) <- seq_along(matches)
    nm <- names(matches)
    if (any(duplicated(nm)))
      stop("List elements in matches must have unique names")
    binned <- lapply(matches, .makeBmBins, binwidth = binwidth, abs = abs)
    bins <- lapply(binned, \(x){
      vapply(splitAsList(x, x$bin), \(df) sum(df$weight), numeric(1))
    })
    df_list <- lapply(nm, \(x){
      df <- data.frame(
        name = x, bin = names(bins[[x]]), total = as.numeric(bins[[x]])
      )
      df$p <- df$total / sum(df$total)
      df$cdf <- cumsum(df$p)
      df
    })
    df <- do.call("rbind", df_list)
    df$name <- factor(df$name, levels = nm)
    name <- sym("name")

  } else {
    binned <- motifTestR:::.makeBmBins(matches, binwidth, abs)
    bins <- vapply(splitAsList(binned, binned$bin), \(x) sum(x$weight), numeric(1))
    df <- data.frame(bin = names(bins), total = as.numeric(bins))
    df$p <- df$total / sum(df$total)
    df$cdf <- cumsum(df$p)
    name <- NULL
  }
  df$bin_start <- as.numeric(gsub("^.(.+),.+", "\\1", df$bin))
  df$bin_end <- as.numeric(gsub("^.+,(.+)\\]", "\\1", df$bin))
  df$bin_centre <- (df$bin_end + df$bin_start) / 2
  df$bin <- factor(df$bin, unique(df$bin))
  x <- sym("bin_centre")
  if (abs) x <- sym("bin_start")

  if (type == "density") {
    p <- ggplot(df, aes({{ x }}, p, colour = {{ name }})) +
      geom_smooth(se = FALSE, method = "loess", formula = 'y ~ x', ...)
  }
  if (type == "cdf") {
    int <- ifelse(abs, 0, 0.5)
    p <- ggplot(df, aes({{ x }}, cdf, colour = {{ name }})) +
      geom_line(...)
  }
  if (type == "heatmap") {
    p <- ggplot(df, aes(bin, name, fill = p)) +
      geom_tile(...) +
      scale_x_discrete(expand = rep_len(0, 4)) +
      scale_y_discrete(expand = rep_len(0, 4)) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  }

  p

}
