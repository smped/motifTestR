#' @title Plot Motif Match Positions
#'
#' @description
#' Plot the distribution of motif matches across sequences
#'
#' @details
#' Multiple options are provided for showing the distribution of PWM matches
#' within a set of sequences, using either the smoothed probability density,
#' the probability CDF or as a heatmap.
#' Distances can be shown as symmetrical around the centre or using absolute
#' distances from the central position within the sequences.
#'
#' Heatmaps are only enabled for comparisons across multiple PWMs, with
#' optional clustering enabled. If adding a dendrogram for clustering,
#' the returned plot object will be a patchwork object.
#'
#' @return A ggplot2 object
#'
#'
#' @param matches Output from \link{getPwmMatches}
#' @param binwidth Width of bins to use when plotting
#' @param abs logical(1) Plot absolute distances from centre
#' @param use_totals logical(1). If TRUE, plots will use total counts. The
#' default (FALSE) plots probabilities.
#' @param type Plot match density, the CDF or a binned heatmap
#' @param geom Type of geom to be used for line plots. Ignored for heatmaps
#' @param cluster logical(1) Cluster motifs when drawing a heatmap. If TRUE a
#' dendrogram will be added to the LHS of the plot
#' @param w Relative width of the dendrogram on (0, 1)
#' @param heat_fill scale_fill_continuous object for heatmaps. If not provided,
#' scale_fill_viridis_c() will be added to the heatmap.
#' @param ... Passed to individual geom* functions
#'
#' @examples
#' ## Load the example PWM
#' data("ex_pwm")
#' esr1 <- ex_pwm$ESR1
#'
#' ## Load the example sequences from the peaks
#' data("ar_er_seq")
#'
#' ## Just the best match
#' bm <- getPwmMatches(esr1, ar_er_seq, best_only = TRUE)
#' plotMatchPos(bm, se = FALSE)
#'
#' ## Matches can also be shown by distance from centre
#' plotMatchPos(bm, abs = TRUE)
#'
#' ## Cumulative Probability plots are also implemented
#' plotMatchPos(bm, type = "cdf", geom = "line", colour = "red") +
#'   geom_abline(intercept = 0.5, slope = 1/ 400)
#'
#' @importFrom S4Vectors splitAsList
#' @importFrom rlang !! sym
#' @import ggplot2
#' @export
plotMatchPos <- function(
        matches, binwidth = 10, abs = FALSE, use_totals = FALSE,
        type = c("density", "cdf", "heatmap"),
        geom = c("smooth", "line", "point", "col"),
        cluster = FALSE, w = 0.1, heat_fill = NULL, ...
){
    .checkMatches(matches)
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
        binned <- .makeBmBins(matches, binwidth, abs)
        bins <- vapply(
            splitAsList(binned, binned$bin), \(x) sum(x$weight), numeric(1)
        )
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
    yval <- ifelse(use_totals, "total", "p")

    geom <- paste0("geom_", match.arg(geom))
    geom_fun <- match.fun(geom)
    if (type == "density") {
        plot_aes <- aes({{ x }}, !!sym(yval), colour = {{ name }})
        if (geom == "geom_col") names(plot_aes)[[3]] <- "fill"
        plot <- ggplot(df, plot_aes) + geom_fun(...)
    }
    if (type == "cdf") {
        plot <- ggplot(df, aes({{ x }}, !!sym("cdf"), colour = {{ name }})) +
            geom_fun(...)
    }
    if (type == "heatmap") {
        if (is.null(name)) stop("Heatmaps not inplemented for a single PWM")
        plot <- .createHeatmap(df, cluster, yval, w, heat_fill, ...)
    }

    plot

}

#' @importFrom stats hclust dist xtabs as.dendrogram as.formula
#' @importFrom patchwork plot_layout
#' @importFrom rlang !! sym
#' @keywords internal
.createHeatmap <- function(df, cluster, yval, w, fill_scale, ...) {

    levels <- unique(df$name)
    if (is.null(fill_scale)) fill_scale <- scale_fill_viridis_c()
    stopifnot(is(fill_scale, "ScaleContinuous"))
    stopifnot(w > 0 & w < 1)
    ylab <- "name"
    ## Create the dendrogram if needed
    if (cluster) {

        if (!requireNamespace('ggdendro', quietly = TRUE))
            stop("Please install 'ggdendro' to use this function.")

        fm <- as.formula(paste(yval, "~ name + bin"))
        mat <- as.matrix(xtabs(fm, df[c("name", "bin", yval)]))
        clust <- hclust(dist(mat), method = "ward.D2")
        dend <- ggdendro::dendro_data(as.dendrogram(clust))
        levels <- dend$labels$label
        ym <- max(dend$segments$xend)
        dend_aes <- aes(
            x = !!sym("yend"), y = !!sym("xend"),
            xend = !!sym("y"), yend = !!sym("x")
        )
        dend_plot <- ggplot(dend$segments) +
            geom_segment(dend_aes) +
            scale_x_reverse(expand = expansion(0)) +
            scale_y_continuous(limits = c(0, ym) + 0.5, expand = expansion(0)) +
            labs(x = "", y = c()) +
            theme_minimal() +
            theme(
                panel.grid = element_blank(),
                axis.text = element_blank(), axis.ticks = element_blank(),
                plot.margin = unit(c(5.5, 0, 5.5, 5.5), "points")
            )
        ylab <- NULL
    }

    # Always create the heatmap...
    df$name <- factor(df$name, levels = levels)
    plot <- ggplot(
        df, aes(!!sym("bin"), !!sym("name"), fill = !!sym(yval))
    ) +
        geom_tile(...) +
        scale_x_discrete(expand = rep_len(0, 4)) +
        scale_y_discrete(expand = rep_len(0, 4), name = ylab) +
        fill_scale +
        theme(
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
        )
    if (!cluster) return(plot)
    dend_plot + plot + plot_layout(widths = c(w, 1 - w))

}
