#' Assign each motif to a cluster
#'
#' Cluster related motifs for testing as a group
#'
#' @details
#' This builds on \link[universalmotif]{compare_motifs}, enabling the
#' assignment of each PWM to a cluster, and subsequent testing of motifs as
#' a cluster, rather than returning individual results
#'
#' @return Named vector with numeric values representing which cluster each
#' motif has been assigned to.
#'
#' @param motifs A list of universalmotifs or a list of PWMs
#' @param type Can be ICM or PPM
#' @param method The moetho to be used for determining similarity/distances
#' @param power Raise correlation matrices to this power before converting to a
#' distance matrix. Only applied if method is either "PCC" or "WPCC"
#' @param agglom Method to be used for agglomeration by \link[stats]{hclust}
#' @param thresh Tree heights below which motifs are formed into a cluster
#' @param plot Show tree produced by \link[stats]{hclust}. If requested the
#' value set by thresh will be shown as a horizontal line
#' @param labels,cex Passed to \link[stats]{plot.hclust}
#' @param linecol Passed to \link[graphics]{abline} as the argument `col`
#' @param ... passed to \link[universalmotif]{compare_motifs}
#'
#' @examples
#' # Load the example motifs
#' data("ex_pwm")
#'
#' # Return a vector with each motif assigned a cluster
#' # The default uses Pearson's Correlation Coefficient
#' clusterMotifs(ex_pwm)
#'
#' # Preview the settings noting that showing labels can clutter the plot
#' # with large numbers of motifs. The defaults for Euclidean distance
#' # show the threshold may need raising
#' clusterMotifs(ex_pwm, plot = TRUE, labels = NULL, method = "EUCL")
#'
#' @importFrom universalmotif compare_motifs create_motif
#' @importFrom stats as.dist hclust cutree
#' @importFrom graphics abline
#' @export
clusterMotifs <- function(
        motifs, type = c("PPM", "ICM"),
        method = c("PCC", "EUCL", "SW", "KL", "ALLR", "BHAT", "HELL", "SEUCL", "MAN", "ALLR_LL", "WEUCL", "WPCC"),
        power = 1, agglom = "complete", thresh = 0.2,
        plot = FALSE, labels = FALSE, cex = 1, linecol = "red", ...
){
    # Convert to universal motif, if a list is passed
    if (all(vapply(motifs, is, logical(1), class2 = "matrix"))) {
        motifs <- lapply(
            names(motifs), \(x) create_motif(motifs[[x]], name = x)
        )
    }
    stopifnot(all(vapply(motifs, is, logical(1), class2 = "universalmotif")))
    method <- match.arg(method)
    type <- match.arg(type)
    if (type == "ICM" & method %in% c("ALLR", "ALLR_LL"))
        stop("Cannot use ICM with ALLR or ALLR_LL")
    is_dist <- method %in% c("EUCL", "KL", "HELL", "SEUCL", "MAN", "WEUCL")
    mat <- compare_motifs(motifs, use.type = type, method = method, ...)
    ## This is really only useful for correlations
    if (power != 1 & method %in% c("PCC", "WPCC")) mat <- mat^power
    ## Make a distance/dissimilarity matrix
    mat <- abs(mat) / max(abs(mat)) # Scale to be in [0,1]
    if (!is_dist) mat <- 1 - abs(mat)
    d <- as.dist(mat)
    cl <- hclust(d, method = agglom)
    if (plot) {
        plot(cl, labels = labels, cex = cex)
        abline(a = thresh, b = 0, col = linecol)
    }
    cl <- cutree(cl, h = thresh)
    cl
}

