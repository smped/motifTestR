#' @importFrom methods slot is
#' @importClassesFrom universalmotif universalmotif
#' @importFrom universalmotif create_motif convert_type
#' @keywords internal
.checkPWM <- function(motif){
    ## Convert any PPM/PCM matrices to PWM
    ## May break a few things downstream
    if (is.matrix(motif)) {
        nsites <- 100
        if (is.integer(motif)) nsites <- max(matrixStats::colSums2(motif))
        motif <- create_motif(
            motif, nsites = nsites, pseudocount = 1, type = "PPM"
        )
    }
    stopifnot(is(motif, "universalmotif"))
    if (slot(motif, "pseudocount") == 0)
        warning("Zero pseudocounts may lead to a PWM with infinite values")
    pwm <- convert_type(motif, "PWM")
    slot(pwm, "motif")
}

#' @keywords internal
.viewMapFromXStringset <- function(stringset){
    stopifnot(is(stringset, "XStringSet"))
    # Form a map from the original sequences to a views object
    map <- data.frame(end = cumsum(width(stringset)), width = width(stringset))
    map$start <- map$end - map$width + 1
    map$names <- names(stringset)
    map
}

#' @keywords internal
#' @importClassesFrom universalmotif universalmotif
.cleanMotifList <- function(x) {
    if (all(vapply(x, is, logical(1), "universalmotif"))) {
        ## Add the name as the motif name, swicthing to altname where required
        n <- length(x)
        nm <- vapply(x, slot, character(1), "name")
        if (length(unique(nm)) < n) {
            message("Found non-unique values in the name slot. Trying altname")
            nm <- unlist(lapply(x, slot, "altname"))
            if (length(unique(nm)) < n) {
                msg <- "Non-unique values in the altname slot. Please resolve"
                stop(msg)
            }
        }
        x <- lapply(x, slot, "motif")
        names(x) <- nm

        # nm <- unlist(lapply(x, slot, "altname")) # Less duplication than name
        # if (all(is.na(nm))) { # Handle the case where altname is empty (NA)
        #     nm <- vapply(x, slot, character(1), "name")
        # }
        # names(x) <- nm
    }
    ## Now check everything is a matrix
    all_mat <- vapply(
        x,
        \(mat) all(is.matrix(mat), rownames(mat) == c("A", "C", "G", "T")),
        logical(1)
    )
    stopifnot(all_mat)
    x
}

#' @keywords internal
.checkMatches <- function(matches) {
    colTypes <- c(
        score = "numeric", direction = "factor",
        start = "integer", end = "integer",
        from_centre = "numeric", seq_width = "integer", match = "DNAStringSet"
    )
    reqdCols <- c("seq", names(colTypes)) # seq can be any type
    msg <- paste(reqdCols, collapse = ",")

    if (is(matches, "DataFrame")) {
        stopifnot(all(reqdCols %in% colnames(matches)))
        types <- vapply(matches, \(x) is(x)[[1]], character(1))
        stopifnot(
            identical(types[names(colTypes)], colTypes)
        )
    } else {
        stopifnot(is(matches, "list") & !is(matches, "data.frame"))
        hasCols <- vapply(
            matches, \(x) all(reqdCols %in% colnames(x)), logical(1)
        )
        if (any(!hasCols))
            stop("All objects must contain the columns: \n", msg)
        correctTypes <- vapply(
            matches,
            \(x) {
                types <- vapply(x, \(cols) is(cols)[[1]], character(1))
                identical(types[names(colTypes)], colTypes)
            }, logical(1)
        )
        if (any(!correctTypes)) stop("All columns are not of the correct type")
    }

    invisible(TRUE)
}

#' @keywords internal
.makeBmBins <- function(matches, binwidth, abs){

    ## Set sequence weights for multiple matches
    reps <- table(as.character(matches$seq)) # Counts reps
    matches$weight <- as.numeric(1 / reps[as.character(matches$seq)]) # 1 / reps

    ## Make the bins
    longest_seq <- max(matches$seq_width)
    stopifnot(binwidth < longest_seq)
    if (abs) {
        ## Set the dist from centre to be +ve only
        matches$from_centre <- abs(matches$from_centre)
        all_breaks <- c(seq(0, longest_seq / 2, by = binwidth), longest_seq / 2)
        seq_per_bin <- vapply(
            all_breaks, \(x) sum(matches$seq_width >= x), integer(1)
        )
    } else {
        ## To define bins, make sure the central bin is around zero and that
        ## bins extends symetrically from this point to the widest sequence
        pos_breaks <- c(
            ## Ensure the upper limit is included, even if it leaves a smaller bin
            seq(binwidth / 2, longest_seq / 2, by = binwidth), longest_seq / 2
        )
        all_breaks <- c(-1 * pos_breaks, pos_breaks)
        pos_bins <- vapply(
            pos_breaks, \(x) sum(matches$seq_width >= x), integer(1)
        )
        seq_per_bin <- c(pos_bins, pos_bins)
    }

    ## Cut into bins
    all_breaks <- unique(sort(all_breaks))
    matches$bin <- cut(matches$from_centre, breaks = all_breaks, include.lowest = TRUE)
    matches
}

#' @importFrom IRanges Views
#' @import Biostrings
#' @keywords internal
.hasPwmMatch <- function(pwm, stringset, rc = TRUE, min_score = "80%", ...) {
    ## Returns a logical vector the same length as the input stringset
    ## Checks & the map
    pwm <- .checkPWM(pwm)
    stopifnot(is(stringset, "XStringSet"))

    ## Handle empty stringsets
    n <- length(stringset)
    seq_with_hits <- NULL

    if (n > 0) {
        map <- .viewMapFromXStringset(stringset)
        # Form the entire XStringSetList into a Views object
        views <- Views(
            unlist(stringset), start = map$start, width = map$width,
            names = map$names
        )
        hits <- matchPWM(pwm, views, min.score = min_score, ...)
        if (rc) {
            rev_pwm <- reverseComplement(pwm)
            hits_rev <- matchPWM(rev_pwm, views, min.score = min_score, ...)
            hits <- c(hits, hits_rev)
        }
        seq_with_hits <- findInterval(start(hits), map$start)
    }
    seq_len(n) %in% seq_with_hits

}
