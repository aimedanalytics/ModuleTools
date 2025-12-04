#' Build dendrogram polylines with internal clustering
#'
#' Computes hierarchical clustering for rows and/or columns from a numeric matrix
#' and converts the resulting dendrograms to NA-separated polyline coordinates
#' suitable for ECharts. The output is exactly:
#' `list(topLines = <cols dendro>, rightLines = <rows dendro>)`.
#'
#' Row dendrogram is rotated so that x = height and y = row positions.
#'
#' @param mat numeric matrix. Rows are genes/features, columns are samples.
#' @param cluster_rows logical. If TRUE and nrow(mat) > 1, compute row clustering
#'   using Euclidean distance and complete linkage. If FALSE, omit row dendrogram.
#' @param cluster_cols logical. If TRUE and ncol(mat) > 1, compute column clustering
#'   using Euclidean distance and complete linkage. If FALSE, omit column dendrogram.
#'
#' @return A list with two elements:
#'   - topLines: list with numeric vectors `x` and `y`, or NA if not computed
#'   - rightLines: list with numeric vectors `x` and `y`, or NA if not computed
#'
#' @importFrom ggdendro dendro_data
#' @examples
#' set.seed(1)
#' m <- matrix(rnorm(100), nrow = 10)
#' heatmap_dendrograms(m, cluster_rows = TRUE, cluster_cols = TRUE)
#'
#' @export
heatmap_dendrograms <- function(mat,
                                         cluster_rows = TRUE,
                                         cluster_cols = TRUE) {
    if (!is.matrix(mat) || !is.numeric(mat)) {
        stop("`mat` must be a numeric matrix.")
    }

    # Compute clustering as requested
    hc_rows <- NULL
    hc_cols <- NULL
    if (isTRUE(cluster_rows) && nrow(mat) > 1) {
        hc_rows <- stats::hclust(stats::dist(mat, method = "euclidean"), method = "complete")
    }
    if (isTRUE(cluster_cols) && ncol(mat) > 1) {
        hc_cols <- stats::hclust(stats::dist(t(mat), method = "euclidean"), method = "complete")
    }

    # Helper: convert dendro segments to NA-separated polylines
    segments_to_polylines <- function(seg) {
        if (is.null(seg) || nrow(seg) == 0) {
            return(list(x = NA_real_, y = NA_real_))
        }
        x  <- as.numeric(seg$x)
        y  <- as.numeric(seg$y)
        xe <- as.numeric(seg$xend)
        ye <- as.numeric(seg$yend)

        n <- nrow(seg)
        X <- numeric(0)
        Y <- numeric(0)
        for (i in seq_len(n)) {
            X <- c(X, x[i],  xe[i],  NA_real_)
            Y <- c(Y, y[i],  ye[i],  NA_real_)
        }
        list(x = X, y = Y)
    }

    # Columns/top dendrogram
    top_lines <- NA
    if (!is.null(hc_cols) && length(hc_cols$order) > 1) {
        dd_cols <- ggdendro::dendro_data(hc_cols)
        top_lines <- segments_to_polylines(dd_cols$segments)
    }

    # Rows/right dendrogram (rotated: x = height, y = row positions)
    right_lines <- NA
    if (!is.null(hc_rows) && length(hc_rows$order) > 1) {
        dd_rows <- ggdendro::dendro_data(hc_rows)
        seg <- dd_rows$segments

        seg_rot <- seg
        seg_rot$x    <- as.numeric(seg$y)
        seg_rot$y    <- as.numeric(seg$x)
        seg_rot$xend <- as.numeric(seg$yend)
        seg_rot$yend <- as.numeric(seg$xend)

        right_lines <- segments_to_polylines(seg_rot)
    }

    list(
        topLines  = top_lines,
        rightLines = right_lines
    )
}
