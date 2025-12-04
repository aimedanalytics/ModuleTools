#' Prepare ordered matrix and long-format data (first column as labels)
#'
#' Expects a data frame where the first column contains row labels and all
#' remaining columns are numeric. Builds a numeric matrix, optionally clusters
#' rows/columns, reorders the matrix and labels, and returns a long-format
#' data.frame alongside axis labels.
#'
#' @param df data.frame whose first column is the row-label column; remaining
#'   columns must be numeric.
#' @param cluster_rows logical. If TRUE and there are at least 2 rows, perform
#'   hierarchical clustering on rows (Euclidean distance, complete linkage).
#' @param cluster_cols logical. If TRUE and there are at least 2 columns, perform
#'   hierarchical clustering on columns (Euclidean distance, complete linkage).
#'
#' @return A list with:
#'   - matrix_long: data.frame with columns `row`, `col`, `value`
#'   - x_labels: character vector of column labels (from numeric columns) in order
#'   - y_labels: character vector of row labels (from first column) in order
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @examples
#' set.seed(1)
#' df <- data.frame(
#'   label = paste0("G", 1:5),
#'   S1 = rnorm(5),
#'   S2 = rnorm(5),
#'   S3 = rnorm(5),
#'   check.names = FALSE
#' )
#' out <- heatmap_matrix_long(df, cluster_rows = TRUE, cluster_cols = TRUE)
#' head(out$matrix_long)
#'
#' @export
heatmap_matrix_long <- function(df,
                                cluster_rows = TRUE,
                                cluster_cols = TRUE) {
    if (!is.data.frame(df)) {
        stop("`df` must be a data.frame or tibble.")
    }
    if (ncol(df) < 2) {
        stop("`df` must have at least 2 columns: one label column and >=1 numeric column.")
    }

    # First column as labels, rest must be numeric
    label_col <- df[[1]]
    if (any(is.na(label_col))) {
        # Allow NA labels but coerce to character for consistent indexing later
        label_col <- as.character(label_col)
    } else {
        label_col <- as.character(label_col)
    }

    num_df <- df[, -1, drop = FALSE]
    # Validate numeric columns
    if (!all(vapply(num_df, is.numeric, logical(1)))) {
        stop("All columns except the first must be numeric.")
    }

    mat <- as.matrix(num_df)

    # Compute clustering
    hc_rows <- if (isTRUE(cluster_rows) && nrow(mat) > 1) {
        stats::hclust(stats::dist(mat, method = "euclidean"), method = "complete")
    } else NULL

    hc_cols <- if (isTRUE(cluster_cols) && ncol(mat) > 1) {
        stats::hclust(stats::dist(t(mat), method = "euclidean"), method = "complete")
    } else NULL

    # Orders
    row_order <- if (!is.null(hc_rows)) hc_rows$order else seq_len(nrow(mat))
    col_order <- if (!is.null(hc_cols)) hc_cols$order else seq_len(ncol(mat))

    # Reorder
    mat_ord <- mat[row_order, col_order, drop = FALSE]

    x_labels <- colnames(mat)[col_order]
    y_labels <- label_col[row_order]

    # Long format
    matrix_long <- as.data.frame(mat_ord) %>%
        tibble::rownames_to_column("row_idx") %>%
        dplyr::mutate(row = y_labels[as.integer(row_idx)]) %>%
        dplyr::select(-row_idx) %>%
        tidyr::pivot_longer(cols = -row, names_to = "col", values_to = "value")

    list(
        matrix_long = matrix_long,
        x_labels = x_labels,
        y_labels = y_labels
    )
}


