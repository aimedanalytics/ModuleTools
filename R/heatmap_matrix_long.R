#' Prepare ordered heatmap matrix and long-format data
#'
#' Takes a data frame containing a `gene_label` column and numeric expression
#' columns (samples), builds a numeric matrix, optionally clusters rows/columns,
#' reorders the matrix and labels, and returns a long-format data.frame suitable
#' for JSON serialization alongside axis labels.
#'
#' @param df data.frame with a column `gene_label` and additional numeric columns
#'   representing samples.
#' @param cluster_rows logical. If TRUE and there are at least 2 rows, perform
#'   hierarchical clustering on rows (Euclidean distance, complete linkage) to
#'   define row order.
#' @param cluster_cols logical. If TRUE and there are at least 2 columns, perform
#'   hierarchical clustering on columns (Euclidean distance, complete linkage) to
#'   define column order.
#'
#' @return A list with:
#'   - matrix_long: data.frame with columns `row`, `col`, `value`
#'   - x_labels: character vector of column (sample) labels in display order
#'   - y_labels: character vector of row (gene) labels in display order
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @examples
#' # Example assuming df has gene_label + sample columns
#' # df <- data.frame(gene_label = paste0("G", 1:5), S1 = rnorm(5), S2 = rnorm(5))
#' # out <- heatmap_matrix_long(df, cluster_rows = TRUE, cluster_cols = TRUE)
#' # str(out$matrix_long)
#'
#' @export
heatmap_matrix_long <- function(df,
                                      cluster_rows = TRUE,
                                      cluster_cols = TRUE) {
    # Basic checks
    if (!is.data.frame(df)) {
        stop("`df` must be a data.frame or tibble.")
    }
    if (!"gene_label" %in% colnames(df)) {
        stop("`df` must contain a `gene_label` column.")
    }

    # Extract numeric matrix (drop gene_label)
    num_cols <- setdiff(colnames(df), "gene_label")
    if (length(num_cols) == 0) {
        stop("No numeric sample columns found (besides `gene_label`).")
    }
    # Coerce to matrix, ensuring numeric
    mat <- as.matrix(df[, num_cols, drop = FALSE])
    if (!is.numeric(mat)) {
        stop("Sample columns must be numeric.")
    }

    # Compute clustering
    hc_rows <- if (isTRUE(cluster_rows) && nrow(mat) > 1) {
        stats::hclust(stats::dist(mat, method = "euclidean"), method = "complete")
    } else NULL

    hc_cols <- if (isTRUE(cluster_cols) && ncol(mat) > 1) {
        stats::hclust(stats::dist(t(mat), method = "euclidean"), method = "complete")
    } else NULL

    # Reorder matrix and labels
    row_order <- if (!is.null(hc_rows)) hc_rows$order else seq_len(nrow(mat))
    col_order <- if (!is.null(hc_cols)) hc_cols$order else seq_len(ncol(mat))

    mat_ord <- mat[row_order, col_order, drop = FALSE]

    x_labels <- colnames(mat)[col_order]
    y_labels <- df$gene_label[row_order]

    # Build long format
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
