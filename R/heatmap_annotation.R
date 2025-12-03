#' Build heatmap annotations for ECharts JSON
#'
#' Constructs annotation metadata for heatmap visualizations based on sample-level
#' metadata columns. For each requested annotation field, this function:
#'   - Converts values to character labels, replacing NA with "(missing)"
#'   - Computes a stable code per unique label (1..K in label appearance order)
#'   - Produces a legend mapping codes to labels
#'   - Produces a row-wise mapping aligned to the provided sample column
#'
#' The returned structure matches the schema expected by the heatmap JSON:
#' a list of annotation blocks, each with:
#'   - name: the original column name
#'   - field: the generated code field name (e.g., "<col>Code")
#'   - legend: list of {code, label}
#'   - rows: list of per-sample entries with {sample, <field>=code}
#'
#' @param meta_data data.frame or tibble with sample-level metadata. Must contain
#'   the sample identifier column given by `sample_column`, and all columns listed
#'   in `col.annotation`.
#' @param sample_column character(1). Name of the column in `meta_data` that
#'   contains sample IDs. Ordering of rows in the output follows the order of
#'   `meta_data[[sample_column]]`.
#' @param col.annotation character vector. Names of columns in `meta_data` to
#'   encode as annotations.
#'
#' @return A list of annotation blocks. Each block is a list with elements:
#'   - name (character)
#'   - field (character)
#'   - legend (list of lists with elements `code` (integer) and `label` (character))
#'   - rows (list of lists; each has `sample` and the code field)
#'
#' @details
#' Label coding preserves the first-seen order of unique labels in the column.
#' Missing values are represented by the label "(missing)" and receive a code
#' like any other label.
#'
#' @examples
#' meta <- data.frame(
#'   sample_id = c("S1","S2","S3","S4"),
#'   group = c("A","A","B","B"),
#'   batch = c("X","Y",NA,"Y"),
#'   stringsAsFactors = FALSE
#' )
#' ann <- heatmap_annotation(
#'   meta_data = meta,
#'   sample_column = "sample_id",
#'   col.annotation = c("group","batch")
#' )
#' str(ann, max.level = 2)
#'
#' @export
heatmap_annotation <- function(meta_data,
                               sample_column,
                               col.annotation) {
    # Basic validations
    if (!is.data.frame(meta_data)) {
        stop("`meta_data` must be a data.frame or tibble.")
    }
    if (!is.character(sample_column) || length(sample_column) != 1) {
        stop("`sample_column` must be a single character string.")
    }
    if (!sample_column %in% colnames(meta_data)) {
        stop(sprintf("`sample_column` '%s' not found in `meta_data`.", sample_column))
    }
    if (!is.character(col.annotation) || length(col.annotation) == 0) {
        stop("`col.annotation` must be a non-empty character vector of column names.")
    }
    missing_cols <- setdiff(col.annotation, colnames(meta_data))
    if (length(missing_cols) > 0) {
        stop(sprintf("Columns not found in `meta_data`: %s",
                     paste(missing_cols, collapse = ", ")))
    }

    samples <- as.character(meta_data[[sample_column]])
    annotation <- vector("list", length(col.annotation))

    for (idx in seq_along(col.annotation)) {
        col_name <- col.annotation[[idx]]
        vals <- meta_data[[col_name]]

        # Normalize labels and build codes
        labels <- ifelse(is.na(vals), "(missing)", as.character(vals))

        # Preserve first-seen order of unique labels
        uniq <- labels[!duplicated(labels)]
        codes <- seq_along(uniq)
        label_to_code <- stats::setNames(codes, uniq)

        # Build legend
        legend <- lapply(seq_along(uniq), function(k) {
            list(
                code  = unname(codes[k]),
                label = uniq[k]
            )
        })

        # Field name and rows
        field_name <- paste0(col_name, "Code")
        rows <- lapply(seq_along(samples), function(r) {
            rec <- list(sample = samples[r])
            rec[[field_name]] <- unname(label_to_code[[labels[r]]])
            rec
        })

        annotation[[idx]] <- list(
            name   = col_name,
            field  = field_name,
            legend = legend,
            rows   = rows
        )
    }

    annotation
}

