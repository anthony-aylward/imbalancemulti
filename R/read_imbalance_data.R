#===============================================================================
# read_imbalance_data.R
#===============================================================================

# Read imbalance data from multiple experiments




# Imports ======================================================================

#' @import parallel




# Functions ====================================================================

#' @title Read Imbalance Data
#'
#' @description Read imbalance data from multiple experiments
#'
#' @param ... paths to imbalance data files
#' @param experiment_names character, names of experiments
#' @param cores integer, number of cores to use
#' @return list containing matrices of allele counts, indexed by variant id
#' @export
read_imbalance_data <- function(..., experiment_names = NULL, cores = 1) {
  files <- unlist(list(...))
  if (!is.null(experiment_names)) names(experiment_names) <- files
  counts <- list()
  for (file in files) {
    if (is.null(experiment_names)) {
      experiment_name <- sub(".tsv", "", basename(file), fixed = TRUE)
    } else {
      experiment_name <- experiment_names[[file]]
    }
    counts_frame <- read.table(
      file,
      header = TRUE,
      row.names = "id",
      stringsAsFactors = FALSE
    )[c("coverage", "ref_count")]
    counts_list <- setNames(
      mclapply(
        rownames(counts_frame),
        function(id) {
          row <- counts_frame[id,]
          matrix(
            c(row[["ref_count"]], row[["coverage"]] - row[["ref_count"]]),
            nrow = 1,
            dimnames = list(experiment_name, c("ref", "alt"))
          )
        },
        mc.cores = cores
      ),
      rownames(counts_frame)
    )
    for (id in names(counts_list)) {
      if (id %in% names(counts)) {
        counts[[id]] <- rbind(counts[[id]], counts_list[[id]])
      } else {
        counts[[id]] <- counts_list[[id]]
      }
    }
  }
  counts
}
