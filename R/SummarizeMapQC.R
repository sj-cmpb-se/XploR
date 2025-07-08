#' Summarise QC Tables for a Sample
#'
#' Reads and summarizes mapping and coverage QC metrics for a sample, and writes a summary table.
#'
#' @param input_dir Path to the directory containing the QC files.
#' @param prefix Sample ID and file prefix.
#' @param output_dir Output directory.
#'
#' @return Invisibly returns the QC summary data frame.
#'
#' @importFrom dplyr filter select
#' @importFrom data.table fread
#' @importFrom utils read.csv write.table
#' @export
SummarizeMapQC <- function(input_dir, prefix, output_dir) {
  header <- c("info", "value", "pct")
  map_file <- file.path(input_dir, paste0(prefix, ".mapping_metrics.csv"))
  map <- utils::read.csv(map_file, header = FALSE)
  map <- map %>% dplyr::filter(V2 == "") %>% dplyr::select(3:5)
  colnames(map) <- header

  cov_file <- file.path(input_dir, paste0(prefix, ".target_bed_coverage_metrics.csv"))
  cov <- utils::read.csv(cov_file, header = FALSE)
  cov <- cov[, 3:5]
  colnames(cov) <- header

  cov2_file <- file.path(input_dir, paste0(prefix, ".target_bed_fine_hist.csv"))
  cov2 <- utils::read.csv(cov2_file, header = TRUE)
  bins <- c(10, 20, 30, 40, 45, 50, 60, 70, 80, 90, 100, 150, 200)
  bin_cov <- function(Coverage, windows) {
    Coverage$Overall <- as.numeric(Coverage$Overall)
    Coverage$pct <- Coverage$Overall / sum(Coverage$Overall)
    pct <- sapply(windows, function(x) {
      c_pct <- round(1 - sum(Coverage[1:x, "pct"]), digits = 2)
      return(c_pct)
    })
    return(pct)
  }
  cov_hist <- data.frame(
    info = paste0("pct over cov X", bins),
    pct = bin_cov(Coverage = cov2, windows = bins)
  )
  cov_hist$value <- NA

  map <- map[c(1,2,4,5,7,8,15,24,26,27,31:36,46,56,63,70:75), ]
  cov <- cov[c(2,3,4,5,34), ]
  qc <- rbind(map, cov, cov_hist)
  qc$sample <- prefix

  outFile <- file.path(output_dir, paste0(prefix, "_QC_summary.tsv"))
  utils::write.table(qc, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(qc)
}
