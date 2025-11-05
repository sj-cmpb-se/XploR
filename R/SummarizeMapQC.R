#' Summarize QC Metrics for a Sample from a DRAGEN run.
#'
#' Reads mapping and coverage quality control (QC) files for a sample, summarizes key metrics, and writes a summary table to the output directory.
#'
#' @param input_dir Character. Path to the directory containing QC files for the sample.
#' @param prefix Character. Sample ID and file prefix.
#' @param platform Character. Sequencing platform, either \code{"wgs"} (whole genome sequencing) or \code{"targeted"}.
#' @param out_dir Character. Output directory for the summary table.
#'
#' @return Invisibly returns a data frame containing the QC summary metrics.
#'
#' @details
#' The function reads mapping and coverage QC files from \code{input_dir}, extracts key metrics, and writes a summary table to \code{out_dir} with the specified \code{prefix}. The summary includes metrics relevant to the specified \code{platform}.
#'
#' @importFrom dplyr filter select
#' @importFrom data.table fread
#' @importFrom utils read.csv write.table
#' @export
SummarizeMapQC <- function(input_dir, prefix, platform, out_dir) {
  header <- c("info", "value", "pct")
  map_file <- file.path(input_dir, paste0(prefix, ".mapping_metrics.csv"))
  map <- utils::read.csv(map_file, header = FALSE)
  map <- map %>% dplyr::filter(V2 == "") %>% dplyr::select(3:5)
  colnames(map) <- header

  if( platform == "targeted"){
    cov_file <- file.path(input_dir, paste0(prefix, ".target_bed_coverage_metrics.csv"))
    cov2_file <- file.path(input_dir, paste0(prefix, ".target_bed_fine_hist.csv"))
  }else if( platform == "wgs"){
    cov_file <- file.path(input_dir, paste0(prefix, ".wgs_coverage_metrics.csv"))
    cov2_file <- file.path(input_dir, paste0(prefix, ".wgs_fine_hist.csv"))
  }

  cov <- utils::read.csv(cov_file, header = FALSE)
  cov <- cov[, 3:5]
  colnames(cov) <- header


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
  map_extract_id <- c("Total input reads","duplicate marked reads","Number of unique reads",
                      "QC-failed reads","Mapped reads","rRNA","Unmapped reads","Singleton","Properly paired reads",
                      "discordant","multiple locations","Reads with MAPQ","Total alignments","Secondary alignments",
                      "chimeric","read length","Insert length","contamination")
  map_extract_id <- paste0(map_extract_id,collapse = "|")

  #map <- map[c(1,2,4,5,7,8,15,24,26,27,31:36,46,56,63,70:75), ]
  map <- map[grep(map_extract_id, map$info),]
  cov_extract_id <- paste0(c("Average alignment coverage over target region","Uniformity",
                           "10x: inf","20x: inf","50x: inf","100x: inf","Aligned reads in target region"),collapse = "|")
  #cov <- cov[c(2,3,4,5,34), ]
  cov <- cov[grep(cov_extract_id, cov$info),]
  cov[grep("PCT|pct",cov$info),"pct"] <-  cov[grep("PCT|pct",cov$info),"value"]
  qc <- rbind(map, cov, cov_hist)
  qc$sample <- prefix

  outFile <- file.path(out_dir, paste0(prefix, "_QC_summary.tsv"))
  utils::write.table(qc, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(qc)
}
