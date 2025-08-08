#' Summarise BAF QC from CNV Segment Annotation
#'
#' Summarises the number and size of "PASS" segments per chromosome from a CNV annotation file.
#'
#' @param annofile Path to the CNV annotation file (e.g., *_CNV_annotation.tsv).
#' @param out_dir Path to QC output file.
#' @param prefix Prefix of QC output file.
#'
#' @return Invisibly returns the summary data frame.
#'
#' @importFrom dplyr select group_by summarise ungroup mutate filter left_join rename all_of arrange
#' @importFrom data.table fread
#' @importFrom utils write.table
#' @export
BafQC <- function(annofile, out_dir, prefix) {
  data <- data.table::fread(annofile)
  select_col <- c("chrom", "loc.start", "loc.end", "FILTER")
  data <- data %>% dplyr::select(dplyr::all_of(select_col))
  data$length <- data$loc.end - data$loc.start

  # summary overall pass count
  summary_pass_count <- data %>%
    dplyr::select(chrom, FILTER) %>%
    dplyr::group_by(chrom, FILTER) %>%
    dplyr::summarise(count = n(), .groups = "drop_last") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(PASS_Seg_Percent = count / sum(count))

  # summary overall pass size
  summary_pass_size <- data %>%
    dplyr::select(chrom, FILTER, length) %>%
    dplyr::group_by(chrom, FILTER) %>%
    dplyr::summarise(size = sum(length), .groups = "drop_last") %>%
    dplyr::ungroup() %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(PASS_Seg_Size_Percent = size / sum(size))

  summary <- summary_pass_count %>%
    dplyr::left_join(summary_pass_size, by = c("chrom", "FILTER")) %>%
    dplyr::group_by(chrom) %>%
    dplyr::mutate(
      Total_segment_count = sum(count),
      Total_segment_size = sum(size)
    ) %>%
    dplyr::filter(FILTER == "PASS") %>%
    dplyr::filter(chrom %in% as.character(1:22)) %>%
    dplyr::select(
      chrom, FILTER, Total_segment_count, count, PASS_Seg_Percent,
      Total_segment_size, size, PASS_Seg_Size_Percent
    ) %>%
    dplyr::rename(
      PASS_Seg_Count = count,
      PASS_Seg_Size = size
    )

  chrom_levels <- c(seq(1,22),"X","Y")
  summary$chrom <- factor(summary$chrom,levels = chrom_levels)
  summary <- summary %>%
    dplyr::arrange( chrom)

  outfile_chr <- paste0(out_dir,"/", prefix,"_PASS_STAT_chr.txt")

  if (nrow(summary) == 0) {
    summary <- data.frame(PASS_Status = "No segment pass the QC!")
  }
  utils::write.table(summary, file = outfile_chr, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(summary)
}
