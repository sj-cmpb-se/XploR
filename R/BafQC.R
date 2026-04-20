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
  model_source <- as.character(data$Model_source[1])
  purity <- data$rho[1]
  scale_factor <- data$mu[1]
  gatk_gender <- data$gatk_gender[1]
  pipeline_gender <- data$pipeline_gender[1]
  select_col <- c("chrom", "loc.start", "loc.end","cov_mad", "FILTER")
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
    dplyr::filter(chrom %in% as.character(1:22)) %>%
    dplyr::select(
      chrom, FILTER, Total_segment_count, count, PASS_Seg_Percent,
      Total_segment_size, size, PASS_Seg_Size_Percent
    ) %>%
    dplyr::rename(
      PASS_Seg_Count = count,
      PASS_Seg_Size = size
    )

  summary_pass <- summary %>%
    filter(FILTER == "PASS")

  ## chromosomes without any passed segments
  chrom_no_pass <- setdiff(c(1:22), unique(summary_pass$chrom))
  summary_fail <- summary %>%
    filter(FILTER == "FAILED" & chrom %in% chrom_no_pass) %>%
    mutate( FILTER = "PASS",
            PASS_Seg_Count = 0,
            PASS_Seg_Percent = 0,
            PASS_Seg_Size = 0,
            PASS_Seg_Size_Percent = 0 )

  summary <- rbind(summary_pass , summary_fail)

  chrom_levels <- c(seq(1,22),"X","Y")
  summary$chrom <- factor(summary$chrom,levels = chrom_levels)
  summary <- summary %>%
    dplyr::arrange( chrom)

  sample_qc <- data.frame(
    Median_chr_seg_num = median(summary$Total_segment_count,na.rm=T),
    Median_chr_cov_MAD = median(data$cov_mad, na.rm =T),
    Fail_chr_count = nrow( summary %>% filter(PASS_Seg_Size_Percent < 0.9) ),
    Model_Source = model_source,
    Purity = purity,
    DiploidCov_scale = scale_factor,
    GATK_gender = gatk_gender,
    Final_gender = pipeline_gender
  )

  sample_qc <- sample_qc %>%
    dplyr::mutate( QC_suggestion = ifelse( Median_chr_seg_num > 9 || Median_chr_cov_MAD > 2 || Fail_chr_count > 2, "Fail","PASS"))

  outfile_chr <- paste0(out_dir,"/", prefix,"_PASS_STAT_chr.txt")
  outfile_sample <- paste0(out_dir,"/", prefix,"_sample_qc.txt")

  if (nrow(summary) == 0) {
    summary <- data.frame(PASS_Status = "No segment pass the QC!")
  }
  utils::write.table(summary, file = outfile_chr, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(sample_qc, file = outfile_sample, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(summary)
}




