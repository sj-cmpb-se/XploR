#' Plot CNV Profile (Copy Number, BAF, Clonality, Quality)
#'
#' Plots the CNV profile using coverage, B-allele, segment, and whitelist data.
#'
#' @param seg Path to the GATK-processed segment file.
#' @param coverage Path to the coverage file.
#' @param ballele Path to the DRAGEN B-allele file.
#' @param ai_binsize Numeric. Bin size for AI plot (default 100000).
#' @param cov_binsize Numeric. Bin size for coverage plot (default 100000).
#' @param whitelist Path to whitelist file.
#' @param gender Character. "male" or "female".
#' @param out_dir Output directory.
#' @param prefix Sample ID or output prefix.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate select
#' @importFrom ggplot2 ggsave
#' @export
runPlotCNV <- function(
    seg, coverage, ballele,
    ai_binsize = 100000, cov_binsize = 100000,
    whitelist, gender = "female",
    out_dir, prefix
) {
  # Read files
  filtered_lines <- grep("^[^@]", readLines(coverage), value = TRUE)
  cov <- read.table(text = paste(filtered_lines, collapse = "\n"), sep = "\t", header = TRUE)
  cov$size <- cov$END - cov$START
  seg <- data.table::fread(seg)
  cov <- CorrectGender(cov = cov, seg = seg)
  purity <- as.numeric(unique(seg$rho))
  purity <- purity[!is.na(purity)]
  sf <- as.numeric(unique(seg$mu))
  sf <- sf[!is.na(sf)]
  cov <- correct_purity(cov = cov, gender = gender, purity = purity, sf = sf)
  smooth_cov <- cov_bin(cov = cov, cov_binsize = cov_binsize)

  # Read BAF file
  ai <- data.table::fread(ballele)
  ai <- ai %>% dplyr::filter(allele1Count + allele2Count >= 10) %>%
    dplyr::select(contig, start, stop, refAllele, allele1, allele2, allele1Count, allele2Count) %>%
    dplyr::filter(!grepl(",", allele1) & !grepl(",", allele2))
  allelePreference <- setNames(0:3, c("A", "T", "G", "C"))
  ai <- ai %>%
    dplyr::mutate(af = ifelse(allelePreference[allele1] > allelePreference[allele2],
                              allele1Count / (allele1Count + allele2Count),
                              allele2Count / (allele1Count + allele2Count)))
  ai <- split(ai, by = "contig")
  rounded_ai <- round_ai(ai = ai)

  # Only keep whitelist region
  whitelist <- read_whitelist(whitelist, gender)
  clean_cov <- keep_whitelist_cov(smooth_cov = smooth_cov, whitelist = whitelist)
  clean_ai <- clean_homalt(rounded_ai = rounded_ai, call_seg = seg)
  final_ai <- smooth_ai(df = clean_ai, ai_binsize = ai_binsize, gender = gender)

  colnames(seg)[1] <- "seqnames"
  final_plot <- plot_cov(
    df_cov = clean_cov,
    df_ai = final_ai,
    call_seg = seg,
    whitelist = whitelist,
    gender = gender,
    prefix = prefix,
    purity = purity
  )

  outFile <- paste0(prefix, "_CNV_DRAGEN_plot.png")
  ggplot2::ggsave(
    filename = outFile,
    plot = final_plot,
    device = "png",
    path = out_dir,
    width = 15, height = 8, dpi = 600, bg = 'white'
  )
  message(paste0("Plotting is DONE. The plot is saved at: ", out_dir, "/", outFile))

  invisible(final_plot)
}
