#' Plot CNV Profile (Copy Number, BAF, Clonality, Quality)
#'
#' Plots the CNV profile using coverage, B-allele, segment, and whitelist data.
#'
#' @param seg Path to the GATK-processed segment file.
#' @param cr Path to the GATK denoised CR file with extension denoisedCR.tsv.
#' @param ballele Path to the DRAGEN B-allele file. Or other souce B-allele file with column chrom ref alt baf
#' @param ai_binsize Numeric. Bin size for AI plot (default 100000).
#' @param cov_binsize Numeric. Bin size for coverage plot (default 100000).
#' @param whitelist Path to whitelist file.
#' @param gender Character. "male" or "female".
#' @param out_dir Output directory.
#' @param prefix Sample ID or output prefix.
#' @param aitype Character. Type of allelic imbalance data: \code{"gatk"} or \code{"other"} or \code{"dragen"}
#'   If \code{"GATK"}, the input must include columns \code{CONTIG}, \code{POSITION}, \code{ALT_COUNT}, \code{REF_COUNT}, \code{REF_NUCLEOTIDE}, and \code{ALT_NUCLEOTIDE}.
#'   If \code{"dragen"}, the input must include columns \code{contig}, \code{start}, \code{stop}, \code{refAllele}, \code{allele1}, \code{allele2}, \code{allele1Count}, \code{allele2Count}, \code{allele1AF}, and \code{allele2AF}.
#'   If \code{"other"}, the input must include columns \code{CONTIG}, \code{POSITION}, \code{ALT_COUNT}, \code{REF_COUNT}, \code{REF_NUCLEOTIDE}, and \code{ALT_NUCLEOTIDE}.
#'
#' @return Invisibly returns the ggplot object.
#'
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate select
#' @importFrom ggplot2 ggsave
#' @export
RunPlotCNV <- function(
    seg, cr, ballele,
    ai_binsize = 100000, cov_binsize = 100000,
    whitelist, gender = NULL,
    out_dir, prefix, aitype
) {
  # Read files
  filtered_lines <- grep("^[^@]", readLines(cr), value = TRUE)
  cov <- read.table(text = paste(filtered_lines, collapse = "\n"), sep = "\t", header = TRUE)
  cov$size <- cov$END - cov$START
  seg <- data.table::fread(seg)
  if(gender == "female"){
    seg <- seg %>% dplyr::filter( chrom != "Y")
  }
  cov <- CorrectGender(cov = cov, seg = seg)
  if(gender == "female"){
    cov <- cov %>% dplyr::filter(CONTIG != "Y")
  }
  purity <- as.numeric(unique(seg$rho))
  purity <- purity[!is.na(purity)]
  sf <- as.numeric(unique(seg$mu))
  sf <- sf[!is.na(sf)]
  cov <- CRCorrectPurity(cr = cov, gender = gender, purity = purity, sf = sf)
  smooth_cov <- CovBin(cov = cov, cov_binsize = cov_binsize)

  # Read BAF file
  if(aitype == "dragen"){
    ai <- data.table::fread(ballele)
    ai <- ai[,c("contig","start","refAllele","allele2","allele1Count","allele2Count")]
    colnames(ai) <- c("CONTIG", "POSITION", "REF_NUCLEOTIDE", "ALT_NUCLEOTIDE", "REF_COUNT", "ALT_COUNT")
  }else if(aitype == "gatk" || aitype == "other"){
    filtered_lines <- grep("^[^@]", readLines(ballele), value = TRUE)
    ai <- read.table(text = paste(filtered_lines, collapse = "\n"), sep = "\t", header = TRUE)
    ai <- ai %>% dplyr::filter(REF_COUNT + ALT_COUNT >= 10) %>%
      dplyr::select(CONTIG, POSITION, REF_NUCLEOTIDE, ALT_NUCLEOTIDE, REF_COUNT, ALT_COUNT) %>%
      dplyr::filter(!grepl(",", REF_NUCLEOTIDE) & !grepl(",", ALT_NUCLEOTIDE))
  }

  allelePreference <- setNames(0:3, c("A", "T", "G", "C"))
  ai <- ai %>%
    dplyr::mutate(af = ifelse(allelePreference[REF_NUCLEOTIDE] > allelePreference[ALT_NUCLEOTIDE],
                              REF_COUNT / (REF_COUNT + ALT_COUNT),
                              ALT_COUNT / (REF_COUNT + ALT_COUNT)))
  ai <- split(ai, f= ai$CONTIG)
  rounded_ai <- RoundAI(ai = ai)

  # Only keep whitelist region
  whitelist <- utils::read.table(whitelist ,header = T, stringsAsFactors = F,quote = "")
  whitelist <- whitelist %>%
    dplyr::filter(chrom %in% c(1:22,"X","Y"))

  if("Y" %in% whitelist$chrom & gender == 'female'){
    stop("Chromosome Y is present in the female whitelist, please double check...")
  } else if( ! "Y" %in% unique(whitelist$chrom) & gender == "male"){
    stop("Chromosome Y is not detected in the male whitelist, please double check ......")
  }
  clean_cov <- KeepWhitelistCov(smooth_cov = smooth_cov, whitelist = whitelist, gender = gender)
  clean_ai <- CleanHomalt(rounded_ai = rounded_ai, call_seg = seg)
  final_ai <- SmoothAI(df = clean_ai, ai_binsize = ai_binsize, gender = gender)

  colnames(seg)[1] <- "seqnames"
  final_plot <- PlotCov(
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
