#' Run AI Segmentation Workflow
#'
#' This function executes the full AI segmentation workflow, including reading input files, merging segments, BAF processing, and writing results.
#'
#' @param seg Path to the GATK segment file.
#' @param cov Path to the GATK coverage count file.
#' @param ai Path to the DRAGEN BAF bedgraph file.
#' @param gender Character, "male" or "female".
#' @param out_dir Output directory path.
#' @param prefix Output file prefix.
#' @param aibinsize Numeric, AI bin size (default: 500000).
#' @param mergeai Numeric, BAF difference threshold for merging segments (default: 0.2).
#' @param mergecov Numeric, CNV difference threshold for merging segments (default: 0.2).
#' @param minaisize Numeric, smallest AI segment size (default: 1000000).
#' @param snpmin Numeric, minimum SNPs for BAF segmentation (default: 7).
#' @param minsnpcallaicutoff Numeric, minimum SNPs for reliable CNLOH/GAINLOH (default: 10).
#' @param mergecovminsize Numeric, minimum size for GATK segment merge (default: 500000).
#' @return Invisibly returns the output file path.
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate rowwise ungroup select
#' @importFrom stringr str_split
#' @importFrom tidyr unnest_wider
#' @export
runAIsegmentation <- function(
    seg, cov, ai, gender, out_dir, prefix,
    aibinsize = 500000,
    mergeai = 0.2,
    mergecov = 0.2,
    minaisize = 1000000,
    snpmin = 7,
    minsnpcallaicutoff = 10,
    mergecovminsize = 500000
) {
  # Read BAF file
  baf <- data.table::fread(ai)
  colnames(baf) <- c("Chromosome", "Start", "End", "baf")
  baf <- baf %>%
    dplyr::filter(!Chromosome %in% c("X", "Y")) %>%
    dplyr::filter(baf != 0 & baf != 1) %>%
    dplyr::mutate(baf = ifelse(baf > 0.5, 1 - baf, baf))

  # Read coverage file
  lines <- readLines(cov)
  filtered_lines <- grep("^\\s*@", lines, invert = TRUE, value = TRUE)
  cov_df <- data.table::fread(text = filtered_lines)

  # Read segment file
  seg_df <- data.table::fread(seg) %>%
    dplyr::mutate(size = End - Start)

  # Prepare options list for internal function compatibility
  opt <- list(
    gender = tolower(gender),
    aibinsize = aibinsize,
    mergeai = mergeai,
    mergecov = mergecov,
    minaisize = minaisize,
    snpmin = snpmin,
    minsnpcallaicutoff = minsnpcallaicutoff,
    mergecovminsize = mergecovminsize
  )

  # Step 1: Gender check
  seg_df <- CheckGender(cov = cov_df, seg = seg_df, gender = opt$gender, opt = opt)
  merge_seg <- seg_df

  # Step 2: Merge GATK segments by size
  if (opt$mergecovminsize < 1e+05) opt$mergecovminsize <- 1e+05
  for (minsize in c(1000, 10000, 50000, seq(from = 100000, to = as.numeric(opt$mergecovminsize), by = 100000))) {
    merge_seg <- merge_seg %>%
      dplyr::ungroup() %>%
      dplyr::filter(size >= minsize | Chromosome %in% c("X", "Y"))
    merge_seg <- CallMerge(data = merge_seg, AIorSeg = "Seg", opt = opt)
  }
  merge_seg <- subset(merge_seg, select = -Call)
  merge_seg <- merge_seg %>% dplyr::mutate(BreakpointSource = "GATK")

  # Step 3: Correct BAF scale
  baf_scale <- 2
  while (baf_scale > 1.005) {
    baf_scale <- EstimateBAFshift(baf)
    if (!is.na(baf_scale) & baf_scale <= 1.2) {
      baf <- CorrectBAF(baf, baf_scale)
      baf_scale <- EstimateBAFshift(baf)
    } else {
      baf_scale <- 1
    }
  }

  # Step 4: Search for breakpoints in each segment
  combindAIseg <- lapply(seq_len(nrow(merge_seg)), function(x) {
    SearchBreakpoint(seg_row = merge_seg[x, ], baf = baf, opt = opt)
  })
  combined_result <- do.call(rbind, combindAIseg)

  # Step 5: Add quality tag
  result <- combined_result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(FILTER = AddQualTag(
      MAF_gmm_weight = MAF_gmm_weight,
      MAF_Probes = MAF_Probes,
      MAF_gmm_G = MAF_gmm_G,
      snpmin = opt$snpmin
    ))

  # Step 6: Write output
  outFile <- file.path(out_dir, paste0(prefix, "_GATK_DRAGEN_merge.tsv"))
  write.table(result, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(outFile)
}
