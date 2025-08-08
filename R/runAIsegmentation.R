#' Run AI Segmentation Workflow
#'
#' This function executes the full AI segmentation workflow, including reading input files, merging segments, BAF processing, and writing results.
#'
#' @param seg Path to the GATK segment file.
#' @param cov Path to the GATK denoised coverage count file.
#' @param ai Path to the BAF file or allelic count file.
#' @param gender Character, "male" or "female".
#' @param out_dir Output directory path.
#' @param prefix Output file prefix.
#' @param aibinsize Numeric, AI bin size (default: 500000).
#' @param mergeai Numeric, MAF difference threshold for merging segments (default: 0.15).
#' @param mergecov Numeric, CNV difference threshold for merging segments (default: 0.2).
#' @param minaisize Numeric, smallest AI segment size (default: 1000000).
#' @param snpmin Numeric, minimum SNPs for MAF segmentation (default: 7).
#' @param minsnpcallaicutoff Numeric, minimum SNPs for reliable CNLOH/GAINLOH (default: 10).
#' @param mergecovminsize Numeric, minimum size for GATK segment merge (default: 500000).
#' @param aitype Character. Type of allelic imbalance data: \code{"gatk"} or \code{"other"} or \code{"dragen"}
#'   If \code{"GATK"}, the input must include columns \code{CONTIG}, \code{POSITION}, \code{ALT_COUNT}, \code{REF_COUNT}, \code{REF_NUCLEOTIDE}, and \code{ALT_NUCLEOTIDE}.
#'   If \code{"dragen"}, the input must include columns \code{contig}, \code{start}, \code{stop}, \code{refAllele}, \code{allele1}, \code{allele2}, \code{allele1Count}, \code{allele2Count}, \code{allele1AF}, and \code{allele2AF}.
#'   If \code{"other"}, the input must include columns \code{CONTIG}, \code{POSITION}, \code{ALT_COUNT}, \code{REF_COUNT}, \code{REF_NUCLEOTIDE}, and \code{ALT_NUCLEOTIDE}.
#' @return Invisibly returns the output file path.
#' @importFrom data.table fread
#' @importFrom dplyr filter mutate rowwise ungroup select
#' @importFrom stringr str_split
#' @importFrom tidyr unnest_wider
#' @export
RunAIsegmentation <- function(
    seg, cov, ai, gender, out_dir, prefix,
    aibinsize = 500000,
    mergeai = 0.2,
    mergecov = 0.2,
    minaisize = 1000000,
    snpmin = 7,
    minsnpcallaicutoff = 10,
    mergecovminsize = 500000,
    aitype="count"
) {
  gender = tolower(gender)
  # Read BAF file eg. GATK
  maf <- ReadAI(aitype = aitype, ballele = ai )

  # Read coverage file
  lines <- readLines(cov)
  filtered_lines <- grep("^\\s*@", lines, invert = TRUE, value = TRUE)
  cov_df <- data.table::fread(text = filtered_lines)

  # Read segment file
  seg_df <- data.table::fread(seg) %>%
    dplyr::mutate(size = End - Start)

  # Step 1: Gender check
  seg_df <- CheckGender(cov = cov_df, seg = seg_df, gender = gender)
  merge_seg <- seg_df

  # Step 2: Merge GATK segments by size
  if ( mergecovminsize < 1e+05) mergecovminsize <- 1e+05
  for (minsize in c(1000, 10000, 50000, seq(from = 100000, to = as.numeric(mergecovminsize), by = 100000))) {
    merge_seg <- merge_seg %>%
      dplyr::ungroup() %>%
      dplyr::filter(size >= minsize | Chromosome %in% c("X", "Y"))
    merge_seg <- CallMerge(data = merge_seg, AIorSeg = "Seg", snpmin = snpmin, mergeai = mergeai, mergecov = mergecov)
  }
  merge_seg <- subset(merge_seg, select = -Call)
  merge_seg <- merge_seg %>% dplyr::mutate(BreakpointSource = "GATK")

 # Step 3: Correct MAF scale
  maf_scale <- 2
  while (maf_scale > 1.005) {
   maf_scale <- EstimateMAFshift(maf)
   if (!is.na(maf_scale) & maf_scale <= 1.2) {
      maf <- CorrectMAF(maf, maf_scale)
      maf_scale <- EstimateMAFshift(maf)
    } else {
      maf_scale <- 1
    }
  }

  # Step 4: Search for breakpoints in each segment
  combindAIseg <- lapply(seq_len(nrow(merge_seg)), function(x) {
    SearchBreakpoint(seg_row = merge_seg[x, ], maf = maf, mergeai = mergeai, snpmin = snpmin, minaisize = minaisize, aibinsize = aibinsize)
  })
  combined_result <- do.call(rbind, combindAIseg)

  # Step 5: Add quality tag
  result <- combined_result %>%
    dplyr::rowwise() %>%
    dplyr::mutate(FILTER = AddQualTag(
      MAF_gmm_weight = MAF_gmm_weight,
      MAF_Probes = MAF_Probes,
      MAF_gmm_G = MAF_gmm_G,
      snpmin = snpmin
    ))
  if( gender == "female"){
    result <- result %>%
      dplyr::filter( Chromosome != "Y")
  }
  # Step 6: Write output
  outFile <- file.path(out_dir, paste0(prefix, "_GATK_AI_segment.tsv"))
  write.table(result, file = outFile, sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(outFile)
}
