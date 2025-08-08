#' Rerun CNV Calling with User-Defined Purity and Size Factor or Diploid Region
#'
#' Allows users to rerun copy number variant (CNV) calling with custom purity and scale factor (size factor) values, or by specifying a diploid region for normalization.
#'
#' When \code{mode = "region"}, \code{chromosome} must be set. If only \code{chromosome} and \code{start} are provided, the region spans from \code{start} to the end of the chromosome. If only \code{chromosome} and \code{end} are provided, the region spans from the start of the chromosome to \code{end}. If both \code{start} and \code{end} are given, the user-defined region is used.
#'
#' When \code{mode = "model"}, at least one of \code{dicovsf} or \code{purity} must be set. Both can be specified as single values or as ranges (e.g., "0.9:1.1").
#'
#' @param seg Character. Path to the AI segment file (e.g., \code{*_GATK_AI_segment.tsv}).
#' @param input Character. Path to the top likelihood row file (e.g., \code{*_top_likelihood_calls.tsv}).
#' @param models Character. Path to the model likelihood file (e.g., \code{*_Models_likelihood.tsv}).
#' @param call Character. Path to the final call file.
#' @param gender Character. Sample gender, either \code{"male"} or \code{"female"}.
#' @param dicovsf Numeric or character. Desired scale factor or range (e.g., \code{0.9}, \code{"0.9:1.1"}). Optional.
#' @param purity Numeric or character. Desired purity or range (e.g., \code{0.5}, \code{"0.5:0.7"}). Optional.
#' @param callcov Numeric. Subclonal events calling cutoff based on coverage.
#' @param chromosome Character. Chromosome for diploid region (e.g., \code{"1"}, \code{"X"}). Optional.
#' @param start Integer. Start position for diploid region. Optional.
#' @param end Integer. End position for diploid region. Optional.
#' @param mode Character. Rerun mode: either \code{"model"} or \code{"region"}.
#' @param out_file Character. Output file path for the final calls.
#'
#' @return Invisibly returns the final call data frame.
#'
#' @details
#' This function enables re-calling of CNV segments using user-specified purity and/or size factor (SF), or by defining a diploid region for normalization. It supports both "model" and "region" rerun modes. The resulting calls are written to \code{out_file}.
#'
#' @importFrom data.table fread
#' @importFrom dplyr rowwise mutate ungroup arrange filter select
#' @importFrom utils write.table
#'
#' @examples
#' \dontrun{
#' # Rerun using a specific purity and size factor
#' RerunCNV(
#'   seg = "sample_GATK_AI_segment.tsv",
#'   input = "sample_top_likelihood_calls.tsv",
#'   models = "sample_Models_likelihood.tsv",
#'   call = "sample_final_call.tsv",
#'   gender = "female",
#'   dicovsf = "0.95:1.05",
#'   purity = "0.6:0.8",
#'   mode = "model",
#'   out_file = "sample_final_call_refined.tsv"
#' )
#'
#' # Rerun using a user-defined diploid region
#' RerunCNV(
#'   seg = "sample_GATK_AI_segment.tsv",
#'   input = "sample_top_likelihood_calls.tsv",
#'   models = "sample_Models_likelihood.tsv",
#'   call = "sample_final_call.tsv",
#'   gender = "male",
#'   chromosome = "3",
#'   start = 1e6,
#'   end = 5e7,
#'   mode = "region",
#'   out_file = "sample_final_call_refined.tsv"
#' )
#' }
#'
#' @export
RerunCNV <- function(
    seg, input, models, call, gender, callcov = 0.3,
    dicovsf = NULL, purity = NULL,
    chromosome = NULL, start = NULL, end = NULL,
    mode = NULL,
    out_file
) {
  # Read data
  top_likelihood_rows <- data.table::fread(input)
  top_likelihood_rows$index <- as.character(top_likelihood_rows$index)
  models_df <- data.table::fread(models)
  call_df <- data.table::fread(call)
  call_df$index <- as.character(call_df$index)
  seg_df <- data.table::fread(seg)
  diploid_cov <- 100
  chr_levels <- c(as.character(1:22), "X", "Y")
  seg_df <- seg_df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Segcov = SegmentMeanToOriCov(
      gender = gender,
      chromosome = Chromosome,
      diploid_cov = diploid_cov,
      SM = Segment_Mean
    )) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(Chromosome = factor(Chromosome, levels = chr_levels)) %>%
    dplyr::arrange(Chromosome, Start) %>%
    dplyr::mutate(index = as.character(row_number()))

  # Print parameters

  print("Rerun Parameters are:")
  print(paste0("Gender:", gender))
  print(paste0("mode: ", mode))
  if(! is.null(dicovsf)){
    print(paste0("dicovsf: ", dicovsf))
  }
  if(! is.null(purity)){
    print(paste0("purity: ", purity))
  }
  if( ! is.null( chromosome)){
    print(paste0("chromosome: ", chromosome))
  }
  if(! is.null(start) ){
    print(paste0("start: ", start))
  }
  if(! is.null(end)){
    print(paset0("end: ", end ))
  }



  # Check parameters
  updated_parameters <- Checkmode(mode = mode, purity = purity, dicovsf = dicovsf, chromosome = chromosome, call = call_df, start = start, end = end)
  dicovsf <- updated_parameters$dicovsf
  purity <- updated_parameters$purity
  mode <- updated_parameters$mode
  # Main: parse parameters and rerun
  final_model <- ParseParm(purity = purity, dicovsf = dicovsf, models = models_df)

  # Extract final call
  raw_call <- ExtractCall(
    df = top_likelihood_rows,
    max_L_mu = final_model$final_sf,
    max_L_rho = final_model$final_purity,
    seg = seg_df
  )

  # Recalculate all values
  refined_call <- RefineCalls(
    df = raw_call,
    max_L_mu = final_model$final_sf,
    max_L_rho = final_model$final_purity,
    gender = gender
  )

  # Refine mismatched coverage/AI profiles
  final_call <- RefineCallsSecond(
    df = refined_call,
    results = call_df,
    final_mu = final_model$final_sf,
    final_rho = final_model$final_purity,
    gender = gender,
    callcov = callcov
  )
  final_call$Model_source <- "User_defined"

  # Write results
  utils::write.table(final_call, file = out_file, row.names = FALSE, quote = FALSE, sep = "\t")
  invisible(final_call)
}
