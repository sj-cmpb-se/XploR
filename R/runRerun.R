#' Rerun CNV Calling with User-Defined Purity and Size Factor
#'
#' Allows the user to rerun copy number variant (CNV) calling using user-defined purity and scale factor (size factor) parameters, or by specifying a diploid region.
#'
#' @param seg Path to the segment file (e.g., *_final_call.tsv).
#' @param input Path to the top likelihood row file.
#' @param models Path to the model file.
#' @param call Path to the final call file.
#' @param gender Character. Sample gender ("male" or "female").
#' @param dicovsf Character or numeric. Desired scale factor, or a range (e.g., "0.9:1.1").
#' @param purity Character or numeric. Desired purity, or a range (e.g., "0.5:0.7").
#' @param chromosome Character. Chromosome number for diploid region (optional).
#' @param start Integer. Start position for diploid region (optional).
#' @param end Integer. End position for diploid region (optional).
#' @param mode Character. Either "model" or "region".
#' @param out_file Output file path for the final calls.
#'
#' @return Invisibly returns the final call data frame.
#'
#' @details
#' This function allows re-calling CNV segments with user-specified purity and/or size factor (SF, scale factor) or by defining a diploid region.
#' It supports both "model" and "region" rerun modes. The result is written to \code{out_file}.
#'
#' @importFrom data.table fread
#' @importFrom dplyr rowwise mutate ungroup arrange filter select
#' @importFrom utils write.table
#' @examples
#' # rerun_cnv_calling(
#' #   seg = "sample_final_call.tsv",
#' #   input = "sample_top_likelihood_calls.tsv",
#' #   models = "sample_Models_likelihood.tsv",
#' #   call = "sample_final_calls.tsv",
#' #   gender = "male",
#' #   purity = "0.5:0.7",
#' #   dicovsf = "0.9:1.1",
#' #   mode = "model",
#' #   out_file = "sample_final_calls_rerun.tsv"
#' # )
#' @export
rerunCalling <- function(
    seg, input, models, call, gender,
    dicovsf = NULL, purity = NULL,
    chromosome = NULL, start = NULL, end = NULL,
    mode = "model",
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

  # Prepare options list for helper functions
  opt <- list(
    gender = gender,
    dicovsf = dicovsf,
    purity = purity,
    chromosome = chromosome,
    start = start,
    end = end,
    mode = mode
  )

  # Check parameters
  opt <- Checkmode(opt = opt, call = call_df)

  # Main: parse parameters and rerun
  final_model <- ParseParm(opt, models_df)

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
    opt = opt
  )

  # Refine mismatched coverage/AI profiles
  final_call <- RefineCallsSecond(
    df = refined_call,
    results = call_df,
    final_mu = final_model$final_sf,
    final_rho = final_model$final_purity,
    opt = opt
  )
  final_call$Model_source <- "User_defined"

  # Write results
  utils::write.table(final_call, file = out_file, row.names = FALSE, quote = FALSE, sep = "\t")
  invisible(final_call)
}
