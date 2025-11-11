#' Check and Set Mode-Specific Parameters for CNV Calling
#'
#' Checks and sets required parameters for "model" or "region" mode in CNV calling, ensuring user input is sufficient and filling in missing values when possible.
#'
#' @param mode Character. Mode of operation; either \code{"model"} or \code{"region"}.
#' @param purity Numeric or NULL. Tumor purity value or range (required for "model" mode if \code{dicovsf} is not provided).
#' @param dicovsf Numeric or NULL. Size factor (diploid coverage scale factor) or range (required for "model" mode if \code{purity} is not provided).
#' @param chromosome Character or NULL. Chromosome number for diploid region (required for "region" mode).
#' @param start Numeric or NULL. Start coordinate of the diploid region. If not set in "region" mode, the program will use the first nucleotide of the detectable region of the chromosome.
#' @param end Numeric or NULL. End coordinate of the diploid region. If not set in "region" mode, the program will use the last nucleotide of the detectable region of the chromosome.
#' @param call Data frame or tibble. Segment call data, required if estimating \code{dicovsf} in "region" mode.
#'
#' @return A named list with updated parameters: \code{mode}, \code{dicovsf}, and \code{purity}.
#'
#' @details
#' \itemize{
#'   \item In "model" mode, at least one of \code{purity} or \code{dicovsf} must be provided.
#'   \item In "region" mode, \code{chromosome} must be provided; if \code{dicovsf} is missing, it is estimated using \code{EstimateCovSF()} with the specified region.
#' }
#'
#' @examples
#' Checkmode("model", purity = 0.7, dicovsf = NULL, chromosome = NULL, start = NULL, end = NULL, call = call_df)
#' Checkmode("region", purity = NULL, dicovsf = NULL, chromosome = "1", start = 100000, end = 200000, call = call_df)
#'
#' @export

Checkmode <- function(mode, purity, dicovsf, chromosome, start, end, call ){

  if( mode == "model"){

    if( is.null( purity ) && is.null( dicovsf )  ){

      stop("Error: Under 'model' mode, please provide at least one of 'purity' or 'dicovsf'.")
    }
  }else if( mode == "region"){

    if( is.null( chromosome ) ){
      stop("Error: Under 'region' mode, please provide at least chromosome without 'chr'.")
    }else{
      if( is.null( dicovsf) ){
        dicovsf <- EstimateCovSF( seg = call , chromosome = chromosome, start = start, end = end )
      }
    }

  }
  updated_parameters <- list(
    mode = mode,
    dicovsf = dicovsf,
    purity = purity
  )
  return( updated_parameters )
}


#' Estimate Coverage Scale Factor for a Chromosome or Region
#'
#' Estimates the coverage scale factor (covsf) for a specified chromosome or subregion, based on segment mean and probe count.
#'
#' @param chromosome Character. Chromosome name (e.g., "1", "X").
#' @param start Numeric or NULL. Start coordinate of the region. If \code{NULL}, uses the first segment start in the chromosome.
#' @param end Numeric or NULL. End coordinate of the region. If \code{NULL}, uses the last segment end in the chromosome.
#' @param seg Data frame or tibble. Segment data, must include columns \code{Chromosome}, \code{Start}, \code{End}, \code{Segment_Mean}, and \code{Num_Probes}.
#'
#' @return Numeric. The estimated coverage scale factor (covsf) for the specified region.
#'
#' @details
#' The function calculates a weighted mean coverage based on the segment means and probe counts within the specified region. If \code{start} or \code{end} are not provided, the entire chromosome is used for estimation.
#'
#' @importFrom dplyr filter mutate
#' @examples
#' # EstimateCovSF("3", start = NULL, end = NULL, seg = seg_df)
#' # EstimateCovSF("7", start = 100000, end = 500000, seg = seg_df)
#'
#' @export

EstimateCovSF <- function( chromosome, start, end, seg ){
  raw_mu <- as.numeric(seg$mu[1])
  diploid_cov <- 100
  tmp <- seg %>% dplyr::filter( Chromosome == chromosome )
  if( is.null(start)){
    print("Start point is not defined in region mode; using the minimum detectable boundary of the chromosome as the end position.")
    start <- min(tmp$Start) }else{start <- as.numeric(start)}
  if( is.null(end)){
    print("End point is not defined in region mode; using the maximum detectable boundary of the chromosome as the end position.")
    end <- min(tmp$End)  }else{ end <- as.numeric(end) }

  covsf_tmp <- tmp %>% dplyr::mutate( select_start_index = ifelse( Start>= start, 1, 0),
                               select_end_index = ifelse( End <= end , 1, 0 )) %>%
    filter( select_start_index ==1 & select_end_index == 1)
  total_size <- sum( covsf_tmp$Num_Probes )
  Segment_cov <- 2^( covsf_tmp$Segment_Mean) * diploid_cov
  dicovsf <-  sum( Segment_cov * ( covsf_tmp$Num_Probes / total_size)) / diploid_cov # this is sf relative to current model, which is raw_mu but we need new size factor relative to raw GATK
  dicovsf <- raw_mu * dicovsf

  return(dicovsf)
}

#' Parse and Select Purity and Scale Factor Parameters from Model Grid
#'
#' Selects the best-matching purity and scale factor (dicovsf) from a grid of models, based on user input (value or range), and returns the optimal combination for downstream analysis.
#'
#' @param purity Numeric, character, or NULL. Desired purity value, or a range in the format \code{"min:max"} (e.g., \code{"0.5:0.7"}). If \code{NULL}, all available purity values in \code{models} are considered.
#' @param dicovsf Numeric, character, or NULL. Desired scale factor (mu), or a range in the format \code{"min:max"} (e.g., \code{"0.9:1.1"}). If \code{NULL}, all available scale factor values in \code{models} are considered.
#' @param models Data frame or tibble. Model grid with columns \code{rho} (purity), \code{mu} (scale factor), \code{total_distance_to_integer}, and \code{total_log_likelihood} generated by initial run.
#'
#' @return A named list with elements:
#'   \item{final_purity}{Selected purity value.}
#'   \item{final_sf}{Selected scale factor (mu) value.}
#'
#' @details
#' - If a range is provided for \code{purity} or \code{dicovsf}, all models within that range are considered.
#' - If a single value is provided, the closest available value in \code{models} is selected.
#' - The model with the lowest \code{total_distance_to_integer} and highest \code{total_log_likelihood} is chosen.
#'
#' @examples
#' # ParseParm(purity = "0.5:0.7", dicovsf = NULL, models = models_df)
#' # ParseParm(purity = 0.6, dicovsf = 1.1, models = models_df)
#'
#' @export
ParseParm <- function( purity, dicovsf, models ){

  if( !is.null( purity )){
    if( grepl(":", purity ) ){
      range <- strsplit(x = purity, split = ":") %>% unlist() %>% as.numeric()
      purity <- models[which(models$rho >= range[1] & models$rho <= range[2]),"rho"]
    }else{
      purity <- models %>% mutate( diff = abs( purity - rho ) ) %>%
        dplyr::arrange(diff)
      purity <- purity[1,"rho"]
    } }else{ purity <- unique(models$rho) }

  if( !is.null(dicovsf)){
    if( grepl(":", dicovsf ) ){
      range <- strsplit(x = dicovsf, split = ":") %>% unlist() %>% as.numeric()
      dicovsf <- models[which(models$mu >= range[1] & models$mu <= range[2]),"mu"]
    }else{
      dicovsf <- models %>% dplyr::mutate( diff = abs( dicovsf - mu ) ) %>% dplyr::arrange(diff)
      dicovsf <- dicovsf[1,"mu"]
    }
  }else{ dicovsf <- unique(models$mu)}
  final_model <- models %>%
    dplyr::filter( mu >= min(dicovsf) & mu <= max(dicovsf)) %>%
    dplyr::filter( rho >= min(purity) & rho <= max(purity) ) %>%
    dplyr::arrange( total_likelihood_cluster,
                    diploid_distance_cluster ,
                    nondiploid_distance_cluster,
                    abs( ploidy - 2),
                    desc(total_log_likelihood_after_refine), desc( rho ) )
  if( nrow(final_model) > 0 ){
    final_model <- list( final_purity = final_model$rho[1],
                         final_sf = final_model$mu[1])
    return(final_model)
  }else{
    stop("Error: The combination is not possible, please choose another one.")

  }


}
