#' Check and Set Mode-Specific Options
#'
#' Checks required options for 'model' or 'region' mode and sets defaults if needed.
#'
#' @param opt List of options. Must include \code{mode} (either "model" or "region"), and possibly \code{purity}, \code{dicovsf}, \code{chromosome}.
#' @param call Data frame or tibble. Segment call data, required if estimating coverage scale factor for region mode.
#'
#' @return The (potentially modified) \code{opt} list.
#'
#' @details
#' - If \code{mode == "model"}: requires at least one of \code{purity} or \code{dicovsf}.
#' - If \code{mode == "region"}: requires \code{chromosome}; if \code{dicovsf} is missing, it is estimated from \code{call}.
#'
#' @importFrom rlang is_null
#' @examples
#' # opt <- list(mode = "model", purity = NULL, dicovsf = NULL)
#' # Checkmode(opt, call)
#'
#' @export

Checkmode <- function(opt, call ){

  if(opt$mode == "model"){

    if( is.null( opt$purity) && is.null( opt$dicovsf )  ){

      stop("Error: Under 'model' mode, please provide at least one of 'purity' or 'dicovsf'.")
    }
  }else if( opt$mode == "region"){

    if( is.null( opt$chromosome ) ){
      stop("Error: Under 'region' mode, please provide at least chromosome without 'chr'.")
    }else{
      if( is.null(opt$dicovsf) ){
        opt$dicovsf <- EstimateCovSF( seg = call , opt = opt )
      }
    }

  }
  return(opt)
}


#' Estimate Coverage Scale Factor for a Chromosome or Region
#'
#' Estimates the coverage scale factor (covsf) for a specified chromosome or region, based on segment mean and probe count.
#'
#' @param opt List of options. Must include \code{chromosome}; may include \code{start} and \code{end} (numeric).
#' @param seg Data frame or tibble. Segment data, must include columns \code{Chromosome}, \code{Start}, \code{End}, \code{Segment_Mean}, and \code{Num_Probes}.
#'
#' @return Numeric. The estimated coverage scale factor (covsf) for the specified region.
#'
#' @details
#' - If only \code{chromosome} is defined, uses the entire chromosome.
#' - If \code{start} is defined, uses regions from \code{start} to the end (or to \code{end} if provided).
#' - If \code{end} is defined, uses regions from chromosome start to \code{end}.
#'
#' @importFrom dplyr filter mutate
#' @examples
#' # opt <- list(chromosome = "3", start = NULL, end = NULL)
#' # EstimateCovSF(opt, seg)
#'
#' @export
EstimateCovSF <- function( opt, seg ){

  diploid_cov <- 100
  tmp <- seg %>% dplyr::filter( Chromosome == opt$chromosome )
  if( is.null(opt$start)){
    start <- min(tmp$Start) }else{start <- as.numeric(opt$start)}
  if( is.null(opt$end)){
    end <- min(tmp$End)  }else{ end <- as.numeric(opt$end) }

  covsf_tmp <- tmp %>% dplyr::mutate( select_start_index = ifelse( Start>= start, 1, 0),
                               select_end_index = ifelse( End <= end , 1, 0 )) %>%
    filter( select_start_index ==1 & select_end_index == 1)
  total_size <- sum( covsf_tmp$Num_Probes )
  Segment_cov <- 2^( covsf_tmp$Segment_Mean) * diploid_cov
  dicovsf <-  sum( Segment_cov * ( covsf_tmp$Num_Probes / total_size)) / diploid_cov

  return(dicovsf)
}

#' Parse Purity and Scale Factor Parameters for Model Selection
#'
#' Parses user-specified purity and scale factor (dicovsf) options, selects matching models, and returns the best match.
#'
#' @param opt List of options. May include \code{purity} (numeric or range as "min:max") and \code{dicovsf} (numeric or range as "min:max").
#' @param models Data frame or tibble. Must include columns \code{rho}, \code{mu}, \code{total_distance_to_integer}, and \code{total_log_likelihood}.
#'
#' @return A list with:
#'   \item{final_purity}{Selected purity value.}
#'   \item{final_sf}{Selected scale factor (mu) value.}
#'
#' @details
#' - If \code{opt$purity} or \code{opt$dicovsf} are ranges (e.g., "0.6:0.8"), selects all models in the range.
#' - Otherwise, selects the model(s) closest to the specified value.
#' - Returns the model with minimal distance to integer and maximal likelihood.
#'
#' @importFrom dplyr mutate arrange filter
#' @importFrom stringr str_split
#' @examples
#' # opt <- list(purity = "0.6:0.8", dicovsf = NULL)
#' # ParseParm(opt, models)
#'
#' @export
ParseParm <- function( opt, models ){

  if( !is.null(opt$purity )){
    if( grepl(":", opt$purity ) ){
      range <- strsplit(x = opt$purity, split = ":") %>% unlist() %>% as.numeric()
      purity <- models[which(models$rho >= range[1] & models$rho <= range[2]),"rho"]
    }else{
      purity <- models %>% mutate( diff = abs( opt$purity - rho ) ) %>%
        dplyr::arrange(diff)
      purity <- purity[1,"rho"]
    } }else{ purity <- unique(models$rho) }

  if( !is.null(opt$dicovsf)){
    if( grepl(":", opt$dicovsf ) ){
      range <- strsplit(x = opt$dicovsf, split = ":") %>% unlist() %>% as.numeric()
      dicovsf <- models[which(models$mu >= range[1] & models$mu <= range[2]),"mu"]
    }else{
      dicovsf <- models %>% dplyr::mutate( diff = abs( opt$dicovsf - mu ) ) %>% dplyr::arrange(diff)
      dicovsf <- dicovsf[1,"mu"]
    }
  }else{ dicovsf <- unique(models$mu)}
  final_model <- models %>%
    filter( mu >= min(dicovsf) & mu <= max(dicovsf)) %>%
    filter( rho >= min(purity) & rho <= max(purity) ) %>%
    arrange( total_distance_to_integer, -total_log_likelihood   )
  if( nrow(final_model) > 0 ){
    final_model <- list( final_purity = final_model$rho[1],
                         final_sf = final_model$mu[1])
    return(final_model)
  }else{
    stop("Error: The combination is not possible, please choose another one.")

  }


}
