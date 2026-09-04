#' Read and Standardize Allelic Imbalance Data
#'
#' Reads allelic imbalance (AI) data from different input formats and standardizes the output.
#'
#' @param aitype Character. Type of AI input file. Supported values are:
#'   \itemize{
#'     \item \code{"gatk"}: Tab-separated file with columns \code{Chromosome}, \code{Pos}, \code{REF_NUCLEOTIDE}, \code{ALT_NUCLEOTIDE}, \code{REF_COUNT}, \code{ALT_COUNT}.
#'     \item \code{"dragen"}: the input must include columns \code{contig}, \code{start}, \code{stop}, \code{refAllele}, \code{allele1}, \code{allele2}, \code{allele1Count}, \code{allele2Count}, \code{allele1AF}, and \code{allele2AF}.
#'     \item \code{"other"}: Tab-separated file with columns \code{Chromosome}, \code{Pos}, \code{REF_NUCLEOTIDE}, \code{ALT_NUCLEOTIDE}, \code{REF_COUNT}, \code{ALT_COUNT}.
#'   }
#' @param ballele Character. Path to the input allelic count file.
#' @param minsnpcov Numeric, minimum coverage of SNPs to included (default: 20).
#' @param gender Character. Sample gender, either \code{"female"} or \code{"male"}. This value is passed to \code{ReadAI()}.
#'
#' @return A data frame with columns \code{Chromosome}, \code{Pos}, \code{ref_count}, \code{alt_count} and \code{maf}, containing only autosomal positions (\code{Chromosome} not X or Y), with minor allele frequency (\code{maf}) between 0 and 0.5.
#'
#' @examples
#' \dontrun{
#' # Read AI data from a count file
#' maf_df <- ReadAI("gatk", "sample.allele_counts", 20, "female")
#'
#' # Read AI data from a DRAGEN file
#' maf_df <- ReadAI("dragen", "sample_dragen.tumor.ballele.counts.gz",20, "female")
#'
#' # Read AI data from a BAF file
#' maf_df <- ReadAI("other", "sample.tsv",20, "female")
#' }
#'
#' @importFrom dplyr filter mutate select
#' @export
ReadAI <- function(aitype, ballele, minsnpcov, gender){
  if( aitype == "gatk"){
    filtered_lines <- grep("^[^@]", readLines(ballele), value = TRUE)
    baf <- read.table(text = paste(filtered_lines, collapse = "\n"), sep = "\t", header = TRUE)
    baf <- baf %>%
      dplyr::filter( REF_COUNT + ALT_COUNT >= minsnpcov ) %>%
      dplyr::mutate( maf = pmin(ALT_COUNT, REF_COUNT)/(ALT_COUNT + REF_COUNT)) %>%
      dplyr::select( -REF_NUCLEOTIDE, -ALT_NUCLEOTIDE ) %>%
      dplyr::relocate( REF_COUNT, ALT_COUNT, maf,.after = dplyr::last_col())

  }else if( aitype == "dragen" ){
    baf <- read.table(ballele, sep = "\t", header = TRUE,stringsAsFactors = F)
    baf <- baf %>%
      dplyr::filter( allele1Count + allele2Count >= minsnpcov ) %>%
      dplyr::mutate( maf = pmin(allele2Count, allele1Count)/(allele2Count + allele1Count))
    baf <- baf[,c("contig","start","allele1Count","allele2Count","maf")]
  }else{
    baf <- read.table(ballele, sep = "\t", header = TRUE,stringsAsFactors = F)
    baf %>%
      dplyr::filter( REF_COUNT + ALT_COUNT >= minsnpcov ) %>%
      dplyr::mutate( maf = pmin(ALT_COUNT, REF_COUNT)/(ALT_COUNT + REF_COUNT))
    baf <- baf[,c("Chromosome","Pos","REF_COUNT","ALT_COUNT","maf")]
  }
  colnames(baf) <- c("Chromosome","Pos","ref_count","alt_count","maf")
  if( gender == "male"){
    maf <- baf %>%
      dplyr::filter(!Chromosome %in% c("X", "Y"))
  }else{
    maf <- baf %>%
      dplyr::filter( Chromosome != "Y")
  }
  return(maf)
}


#' Check Whether Two AI Segments Should Be Merged
#'
#' This function determines whether two adjacent AI (allelic imbalance) segments should be merged,
#' based on the difference in their GMM means and quality thresholds.
#'
#' @param cur_row A data frame row (as a list or tibble row) representing the current segment. Must include relevant columns (e.g., \code{gmm_mean}, \code{gmm_weight}, \code{nonzero_count}).
#' @param next_row A data frame row (as a list or tibble row) representing the next segment. Must include relevant columns (e.g., \code{gmm_mean}, \code{gmm_weight}, \code{nonzero_count}).
#' @param mergeai Numeric. Threshold for the difference in MAF (gmm_mean) between adjacent segments to allow merging.
#' @param snpmin Numeric. Minimum SNP count required for a segment to be considered as a separate segment.
#'
#' @return Logical value: \code{TRUE} if the two segments should be merged, \code{FALSE} otherwise.
#'
#' @export

MergeAICheck <- function(cur_row,next_row, mergeai, snpmin){
  # Merge if
  # ai segment diff <= mergeai or one of the ai segment quality is fail
  if(abs( as.numeric(cur_row$gmm_mean) - as.numeric(next_row$gmm_mean)) <= mergeai ){
    re <- TRUE
  }else if( !is.na(cur_row$gmm_weight) & !is.na(next_row$gmm_weight)){
    if(min(cur_row$gmm_weight , next_row$gmm_weight ) <= 0.2||
       min(cur_row$nonzero_count, next_row$nonzero_count) < 2*snpmin ){ re <- TRUE }else{re <- FALSE}
  } else{ re <- FALSE}
  return(re)
}


#' Cluster Adjacent Values Based on a Threshold
#'
#' Groups adjacent values into clusters if their difference is less than a specified threshold,
#' and computes the mean and total weight for each cluster.
#'
#' @param values Numeric vector of values to cluster (e.g., means from GMM components).
#' @param weights Numeric vector of weights corresponding to \code{values}.
#' @param threshold Numeric value specifying the maximum allowed difference between adjacent values to be grouped into the same cluster. Default is 0.01.
#'
#' @return A list with two elements:
#'   \item{means}{Numeric vector of cluster means.}
#'   \item{weights}{Numeric vector of aggregated cluster weights.}
#'
#' @examples
#' ClusterAdjacent(c(0.1, 0.11, 0.5), c(0.3, 0.4, 0.3))
#'
#' @export
ClusterAdjacent <- function(values, weights, threshold = 0.01) {
  # Initialize cluster ID vector
  cluster_ids <- integer(length(values))
  current_id <- 1
  cluster_ids[1] <- current_id

  # Assign cluster IDs based on adjacent differences
  for (i in 2:length(values)) {
    if (values[i] - values[i-1] < threshold) {
      cluster_ids[i] <- current_id
    } else {
      current_id <- current_id + 1
      cluster_ids[i] <- current_id
    }
  }

  # Calculate new means for each cluster
  cluster_means <- tapply(values, cluster_ids, mean)

  # Aggregate weights for each cluster
  cluster_weights <- tapply(weights, cluster_ids, sum)

  # Return results as list
  list(means = as.numeric(cluster_means),
       weights = as.numeric(cluster_weights))
}

#' Estimate Tumor MAF from BAF Using Normal Controls
#'
#' Estimates tumor minor allele frequency (MAF) from raw B-allele frequency
#' (BAF) values using one- and two-component Gaussian models. Region-matched
#' normal BAF values are used to define an empirical BIC cutoff for technical
#' variation.
#'
#' Tumor and normal BAF values are matched to the same number of observations.
#' A two-peak tumor model is accepted only when its BIC gain exceeds the normal
#' cutoff and the fitted peaks meet minimum weight, separation, and center
#' requirements. Otherwise, the region is classified as balanced with MAF = 0.5.
#'
#' @param tumor_baf_values Numeric vector of tumor BAF values.
#' @param normal_baf_values Numeric vector of region-matched normal BAF values.
#' @param baf_lower Lower BAF bound used for fitting. Default 0.20.
#' @param baf_upper Upper BAF bound used for fitting. Default 0.80.
#' @param min_n Minimum number of usable BAF values. Default 50.
#' @param n_resample Number of matched subsampling iterations. Default 100.
#' @param normal_bic_quantile Quantile of normal BIC gain used as cutoff.
#'   Default 1 (maximum).
#' @param bic_margin Additional BIC gain required above the normal cutoff.
#' @param min_component_weight Minimum weight of each two-peak component.
#' @param min_separation_sd Minimum peak separation in units of fitted SD.
#' @param center_lower Lower allowed midpoint of the two peaks.
#' @param center_upper Upper allowed midpoint of the two peaks.
#' @param seed Random seed for reproducible subsampling.
#'
#' @return Named numeric vector containing:
#'   \item{baf_maf}{Final estimated MAF.}
#'   \item{baf_G}{Selected number of peaks, 1 or 2.}
#'   \item{baf_peak1, baf_peak2}{Estimated BAF peak locations.}
#'   \item{baf_weight1, baf_weight2}{Estimated component weights.}
#'   \item{baf_center}{Midpoint of the fitted peaks.}
#'   \item{baf_d}{Distance from the midpoint to each peak.}
#'   \item{baf_sigma}{Common fitted component SD.}
#'   \item{baf_separation_sd}{Peak separation divided by fitted SD.}
#'   \item{baf_bic_gain}{Median tumor BIC gain from matched subsampling.}
#'   \item{baf_bic_gain_full}{BIC gain using all tumor BAF values.}
#'   \item{normal_bic_cutoff}{Empirical normal BIC cutoff.}
#'   \item{normal_bic_median}{Median normal BIC gain.}
#'   \item{baf_bic_excess}{Tumor BIC gain minus normal cutoff.}
#'   \item{baf_n_tumor, baf_n_normal, baf_n_match}{Numbers of BAF values used.}
#'
#' @details
#' The two-component model constrains peaks to \code{center - d} and
#' \code{center + d}. MAF is estimated as \code{0.5 - d}. Normal controls
#' provide a region-specific empirical baseline for apparent bimodality.
#'
#' @importFrom stats dnorm mean median optim plogis qlogis quantile sample sd var
#' @export

EstimateMAFfromBAF <- function(
    tumor_bafs,
    normal_bafs,
    baf_lower = 0.20,
    baf_upper = 0.80,
    min_n = 10L,
    n_resample = 100L,
    normal_bic_quantile = 1,
    bic_margin = 0,
    min_component_weight = 0.25,
    min_separation_sd = 2,
    center_lower = 0.42,
    center_upper = 0.58,
    seed = 1L
) {

  # ============================================================
  # Prepare BAF values
  # ============================================================

  tumor_baf <- as.numeric(tumor_bafs)
  normal_baf <- as.numeric(normal_bafs)

  tumor_baf <- tumor_baf[ is.finite(tumor_baf) &  tumor_baf > baf_lower & tumor_baf < baf_upper ]
  normal_baf <- normal_baf[ is.finite(normal_baf) & normal_baf > baf_lower & normal_baf < baf_upper]

  n_tumor <- length(tumor_baf)
  n_normal <- length(normal_baf)

  n_match <- min( n_tumor, n_normal )


  # ============================================================
  # Empty result
  # ============================================================

  empty_result <- c(
    "baf_maf" = NA_real_,
    "baf_G" = NA_real_,
    "baf_peak1" = NA_real_,
    "baf_peak2" = NA_real_,
    "baf_weight1" = NA_real_,
    "baf_weight2" = NA_real_,
    "baf_center" = NA_real_,
    "baf_d" = NA_real_,
    "baf_sigma" = NA_real_,
    "baf_separation_sd" = NA_real_,
    "baf_bic_gain" = NA_real_,
    "baf_bic_gain_full" = NA_real_,
    "normal_bic_cutoff" = NA_real_,
    "normal_bic_median" = NA_real_,
    "baf_bic_excess" = NA_real_,
    "baf_n_tumor" = n_tumor,
    "baf_n_normal" = n_normal,
    "baf_n_match" = n_match
  )


  if ( n_tumor < min_n || n_normal < min_n ||n_match < min_n) {
    return(empty_result)
  }


  # ============================================================
  # Internal function:
  # fit one-peak and constrained two-peak models
  # ============================================================

  fit_baf_models <- function(baf) {

    baf <- as.numeric(baf)
    baf <- baf[is.finite(baf)]

    n_baf <- length(baf)

    if (
      n_baf < min_n ||
      !is.finite(stats::sd(baf)) ||
      stats::sd(baf) < 1e-6
    ) {
      return(NULL)
    }


    # ----------------------------------------------------------
    # H0: one Gaussian peak
    # ----------------------------------------------------------

    mu0 <- mean(baf)
    # MLE estimate of sigma
    sigma0 <- sqrt( mean((baf - mu0)^2))
    sigma0 <- max(sigma0, 0.005)

    logLik_one <- sum( stats::dnorm( baf, mean = mu0, sd = sigma0, log = TRUE))

    # Two parameters:
    # mean + sigma
    bic_one <- -2 * logLik_one + 2 * log(n_baf)


    # ----------------------------------------------------------
    # H1: constrained two-component Gaussian mixture
    #
    # Peak 1 = center - d
    # Peak 2 = center + d
    #
    # Both peaks share sigma.
    # Weights can differ.
    # ----------------------------------------------------------

    negative_log_likelihood <- function(par) {

      center <- par[1]
      d <- par[2]
      sigma <- exp(par[3])

      weight1 <- stats::plogis( par[4] )
      weight2 <- 1 - weight1

      peak1 <- center - d
      peak2 <- center + d
      log_density1 <- log(weight1) + stats::dnorm( baf, mean = peak1, sd = sigma, log = TRUE)
      log_density2 <- log(weight2) + stats::dnorm( baf, mean = peak2,sd = sigma,log = TRUE)

      # Stable log-sum-exp
      max_log_density <- pmax( log_density1,log_density2)
      log_mixture_density <-  max_log_density + log( exp(log_density1 - max_log_density) + exp(log_density2 - max_log_density) )

      negative_log_likelihood_value <- -sum(log_mixture_density)

      if (!is.finite(negative_log_likelihood_value)) {
        return(1e100)
      }

      return(negative_log_likelihood_value)
    }


    # ----------------------------------------------------------
    # Multiple starting values
    # ----------------------------------------------------------

    starting_center <- stats::median(baf)
    starting_center <- min(max(starting_center, 0.45),0.55)
    d_starts <- c( 0.01,0.02,0.03,0.04,0.05,0.07, 0.10 )
    weight_starts <- c(0.35, 0.50,0.65)
    fits <- list()
    counter <- 1L

    for (starting_d in d_starts) {
      for (starting_weight in weight_starts) {
        starting_sigma <- sqrt( max( stats::var(baf) - starting_d^2, 0.01^2 ) )
        starting_sigma <- min( max( starting_sigma, 0.005), 0.20 )
        fit <- tryCatch(
          stats::optim( par = c(
            starting_center,
            starting_d,
            log(starting_sigma),
            stats::qlogis( starting_weight )
          ),
          fn = negative_log_likelihood,
          method ="L-BFGS-B",
          lower = c( 0.45, 0.001, log(0.005), stats::qlogis(0.05)),
          upper = c( 0.55, 0.15, log(0.20), stats::qlogis(0.95)),
          control = list( maxit = 1000 )),
          error = function(e) { message( "optim error: ",conditionMessage(e))
            NULL }
        )


        if ( !is.null(fit) && fit$convergence == 0 && is.finite(fit$value) ) {
          fits[[counter]] <- fit
          counter <- counter + 1L
        }}
    }

    if (length(fits) == 0) { return(NULL) }


    # ----------------------------------------------------------
    # Select best two-component solution
    # ----------------------------------------------------------

    fit_values <- vapply(fits, function(x) { x$value },numeric(1))
    best_fit <- fits[[which.min(fit_values)]]

    # ----------------------------------------------------------
    # Extract parameters
    # ----------------------------------------------------------

    center <- best_fit$par[1]
    d <- best_fit$par[2]
    sigma <- exp(best_fit$par[3] )
    weight1 <- stats::plogis(best_fit$par[4])
    weight2 <- 1 - weight1
    peak1 <- center - d
    peak2 <- center + d
    separation_sd <- (peak2 - peak1) /sigma


    # ----------------------------------------------------------
    # BIC for two-component model
    #
    # Four parameters:
    # center + d + sigma + weight
    # ----------------------------------------------------------

    logLik_two <- -best_fit$value
    bic_two <- -2 * logLik_two + 4 * log(n_baf)
    # Positive value favors two components
    bic_gain <- bic_one - bic_two


    return(
      list(
        mu0 = mu0,
        sigma0 = sigma0,
        bic_one = bic_one,
        bic_two = bic_two,
        bic_gain = bic_gain,
        peak1 = peak1,
        peak2 = peak2,
        weight1 = weight1,
        weight2 = weight2,
        center = center,
        d = d,
        sigma = sigma,
        separation_sd = separation_sd,
        n = n_baf
      )
    )
  }


  # ============================================================
  # Fit tumor BAF
  # ============================================================

  tumor_full_fit <- fit_baf_models( baf = tumor_baf)
  if (is.null(tumor_full_fit)) {
    return(empty_result)
  }
  tumor_bic_gain_full <- tumor_full_fit$bic_gain


  # ============================================================
  # Match tumor and normal SNP counts for BIC comparison.
  #
  # BIC depends on sample size, so tumor and normal must contain
  # the same number of BAF observations.
  # ============================================================

  n_resample <- max(1L,as.integer(n_resample))


  # Use fixed seed so clinical results are reproducible.
  if (!is.null(seed)) {
    set.seed(seed)
  }
  tumor_bic_values <- rep(NA_real_, n_resample )
  normal_bic_values <- rep(NA_real_,n_resample )
  for (i in seq_len(n_resample)) {
    # Tumor matched sample

    if (n_tumor == n_match) {
      tumor_sample <- tumor_baf
    } else {
      tumor_sample <- sample( tumor_baf, size = n_match, replace = FALSE)
    }

    # Normal matched sample
    if (n_normal == n_match) {
      normal_sample <- normal_baf
    } else { normal_sample <- sample( normal_baf, size = n_match, replace = FALSE ) }

    tumor_fit_i <- fit_baf_models(tumor_sample)
    normal_fit_i <-fit_baf_models(normal_sample)

    if (!is.null(tumor_fit_i)) { tumor_bic_values[i] <- tumor_fit_i$bic_gain }
    if (!is.null(normal_fit_i)) { normal_bic_values[i] <- normal_fit_i$bic_gain }
  }


  # ============================================================
  # Remove failed resamples
  # ============================================================

  tumor_bic_values <- tumor_bic_values[ is.finite( tumor_bic_values )]
  normal_bic_values <- normal_bic_values[is.finite(normal_bic_values)]

  if ( length(tumor_bic_values) == 0 || length(normal_bic_values) == 0) {
    return(empty_result)
  }

  # ============================================================
  # Tumor BIC used for comparison
  #
  # Median across matched resamples makes the tumor estimate
  # less dependent on a particular random sample.
  # ============================================================

  tumor_bic_gain <- stats::median( tumor_bic_values, na.rm = TRUE)

  # ============================================================
  # Empirical normal BIC cutoff
  #
  # normal_bic_quantile = 1:
  # use the maximum normal BIC gain.
  # ============================================================

  normal_bic_cutoff <- as.numeric(
    stats::quantile(
      normal_bic_values,
      probs = normal_bic_quantile,
      na.rm = TRUE,
      names = FALSE
    )
  )
  normal_bic_median <-
    stats::median(
      normal_bic_values,
      na.rm = TRUE
    )
  bic_excess <- tumor_bic_gain -  normal_bic_cutoff
  # ============================================================
  # Other tumor peak-shape criteria
  # ============================================================

  min_weight <- min( tumor_full_fit$weight1, tumor_full_fit$weight2)

  # ============================================================
  # Final two-peak decision
  #
  # The tumor must:
  # 1. Have greater two-peak BIC evidence than the normal empirical cutoff.
  # 2. Have two substantial components. min peak weight should be over 0.25.
  # 3. Have sufficiently separated peaks. sepatration_ad should be 2 fold of peak variance ??
  # 4. Have peaks centered around the expected heterozygous BAF region. (excluded)
  # 5. Have one peak below and one peak above 0.5.
  # ============================================================

  two_peak_supported <- is.finite(tumor_bic_gain) &&
    is.finite(normal_bic_cutoff) &&
    tumor_bic_gain >  normal_bic_cutoff + bic_margin &&
    min_weight >= min_component_weight &&
    tumor_full_fit$separation_sd >= min_separation_sd &&
    tumor_full_fit$peak1 < 0.5 &&
    tumor_full_fit$peak2 > 0.5



  if (!two_peak_supported) {

    return(
      c(
        "baf_maf" = 0.5,
        "baf_G" = 1,
        "baf_peak1" = tumor_full_fit$mu0,
        "baf_peak2" = NA_real_,
        "baf_weight1" = 1,
        "baf_weight2" = 0,
        "baf_center" = tumor_full_fit$mu0,
        "baf_d" = 0,
        "baf_sigma" = tumor_full_fit$sigma0,
        "baf_separation_sd" = 0,
        "baf_bic_gain" = tumor_bic_gain,
        "baf_bic_gain_full" = tumor_bic_gain_full,
        "normal_bic_cutoff" = normal_bic_cutoff,
        "normal_bic_median" =normal_bic_median,
        "baf_bic_excess" = bic_excess,
        "baf_n_tumor" =n_tumor,
        "baf_n_normal" = n_normal,
        "baf_n_match" = n_match
      )
    )
  }else{
    baf_maf <- 0.5 - tumor_full_fit$d
    baf_maf <- max( 0,min( baf_maf,0.5 ))

    return(
      c(
        "baf_maf" = baf_maf,
        "baf_G" = 2,
        "baf_peak1" = tumor_full_fit$peak1,
        "baf_peak2" = tumor_full_fit$peak2,
        "baf_weight1" = tumor_full_fit$weight1,
        "baf_weight2" = tumor_full_fit$weight2,
        "baf_center" = tumor_full_fit$center,
        "baf_d" = tumor_full_fit$d,
        "baf_sigma" = tumor_full_fit$sigma,
        "baf_separation_sd" = tumor_full_fit$separation_sd,
        "baf_bic_gain" = tumor_bic_gain,
        "baf_bic_gain_full" = tumor_bic_gain_full,
        "normal_bic_cutoff" = normal_bic_cutoff,
        "normal_bic_median" = normal_bic_median,
        "baf_bic_excess" = bic_excess,
        "baf_n_tumor" = n_tumor,
        "baf_n_normal" = n_normal,
        "baf_n_match" =n_match
      )
    )

  }

}

#' Estimate MAF Cluster Mean and Weight Using Gaussian Mixture Modeling
#'
#' Fits a Gaussian mixture model (GMM) to Minor allele frequency (MAF) values and returns the most prominent cluster's mean and weight.
#' Clusters with similar means (within a threshold) are merged. If there are insufficient values or nearly no variance, a simple mean is returned.
#'
#' @param maf_values Numeric vector of MAF values.
#'
#' @return Named numeric vector with elements:
#'   \item{gmm_mean}{The mean of the most prominent cluster (capped at 0.5).}
#'   \item{gmm_weight}{The mixture weight of the most prominent cluster.}
#'   \item{gmm_G}{Number of clusters after merging adjacent means.}
#'
#' @details
#' Uses \code{\link[mclust]{Mclust}} for Gaussian mixture modeling.
#'
#' @import mclust
#' @importFrom mclust Mclust mclustBIC
#' @export
EstimateMAFbyGMM <- function(maf_values){
  logit  <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  ilogit <- function(x) plogis(x)
  n_mafs <- length(maf_values)
  sd <- sd( maf_values, na.rm = T )
  if( n_mafs >= 3 && sd >= 1e-6 ){
    maf_values_trans <- logit(maf_values * 2)
    gmm_model <- mclust::Mclust(maf_values_trans, G = 1:5)
    means <- ilogit(gmm_model$parameters$mean) / 2
    weights <- gmm_model$parameters$pro
    if(gmm_model$G >= 2){
      re <- ClusterAdjacent(values = means, weights = weights, threshold = 0.01)
      means <- re$means
      weights <- re$weights
      num_components <- length(re$means)
    }else{
      num_components <- gmm_model$G
    }

    ## cluster the means if the difference is small


    sorted_indices <- order(-weights)
    gmm_mean <- means[sorted_indices][1] %>% as.numeric()
    #gmm_mean <- ifelse( gmm_mean > 0.48, 0.5, gmm_mean )
    gmm_mean <- min( gmm_mean, 0.5 )
    # gmm_variance <- variances[sorted_indices][1]
    gmm_weight <- weights[sorted_indices][1]

    Estimated_maf <- c(
      "gmm_mean" = gmm_mean,
      #"gmm_variance" = gmm_variance,
      "gmm_weight" = gmm_weight,
      "gmm_G" = num_components
    )}else{
      gmm_mean <- mean(maf_values,na.rm = T)
      #gmm_mean <- ifelse( gmm_mean > 0.48, 0.5, gmm_mean )
      gmm_mean <- min( gmm_mean, 0.5 )
      Estimated_maf <- c(
        "gmm_mean" =  gmm_mean,
        #"gmm_variance" = sd,
        "gmm_weight" = 1,
        "gmm_G" = 0)

    }

  return(Estimated_maf)

}


#' Merge Adjacent AI Segments Based on Similarity Criteria
#'
#' Iteratively merges adjacent rows in a data frame of AI segments if they meet criteria defined by \code{MergeAICheck}.
#' For merged segments, MAF values are re-estimated and segment information is updated.
#'
#' @param df A data frame or tibble of AI segments, with columns such as Chromosome, bin, Start, End, snp_count, each_maf, etc.
#' @param mergeai Numeric. Threshold for the difference in MAF (gmm_mean) between adjacent segments to allow merging.
#' @param snpmin Numeric. Minimum SNP count required for a segment to be considered as a separate segment.
#' @param tmp_maf A data frame or tibble containing maf values and genomic coordinates (must include Chromosome, Start, End, maf).
#'
#' @return A data frame or tibble with merged AI segments, updated maf estimates, and segment information.
#'
#' @details
#' This function uses \code{\link{MergeAICheck}} to determine if two adjacent segments should be merged.
#' For merged segments, maf values are re-estimated using \code{\link{EstimateMAFbyGMM}}.
#'
#' @importFrom dplyr filter arrange
#' @importFrom tibble tibble
#' @importFrom tidyr unnest_wider
#' @export
MergeAIRow <- function(df, mergeai, snpmin, tmp_maf  ) {

  if(nrow(df) > 1 ){
    i <- 1
    while ( i < ( nrow(df) ) ) {
      cur_row <- df[i,]
      next_row <- df[ i+1,]

      if (  MergeAICheck(cur_row = cur_row, next_row = next_row, mergeai = mergeai, snpmin = snpmin) ) {
        #print( paste0("i is :", i ))
        #print(paste0("Cur_row is :", cur_row))
        #print(paste0("Next_row is:", next_row))
        #cur_row_mafs = str_split(cur_row$each_maf, pattern = ";", simplify = TRUE) %>% as.numeric() %>% na.omit()
        #next_row_mafs = str_split(next_row$each_maf, pattern = ";", simplify = TRUE ) %>% as.numeric() %>% na.omit()
        cur_next_mafs <- tmp_maf %>%
          dplyr::filter( Chromosome == cur_row$Chromosome) %>%
          dplyr::filter( Pos >= cur_row$Start & Pos <= next_row$End)

        new_df <- tibble::tibble(
          Chromosome = cur_row$Chromosome,
          bin = paste(c(cur_row$bin, next_row$bin),collapse = ";"),
          Start = cur_row$Start,
          End = next_row$End,
          size = End - Start,
          snp_count = cur_row$snp_count + next_row$snp_count,
          #Estimated_maf= list(EstimateMAFbyGMM(maf_values = c(cur_row_mafs, next_row_mafs))) ,
          Estimated_maf = list(EstimateMAFbyGMM(maf_values = cur_next_mafs$maf )),
          nonzero_count = cur_row$nonzero_count + next_row$nonzero_count,
          each_maf = paste( cur_next_mafs$maf , collapse = ";")
        ) %>%
          tidyr::unnest_wider(col = Estimated_maf)
        df <- rbind(new_df, df[-c(i,i+1),])
        df <- df %>% arrange(by=Chromosome,Start,na.last = T)
        i <- 1
      }else{
        i <- i+1}
    }
  }
  return(df)
}


#' Merge Segments or AI Rows by Chromosome
#'
#' This function splits the input data by chromosome and merges segments or AI rows within each chromosome using either \code{\link{MergeAIRow}} or \code{\link{MergeSegRow}}.
#'
#' @param data A data frame or tibble containing segment or AI information. Must include a \code{Chromosome} column.
#' @param AIorSeg Character string, either \code{"AI"} to merge AI rows or \code{"Seg"} to merge segments.
#' @param mergeai Numeric. Threshold for the difference in MAF (gmm_mean) between adjacent segments to allow merging.
#' @param snpmin Numeric. Minimum SNP count required for a segment to be considered as a separate segment.
#' @param tmp_maf A data frame for use with \code{MergeAIRow}. Not used if \code{AIorSeg == "Seg"}.
#' @param merge_cov Numeric. Threshold for the difference in segment mean (gmm_mean) between adjacent segments to allow merging.
#'
#' @return A data frame or tibble with merged segments or AI rows for each chromosome.
#'
#' @details
#' Uses \code{\link{MergeAIRow}} if \code{AIorSeg == "AI"}, otherwise \code{\link{MergeSegRow}}.
#'
#' @importFrom dplyr arrange
#' @export
CallMerge <- function(data, AIorSeg, tmp_maf, snpmin, mergeai, mergecov){

  if(AIorSeg == "AI"){
    data<- data %>%
      dplyr::arrange(by=Chromosome, Start, na.last = T)
    bychr <- split(data,f=data$Chromosome)
    bychr <- lapply(bychr,function(x) MergeAIRow(df = x, tmp_maf = tmp_maf, snpmin = snpmin, mergeai = mergeai))
    }
  if(AIorSeg == "Seg"){
    data<- data %>%
      dplyr::arrange(by=Chromosome, Start, na.last = T)
    bychr <- split(data,f=data$Chromosome)
    bychr <- lapply(bychr,function(x) MergeSegRow(df = x, mergecov = mergecov ))
  }
  processed_df <- do.call(rbind,bychr)
  return(processed_df)
}


#' Bin MAF Data and Estimate MAF Metrics per Bin
#'
#' Bins input MAF data into fixed-size windows or fixed-snp-count window and calculates summary metrics (e.g., estimated MAF, nonzero counts) for each bin.
#'
#' @param data A data frame or tibble containing at least \code{Chromosome}, \code{Pos}, \code{Start}, \code{End} and \code{maf} columns.
#' @param datatype If individual tumor data then choose "tumor", if aggregated panel of normal samples choose "pon".
#' @param maxgap Maximum gap size inside a bin. If exceed then start another bin.( default: 2000000)
#' @param snpnum SNP number in each bin.( default: 20 )
#' @param maxbinsize Maximum bin size.( default: 1000000 )
#' @param minbinsize Minimum bin size.( default: 500000 )
#' @param minsnpcov Minimum coverage of SNP sites to be included. ( default: 20 )
#'
#' @return A tibble with one row per bin and columns:
#'   \item{Chromosome}{Chromosome identifier.}
#'   \item{bin}{Bin start coordinate.}
#'   \item{Start}{Minimum \code{End} value in the bin.}
#'   \item{End}{Maximum \code{End} value in the bin.}
#'   \item{nonzero_count}{Number of nonzero maf values in the bin.}
#'   \item{each_maf}{Semicolon-separated string of nonzero maf values in the bin.}
#'   \item{gmm_mean, gmm_weight, gmm_G}{Unnested maf or BAF metrics.}
#'
#' @importFrom dplyr mutate group_by summarise ungroup filter n
#' @importFrom tidyr unnest_wider
#' @export
BinMaf <- function(data, datatype,
                   maxgap = 2000000, snpnum = 20,
                   maxbinsize = 1000000, minbinsize = 500000,
                   minsnpcov = 20 ){

  if (datatype == "pon"){

    data_tmp <- data %>%
      dplyr::filter(alt_count != 0 & ref_count != 0) %>%
      dplyr::filter(alt_count + ref_count >= minsnpcov) %>%
      dplyr::group_by(Chromosome, Pos) %>%
      dplyr::ungroup()

    choosesample <- data_tmp %>%
      dplyr::group_by(sampleID) %>%
      dplyr::summarise( snp_count = dplyr::n(),
                 median_cov = median( alt_count + ref_count, na.rm = T)) %>%
      dplyr::arrange( -median_cov, -snp_count)

    data_tmp <- data_tmp %>%
      dplyr::filter( sampleID == choosesample$sampleID[1]) %>%
      dplyr::select( -sampleID)

  }else{ data_tmp <- data}

  ## Generate AI bin boundary files based on SNP counts
  boundary <- data_tmp %>%
      dplyr::group_by(Chromosome) %>%
      dplyr::arrange(Pos, .by_group = TRUE) %>%
      dplyr::mutate(
        gap = dplyr::lead(Pos) - Pos,
        gap = ifelse(is.na(gap), 0, gap)
      ) %>%
      dplyr::group_modify(~ {
        df <- .x
        n <- nrow(df)
        if (n == 0) return(df)

        current_bin_id <- 1
        start_index <- 1
        bin <- integer(n)

        for (i in 1:n) {
          bin[i] <- current_bin_id
          bin_count    <- i - start_index + 1
          current_span <- df$Pos[i] - df$Pos[start_index]

          # look-ahead
          if (i < n) {
            next_count <- bin_count + 1
            next_span  <- df$Pos[i + 1] - df$Pos[start_index]
            next_gap   <- df$Pos[i + 1] - df$Pos[i]
          } else {
            next_count <- Inf; next_span <- Inf; next_gap <- 0
          }

          # close BEFORE overshoot; always respect maxbinsize
          if ( current_span >= maxbinsize ||
               (current_span >= minbinsize &&
                (next_count > snpnum || next_span > maxbinsize || next_gap > maxgap)) ) {

            current_bin_id <- current_bin_id + 1
            start_index <- i + 1
          }
        }

        df$bin <- bin
        df
      }) %>%
      dplyr::ungroup()

    ## summarize boundaries
    boundary <- boundary  %>%
      dplyr::group_by(Chromosome, bin) %>%
      dplyr::summarise(
        Start = min(Pos),
        End = max(Pos),
        size = max(Pos) - min(Pos),
        snp_count = dplyr::n())

    ## overlap data with boundary to get snps in each bin
    data <- data %>%
      dplyr::filter( alt_count >0 & ref_count >0)
    data.table::setDT(data)
    data.table::setDT(boundary)
    data$Pos_end <- data$Pos
    data.table::setkey( data , Chromosome, Pos, Pos_end )
    data.table::setkey( boundary, Chromosome, Start, End )
    binned <- data.table::foverlaps(
      data, boundary,
      by.x = c("Chromosome", "Pos","Pos_end"),
      by.y = c("Chromosome","Start","End"),
      type = "within", nomatch = 0
    )

  if( datatype == "tumor"){
    binned <- binned %>%
      dplyr::group_by(Chromosome, bin, Start, End, size, snp_count) %>%
      dplyr::summarise( Estimated_maf = list(EstimateMAFbyGMM(maf_values = maf[maf != 0] )),
      nonzero_count = sum(maf != 0),
      each_maf = paste(maf[maf != 0 ], collapse = ";"),
      .groups = "drop"
      ) %>%
  tidyr::unnest_wider(col = "Estimated_maf")%>%
  dplyr::filter( !is.na( gmm_mean))
  }else{
    ## summarise snp sites in each bin for each sample
    binned <- binned %>%
      dplyr::mutate( bin = paste0("bin_",Chromosome,"_", bin)) %>%
      dplyr::group_by(sampleID, Chromosome, bin, Start, End, size) %>%
      dplyr::summarise(
        snp_count = dplyr::n(),
        psb_snp_sum_depth = sum(alt_count) + sum(ref_count),
        psb_snp_sum_alt = sum(alt_count),
        psb_snp_baf = psb_snp_sum_alt/psb_snp_sum_depth,
        psb_snp_median_depth = median( alt_count + ref_count, na.rm = TRUE),
        psb_snp_median_baf = median( alt_count/(alt_count + ref_count)),
        psb_snp_median_maf = median( maf, na.rm = TRUE),
        psb_snp_mafs = paste( maf, collapse = ","),
        psb_snp_bafs = paste( baf, collapse = ","),
        .groups = "drop"
      )
  }

  return(binned)

}



#' Refine Breakpoints Within a Segment Using MAF Data
#'
#' Refines breakpoints within a segment using minor allele frequency (MAF) data. If enough informative MAF sites are present, the segment is binned and can be split into finer regions using either stepwise merging or CBS (circular binary segmentation). Optionally, PON-based bias correction is applied to the resulting segments.
#'
#' @param seg_row Data frame row (list or tibble row) representing a single segment. Must have columns: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean, Segment_Mean_raw, Count, Baseline_cov, gatk_gender, pipeline_gender, size.
#' @param maf Data frame or tibble containing MAF data. Must include columns: Chromosome, Pos, maf.
#' @param mergeai Numeric. Threshold for the difference in MAF (gmm_mean) between adjacent segments to allow merging under \code{"merge"} mode segmentation.
#' @param snpmin Numeric. Minimum SNP count required for a segment to be considered as a separate segment under \code{"merge"} mode segmentation.
#' @param maxgap Numeric. Maximum allowed gap between SNPs within a bin.
#' @param snpnum Integer. Target number of SNPs per bin.
#' @param maxbinsize Numeric. Maximum allowed bin size (bp).
#' @param minbinsize Numeric. Minimum allowed bin size (bp). The minimum segment size under \code{"merge"} mode is 2*minbinsize.
#' @param minsnpcov Integer. Minimum coverage of SNP sites to be included.
#' @param segmethod Character. Segmentation method to use: if \code{"merge"}, perform stepwise merging; if \code{"cbs"}, perform CBS (circular binary segmentation).
#' @param cbssmooth Character. If using the \code{"cbs"} segmentation method, set to \code{"yes"} to apply smoothing before segmentation, or \code{"no"} to skip smoothing.
#' @param pon_ref Data frame. Panel of normal reference for bias correction (required for bias correction step).
#' @param gender Character. If \code{"female"}, the X chromosome will also be proceed.
#'
#' @return A data frame with the refined segment(s), including updated breakpoints, MAF metrics, and a \code{BreakpointSource} column indicating whether breakpoints were post-processed or from GATK.
#'
#' @details
#' The function first bins the MAF data within the segment. If \code{segmethod = "merge"}, segments are merged stepwise based on the MAF difference and SNP count. If \code{segmethod = "cbs"}, CBS segmentation is performed on the binned MAF values, with optional smoothing. After segmentation, bias correction using the panel of normal can be applied. The function returns refined segments with updated metrics and a \code{BreakpointSource} label.
#'
#' @importFrom dplyr filter select mutate group_by summarise row_number relocate
#' @importFrom tidyr unnest_wider
#' @importFrom DNAcopy CNA smooth.CNA segment
#' @export
SearchBreakpoint <- function(seg_row, maf, pon_ref, gender, out_dir, prefix,
                             mergeai = 0.05,
                             snpmin = 3,
                             maxgap = 1000000, snpnum = 20,
                             maxbinsize = 1000000,minbinsize = 500000,
                             minsnpcov = 20,
                             segmethod = "cbs",
                             cbssmooth = "no"){
  mergecbs <- function(df) {

    if(nrow(df) > 1 ){
      i <- 1
      while ( i < ( nrow(df) ) ) {

        cur_row <- df[i,]
        next_row <- df[ i+1,]

        if (  abs(cur_row$seg.mean - next_row$seg.mean) <= 0.02 ) {
          new_df <- data.frame(
            startRow = cur_row$startRow,
            endRow = next_row$endRow,
            chrom = cur_row$chrom,
            seg.mean = ifelse( cur_row$n_bins > next_row$n_bins,
                               cur_row$seg.mean,  next_row$seg.mean),
            n_bins = next_row$endRow - cur_row$startRow + 1
          )
          df <- rbind(new_df, df[-c(i,i+1),])
          df <- df %>% dplyr::arrange(by=startRow,na.last = T)
          i <- 1
        }else{i <- i+1}
      }
    }
    return(df)
  }

  tmp_seg <- data.frame(
    Sample = seg_row$Sample,
    Chromosome = seg_row$Chromosome,
    Start = seg_row$Start,
    End = seg_row$End,
    Num_Probes = seg_row$Num_Probes,
    Segment_Mean = seg_row$Segment_Mean,
    gatk_SM_raw = seg_row$Segment_Mean_raw,
    gatk_count = seg_row$Count,
    gatk_baselinecov = seg_row$Baseline_cov,
    gatk_gender = seg_row$gatk_gender,
    pipeline_gender = seg_row$pipeline_gender,
    cov_mad = seg_row$MAD,
    MAF = NA,
    MAF_Probes = 0,
    MAF_gmm_G = NA,
    MAF_gmm_weight = NA,
    #MAF_gmm_variance = NA,
    size = seg_row$size)
if( seg_row$Chromosome != "Y" || ( seg_row$Chromosome != "X" && gender != "male" ) ) {
  tmp_maf <- maf %>%
    dplyr::filter( Chromosome == seg_row$Chromosome & Pos >= seg_row$Start & Pos <= seg_row$End ) %>%
    dplyr::filter( maf != 0 ) %>%
    dplyr::mutate( BAF = alt_count/(alt_count + ref_count))
  n_maf <- nrow(tmp_maf)
  if( is.null(n_maf) ){ n_maf <- 0}
  if( n_maf > snpmin ){
    binned_data <- BinMaf(data = tmp_maf,
                          datatype = "tumor",
                          maxgap = maxgap,
                          snpnum = snpnum,
                          maxbinsize = maxbinsize,
                          minbinsize = minbinsize,
                          minsnpcov = minsnpcov)
    max_nonzero_n <- max(binned_data$nonzero_count,na.rm = T)


    if(segmethod == "merge"){

      for ( min_prob in c(0: min( snpmin, max_nonzero_n - 1 ) ) ){
        binned_data <- binned_data %>% dplyr::filter(nonzero_count >= min_prob )
        binned_data <- CallMerge(data = binned_data, AIorSeg = "AI", mergeai = mergeai, snpmin = snpmin, tmp_maf = tmp_maf)
      }
      binned_data <- binned_data %>% dplyr::filter( End - Start >= 2*minbinsize)
      if( nrow(binned_data) > 0 ){
        merge_ai <- CallMerge(data = binned_data, AIorSeg = "AI", mergeai= mergeai, snpmin = snpmin, tmp_maf = tmp_maf)
        merge_ai <- merge_ai %>%
          dplyr::select(-bin) %>%
          dplyr::mutate(size = End - Start)
        tmp_seg <- data.frame(
          Sample = seg_row$Sample,
          Chromosome = seg_row$Chromosome,
          Start = merge_ai$Start,
          End = merge_ai$End,
          Num_Probes = round(seg_row$Num_Probes * as.numeric(merge_ai$size/sum(merge_ai$size)),digits = 0),
          Segment_Mean = seg_row$Segment_Mean,
          gatk_SM_raw = seg_row$Segment_Mean_raw,
          gatk_count = seg_row$Count,
          gatk_baselinecov = seg_row$Baseline_cov,
          gatk_gender = seg_row$gatk_gender,
          pipeline_gender = seg_row$pipeline_gender,
          cov_mad = seg_row$MAD,
          MAF = merge_ai$gmm_mean,
          MAF_Probes = merge_ai$nonzero_count,
          MAF_gmm_G = merge_ai$gmm_G,
          MAF_gmm_weight = merge_ai$gmm_weight,
          #MAF_gmm_variance = merge_ai$gmm_variance,
          size = merge_ai$End - merge_ai$Start
        )

        tmp_seg[1,"Start"] <- min( seg_row$Start, merge_ai$Start ,na.rm = T)
        tmp_seg[nrow(tmp_seg),"End"] <- max( seg_row$End, merge_ai$End, na.rm = T)

      }
    }


    if(segmethod == "cbs"){
      # Prepare data for CBS (use gmm_mean as the segmentation track)
      maf_seg_data <- data.frame(
        chrom = binned_data$Chromosome,
        maploc = binned_data$Start,    # or bin center if you prefer
        maf = binned_data$gmm_mean,
        weights = ifelse(
          binned_data$nonzero_count >= 2 * snpmin, 1, 0.5
          )
        )

      # Create CNA object
      maf_CNA <- DNAcopy::CNA(
        genomdat = maf_seg_data$maf,
        chrom = maf_seg_data$chrom,
        maploc = maf_seg_data$maploc,
        data.type = "logratio"
      )

        maf_CNA_smoothed <- DNAcopy::smooth.CNA(maf_CNA, smooth.region = 5)
        maf_cbs <- DNAcopy::segment(maf_CNA_smoothed,weights = maf_seg_data$weights, min.width = 3,verbose = 1)

        # Combine segRows with segmentation information
        seg_info <- maf_cbs$segRows %>%
          dplyr::mutate(
            chrom = maf_cbs$output$chrom,
            seg.mean = maf_cbs$output$seg.mean,
            n_bins = endRow - startRow + 1
          )
        last_bin <- nrow(maf_seg_data)
         seg_info <- seg_info %>%
           dplyr::filter(n_bins > 3 )
         seg_info <-   mergecbs(df = seg_info)
         ## double check firs and last bin

         seg_info[1,1] <- 1
         seg_info[nrow(seg_info),"endRow"] <- last_bin

         seg_info$n_bins <- seg_info$endRow - seg_info$startRow + 1
         seg_info <- seg_info %>%
           dplyr::mutate( seg_id = dplyr::row_number() )


        # Assign merged segment index to each bin
        binned_data$merge_index <- rep(
          seg_info$seg_id,
          times = seg_info$n_bins
        )
      tmp_seg <- binned_data %>%
        dplyr::group_by( Chromosome, merge_index) %>%
        dplyr::summarise(
          Start = min( Start),
          End = max(End),
          size = max(End) - min( Start),
          Estimated_maf = list(EstimateMAFbyGMM(maf_values = as.numeric( unlist(strsplit(each_maf, ";")) ))),
          nonzero_count = sum(nonzero_count),
          each_maf = paste(  as.numeric( unlist(strsplit(each_maf,";")) ) , collapse = ";")
        ) %>%
        tidyr::unnest_wider(col = Estimated_maf) %>%
        dplyr::mutate( Sample = seg_row$Sample,
                Num_Probes = round(seg_row$Num_Probes * as.numeric(size/sum(size)) , digits =  0 ),
                Segment_Mean = seg_row$Segment_Mean,
                gatk_SM_raw = seg_row$Segment_Mean_raw,
                gatk_count = seg_row$Count,
                gatk_baselinecov = seg_row$Baseline_cov,
                gatk_gender = seg_row$gatk_gender,
                pipeline_gender = seg_row$pipeline_gender,
                cov_mad = seg_row$MAD,
                MAF = gmm_mean,
                MAF_Probes = nonzero_count,
                MAF_gmm_G = gmm_G,
                MAF_gmm_weight = gmm_weight
        ) %>%
        dplyr::select( -gmm_mean, -nonzero_count, -gmm_G, -gmm_weight, -merge_index, -each_maf) %>%
        dplyr::relocate( Sample, Chromosome, Start, End, Num_Probes, Segment_Mean, gatk_SM_raw,
                  gatk_count, gatk_baselinecov, gatk_gender, pipeline_gender, cov_mad,
                  MAF, MAF_Probes, MAF_gmm_G, MAF_gmm_weight, size )

    }


  tmp_seg<- RefineMAF( tmp_seg, pon_ref, tmp_maf, out_dir, prefix )
    }else{ tmp_seg$gmm_mean_corr <- NA
          tmp_seg$balance_tag <- NA
  }
  if(nrow(tmp_seg) > 1){
    tmp_seg$BreakpointSource <- "Postprocess"}else{tmp_seg$BreakpointSource <- "GATK"}
}else{
    tmp_seg$gmm_mean_corr <- NA
    tmp_seg$balance_tag <- NA
 tmp_seg["BreakpointSource"] <- "GATK"
}
  return(tmp_seg)
}


#' Assign Quality Tag to a Segment Based on maf GMM Metrics and coverage variance
#'
#' Assigns a quality tag ("PASS" or "FAILED") to a segment based on the maf GMM weight and probe count thresholds.
#'
#'A segment is marked as "PASS" if the maf GMM weight is at least 0.5, the probe count is at least twice the minimum SNP count (\code{snpmin}), coverage MAD is less than 2.
#'Otherwise marked as "FAILED". If the maf GMM weight is 0.5 to 0.8, the balance tag is balanced, but the MAF is < 0.48, it is also FAILED.
#' @param Chromosome Character.
#' @param MAF Numeric. Estimated MAF value.
#' @param MAF_gmm_weight Numeric. The mixture weight of the most prominent maf GMM cluster.
#' @param MAF_Probes Integer. The number of nonzero maf probes in the segment.
#' @param MAF_gmm_G Integer. Number of maf GMM clusters.
#' @param snpmin Integer. The minimum nonzero SNP count.
#' @param balance_tag Character. Balanced test result of BAF of the segment.
#' @param sampletype Character. ffpe or ff.
#' @param cov_mad Numeric. MAD of coverage.
#' @return Character. "PASS" if the segment passes quality checks, "FAILED" otherwise.
#'
#' @examples
#' AddQualTag(Chromosome = "1", MAF_gmm_weight = 0.4, MAF_Probes = 20, MAF_gmm_G = 2, snpmin = 7, balance_tag = "balanced", cov_mad = 1, sampletype="ffpe")
#' AddQualTag(Chromosome = "X", MAF_gmm_weight = 0.2, MAF_Probes = 10, MAF_gmm_G = 1, snpmin = 7, balance_tag = "balanced", cov_mad = 1, sampletype="ff")
#'
#' @export
AddQualTag <- function(Chromosome, MAF, MAF_gmm_weight, MAF_Probes, MAF_gmm_G, snpmin, balance_tag, cov_mad, sampletype){

  ## MAF quality tag
  balance_tag <- tolower(trimws(as.character(balance_tag)))
  sampletype <- tolower(trimws(as.character(sampletype)))
  MAF <- as.numeric(MAF)

  maf_tag <- "FAILED"

  if (
    !is.na(MAF_gmm_G) &&
    !is.na(MAF_Probes) &&
    !is.na(snpmin) &&
    !is.na(MAF_gmm_weight) &&
    !is.na(balance_tag) &&
    !is.na(MAF)
  ) {

    if (MAF_Probes >= 2 * snpmin && MAF_gmm_weight >= 0.3) {

      if (balance_tag == "balanced") {

        if (MAF >= 0.48) {
          maf_tag <- "PASS"
        }

      } else {

        if ( abs(MAF - 0.5) > 0.00001 ) {
          maf_tag <- "PASS"
        }
      }
    }
  }


  if(FALSE){
    if(!is.na(MAF_gmm_G)){
      if( MAF_Probes >= 2*snpmin ){
        if( MAF_gmm_weight >= 0.3 ){
          if( balance_tag == "balanced" ){
            if( MAF < 0.48 ){ maf_tag <- "FAILED"
            }else{
              maf_tag <- "PASS"
            }
          }else{
            if(MAF == 0.5){
              maf_tag <- "FAILED"
            }else{
              maf_tag <- "PASS"
            }
          }
        }else{ maf_tag <- "FAILED"}
      }else{ maf_tag <- "FAILED" }
    }else{ maf_tag <- "FAILED"}
  }



  ## coverage quality tag
  if(sampletype == "ff"){
    if(cov_mad < 1.5){ cov_tag <- "PASS"}else{
      cov_tag <- "FAILED"
    }
    }

high_gc_chrom <- c("9","16","17","21","22", "19","Y")
  if(sampletype == "ffpe"){
    if( ! Chromosome %in% high_gc_chrom ){
      if( cov_mad <= 3){ cov_tag <- "PASS" }else{ cov_tag <- "FAILED" }
      }else if( Chromosome %in% c("9","16","17","21","22") ){
        if( cov_mad <= 4 ){ cov_tag <- "PASS"}else{ cov_tag <- "FAILED" }
      }else if( Chromosome %in% c("19","Y") ) {
        if( cov_mad <= 6 ){ cov_tag <- "PASS"}else{ cov_tag <- "FAILED" }
      }
    }


  if( maf_tag == "PASS" && cov_tag == "PASS" ) {
  Qualtag <- "PASS"
}else{
    Qualtag <- "FAILED" }


 if(Chromosome %in% c("X","Y")){
   Qualtag <- "EXCLUDE"
 }
  return(Qualtag)
}

#' Plot Tumor and Normal BAF Density
#'
#' Generates tumor and normal BAF density plots for each segment with at least
#' 10 valid BAF values in both samples. Tumor BAF is shown in red and normal
#' BAF in blue, with gray BAF guide lines.
#'
#' @param seg_corr Data frame containing segment coordinates and BAF values.
#' @param out_dir Character. Output directory.
#' @param prefix Character. Sample ID used in output filenames.
#'
#' @return NULL. PNG files are saved to \code{BAF_screenshots}.
#'
#' @export
PlotBAF <- function( seg_corr, out_dir, prefix ){

  baf_out_dir <- file.path(out_dir, "BAF_screenshots")

  if (!dir.exists(baf_out_dir)) {
    dir.create(
      baf_out_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }


  for ( x in c(1:nrow(seg_corr)) ) {
    chr <- as.character(seg_corr$Chromosome[x])
    seg_start <- as.numeric(seg_corr$Start[x])
    seg_end <- as.numeric(seg_corr$End[x])
    tumor_baf <- as.numeric(unlist(strsplit(seg_corr$each_BAFs[x], ";")))
    normal_baf <- as.numeric(unlist(strsplit(seg_corr$pon_bafs[x], ",")))
    # Tumor BAF
    tumor_baf <- tumor_baf[
      is.finite(tumor_baf) &
        !is.na(tumor_baf) &
        tumor_baf >= 0 &
        tumor_baf <= 1
    ]

    # Normal BAF
    normal_baf <- normal_baf[
      is.finite(normal_baf) &
        !is.na(normal_baf) &
        normal_baf >= 0 &
        normal_baf <= 1
    ]

    n_tumor_baf <- length(tumor_baf)
    n_normal_baf <- length(normal_baf)

    target_normal_n <- if (n_tumor_baf < 500) {
      500
    } else {
      n_tumor_baf
    }

    if (n_normal_baf > target_normal_n) {
      normal_baf <- sample(
        normal_baf,
        size = target_normal_n,
        replace = FALSE
      )
    }

    n_normal_baf <- length(normal_baf)


    if(n_tumor_baf >= 10 & n_normal_baf >= 10){
      out_file <- file.path(
        baf_out_dir,
        paste0(
          prefix, "_",
          chr, "_",
          format(seg_start, scientific = FALSE, trim = TRUE), "_",
          format(seg_end, scientific = FALSE, trim = TRUE),
          "_BAF_distribution.png"
        ))
      baf_guides <- seq(0.1, 1.0, by = 0.1)
      baf_guides <- baf_guides[baf_guides != 0.5]

      # ------------------------------------------------------------
      # Combine tumor and normal BAF
      # ------------------------------------------------------------

      baf_df <- data.frame(
        BAF = c(tumor_baf, normal_baf),
        Sample = c(
          rep("Tumor", length(tumor_baf)),
          rep("Normal", length(normal_baf))
        )
      )

      # ------------------------------------------------------------
      # Guide lines
      # ------------------------------------------------------------

      # Dashed guide lines:
      # 0.1, 0.2, 0.3, 0.4, 0.6, ..., 1.0
      # 0.5 is excluded because it will be a solid gray line
      baf_guides <- seq(0.1, 1.0, by = 0.1)
      baf_guides <- baf_guides[abs(baf_guides - 0.5) > 1e-8]

      # ------------------------------------------------------------
      # Density plot
      # ------------------------------------------------------------

      p <- ggplot2::ggplot(
        baf_df,
        aes(
          x = BAF,
          color = Sample
        )
      ) +

        # Dashed gray guide lines
        ggplot2::geom_vline(
          xintercept = baf_guides,
          linetype = "dashed",
          linewidth = 0.4,
          color = "gray70"
        ) +

        # BAF = 0.5 solid gray line
        ggplot2::geom_vline(
          xintercept = 0.5,
          linetype = "solid",
          linewidth = 0.7,
          color = "gray50"
        ) +

        # Tumor and normal density curves
        ggplot2::geom_density(
          linewidth = 1.2,
          adjust = 1,
          na.rm = TRUE
        ) +

        # Tumor = red; Normal = blue
        ggplot2::scale_color_manual(
          values = c(
            "Tumor" = "red",
            "Normal" = "blue"
          )
        ) +

        ggplot2::scale_x_continuous(
          limits = c(0, 1),
          breaks = seq(0, 1, by = 0.1),
          expand = c(0, 0)
        ) +

        ggplot2::labs(
          title = paste0("Tumor n=", n_tumor_baf, "    ", "Normal n=", n_normal_baf ),
          x = "BAF",
          y = "Density",
          color = NULL
        ) +

        ggplot2::theme_classic(base_size = 12) +
        ggplot2::theme(
          legend.position = "top",
          plot.title = ggplot2::element_text(
            size = 10,
            hjust = 0.5
          )
        )

      # ------------------------------------------------------------
      # Save PNG
      # ------------------------------------------------------------

      ggplot2::ggsave(
        filename = out_file,
        plot = p,
        device = "png",
        width = 5,
        height = 4,
        units = "in",
        dpi = 120
      )

    }
    }

  }


#' Refine MAF estiamtion results
#'
#' Adjusts the minor allele frequency (MAF) values in each segment
#'
#' @param tmp_seg Data frame or data.table. Segmented data to be corrected, must include columns: Chromosome, Start, End, Sample, Num_Probes, Segment_Mean, gatk_SM_raw, gatk_count, gatk_baselinecov, gatk_gender, pipeline_gender, MAF, MAF_Probes, MAF_gmm_G, MAF_gmm_weight, size.
#' @param pon_ref Data frame or data.table. Panel of normal reference, must include columns: Chromosome, Start, End, pon_bafs (comma-separated string of PoN BAF values).
#' @param tmp_maf Data frame or data.table. Per-SNP MAF data, must include columns: Pos, maf.
#'
#' @return A data frame with bias-corrected MAF values for each segment. The \code{MAF} column is updated with the corrected value, and columns \code{gmm_mean_corr}, \code{each_mafs}, and \code{pon_mafs} are removed.
#'
#' @details
#' For each segment, the function performs an interval join with the PoN reference to obtain the PoN MAF distribution. If the segment MAF is not significantly different from the PoN (by Wilcoxon test or threshold), applies a logit-based correction. Otherwise, the original segment MAF is retained. The function uses the panel median MAF for centering and clamps corrected values to [0, 0.5].
#'
#' @importFrom data.table setDT setkey foverlaps
#' @importFrom dplyr group_by summarise mutate select relocate filter
#' @export
RefineMAF <- function( tmp_seg, pon_ref, tmp_maf, out_dir, prefix ){
  #logit  <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  #ilogit <- function(x) plogis(x)
  extract_BAF <- function( tmp_maf, Start, End){

    each_BAFs <- tmp_maf %>%
      dplyr::filter( Start <= Pos & End >= Pos ) %>%
      dplyr::filter( maf != 0 )

    each_BAFs <- paste0( each_BAFs$BAF, collapse = ";")
    return(each_BAFs)
  }

  find_peaks <- function(BAFs){
    baf <- BAFs
    # Estimate density
    dens <- density(baf, bw = "nrd0")  # 'bw' can be tuned

    # Find local maxima (peaks)
    peaks_idx <- which(diff(sign(diff(dens$y))) == -2) + 1
    peak_positions <- dens$x[peaks_idx]
    peak_heights <- dens$y[peaks_idx]

    # Assign points to nearest peak (optional)
    assignments <- sapply(baf, function(x) peak_positions[which.min(abs(x - peak_positions))])
    proportions <- table(assignments) / length(baf)

    # Results
    re <- data.frame(
      peak_mean = peak_positions,
      proportion = as.numeric(proportions[match(peak_positions, names(proportions))])
    )
    return(re)
  }

  test_balance_KDE <- function( each_BAFs){
    each_BAFs <- each_BAFs[ !is.na(each_BAFs)]
    kde_re <- find_peaks(each_BAFs)
    kde_re <- kde_re %>% dplyr::filter( proportion >= 0.25)
    if(nrow(kde_re) == 1 ){
      balance <- "balanced"
    }else{
      kde_re <- kde_re %>% arrange( desc(proportion))
      peak_diff <- abs(kde_re$peak_mean[1] - kde_re$peak_mean[2])
      if( peak_diff <= 0.1 || max(kde_re$peak_mean) <= 0.52 ){
        balance <- "balanced"
      }else{
        balance <- "imbalanced"
      }

    }
    return(balance)
  }
  adjustmaf <- function( pon_bafs, tumor_bafs, gmm_mean){
    gmm_mean <- gmm_mean[1]
    pon_bafs <- as.numeric(unlist(strsplit(pon_bafs, ",")))
    tumor_bafs <- as.numeric(unlist(strsplit(tumor_bafs, ";")))

    ## BAF balance test by KDE if the peak diff is smaller than 0.1 then will be balanced otherwise imbalanced.
    if( !is.na(gmm_mean) ){
      if(gmm_mean >= 0.35 && gmm_mean < 0.5 ){
        balanced_round1 <- test_balance_KDE( tumor_bafs )
      }else if( gmm_mean < 0.35 ){
        balanced_round1 <- "imbalanced"
      }else if( gmm_mean == 0.5){
        balanced_round1 <- "balanced"
      }else{
        balanced_round1 <- NA }
    }else{ balanced_round1 <- NA }
    # ============================================================
    # NEW:
    # RAW BAF-BASED ESTIMATION
    #
    # Only run this additional estimator when the original
    # GMM estimate is in the difficult near-balanced region. gmm_mean >= 0.35
    # ============================================================

    if ( is.finite(gmm_mean) && gmm_mean > 0.35) {
      baf_estimate <-
        EstimateMAFfromBAF(
          tumor_bafs = tumor_bafs,
          normal_bafs = pon_bafs,
          baf_lower = 0.15,
          baf_upper = 0.85,
          min_n = 10L,
          n_resample = 100L,
          normal_bic_quantile = 0.95,
          bic_margin = 1,
          min_component_weight = 0.25,
          min_separation_sd = 2,
          center_lower = 0.42,
          center_upper = 0.58,
          seed = 1L
        )
    if (is.na(baf_estimate["baf_peak1"]) && is.na(baf_estimate["baf_peak"]) ){
      balanced_round2 <- NA
    } else if (!is.na(baf_estimate["baf_peak1"]) && is.na(baf_estimate["baf_peak"])) {
      balanced_round2 <- "balanced"
    } else {
      balanced_round2 <- "imbalanced"
    }

    if (is.na(balanced_round1) && is.na(balanced_round2)) {
      balance_tag <- NA_character_
    } else if (is.na(balanced_round1)) {
      balance_tag <- balanced_round2
    } else if (is.na(balanced_round2)) {
      balance_tag <- balanced_round1
    } else if (balanced_round1 == "balanced") {
      if (balanced_round2 == "imbalanced") {
        balance_tag <- "imbalanced"
      } else {
        balance_tag <- "balanced"
      }
    } else if (balanced_round1 == "imbalanced") {
      balance_tag <- "imbalanced"
    }

    return( list(gmm_mean_corr = baf_estimate["baf_maf"], balance_tag = balance_tag) ) } else{
      return( list(gmm_mean_corr = gmm_mean, balance_tag = "imbalanced") ) }

    }


  data.table::setDT(tmp_seg)
  data.table::setDT(pon_ref)
  data.table::setkey( pon_ref, Chromosome, Start, End )
  data.table::setkey( tmp_seg, Chromosome, Start, End )
  seg_corr <- data.table::foverlaps(
    pon_ref, tmp_seg,
    by.x = c("Chromosome", "Start","End"),
    by.y = c("Chromosome","Start","End"),
    type = "any", nomatch = 0
  ) %>%
    dplyr::group_by( Chromosome, Sample, Start, End, Num_Probes, Segment_Mean, gatk_SM_raw, gatk_count, gatk_baselinecov, gatk_gender,
              pipeline_gender, cov_mad, MAF, MAF_Probes, MAF_gmm_G, MAF_gmm_weight, size ) %>%
    dplyr::summarise( pon_mafs = paste( pon_mafs, collapse = ","),
                      pon_bafs = paste( pon_bafs, collapse = ","))


  if( nrow(seg_corr) > 0){
    seg_corr <- seg_corr %>%
      dplyr::rowwise() %>%
      dplyr::mutate( each_BAFs = extract_BAF( tmp_maf, Start, End))
    ## Plot BAF distribution

    PlotBAF( seg_corr, out_dir, prefix )
    ##
    seg_corr <- seg_corr %>%
      dplyr::mutate( gmm_mean_corr = list(adjustmaf(pon_bafs = pon_bafs, tumor_bafs = each_BAFs , gmm_mean=MAF))) %>%
      tidyr::unnest_wider( gmm_mean_corr) %>%
      dplyr::select( -pon_mafs, -each_BAFs, -pon_bafs)
  }

  return(seg_corr)
}


#' Estimate Beta-Binomial K Parameter for Each Segment
#'
#' Calculates the beta-binomial concentration parameter (\code{K}) for each segment based on local sequencing depth and a fitted dispersion model.
#'
#' @param result Data frame. Segment-level results, must include columns \code{Chromosome}, \code{Start}, and \code{End}.
#' @param theta_fit List or data frame. Fitted dispersion model, must include \code{breaks} (numeric vector of depth bin edges) and \code{theta_table} (data frame with columns \code{depth_bin} and \code{theta}).
#' @param maf Data frame. Allelic count data for all SNPs, must include columns \code{chr}, \code{Pos}, \code{ref_count}, and \code{alt_count}.
#'
#' @return The input \code{result} data frame with additional columns:
#'   \item{depth}{Median SNP depth for the segment.}
#'   \item{depth_bin}{Depth bin index for the segment.}
#'   \item{theta}{Estimated dispersion parameter for the segment.}
#'   \item{K}{Estimated beta-binomial concentration parameter for the segment.}
#'
#' @details
#' For each segment, the function extracts all SNPs within the segment boundaries from the \code{maf} data, calculates the median depth, assigns a depth bin according to \code{theta_fit$breaks}, joins the corresponding theta value from \code{theta_fit$theta_table}, and computes \code{K} using the formula \eqn{K = \frac{\text{depth}}{1 + \text{depth} \cdot \theta} - 1}.
#'
#' @importFrom dplyr filter mutate left_join rowwise
#' @export
EstimateK <- function( result, theta_fit, maf){
  extract_depth <- function( maf, Chromosome, Start, End){
    colnames(maf)[1] <- "chr"
    tmp_maf <- maf %>%
      dplyr::filter( chr == Chromosome & Pos <= End & Pos >= Start) %>%
      dplyr::mutate( depth = ref_count + alt_count)

    depth <- median(tmp_maf$depth)
    return(depth)
  }
  result <- result %>%
    dplyr::rowwise() %>%
    dplyr::mutate( depth = extract_depth( maf, Chromosome, Start, End)) %>%
    dplyr::mutate( depth_bin = cut(depth, breaks = theta_fit$breaks, include.lowest = TRUE, labels = FALSE) ) %>%
    dplyr::left_join(theta_fit$theta_table, by = c("depth_bin")) %>%
    dplyr::mutate( K = depth/(1 + (depth - 1)*theta) - 1 )
  return(result)

}



