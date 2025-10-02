#' Convert Segment Mean to Original Coverage
#'
#' Calculates the original coverage for a segment based on its segment mean (log2 ratio), gender, chromosome, and assumed diploid coverage.
#'
#' @param SM Numeric. The segment mean (log2 ratio) from CNV analysis.
#' @param gender Character. Sample gender, either \code{"male"} or \code{"female"}.
#' @param chromosome Character or factor. Chromosome name (e.g., "1", "2", ..., "X", "Y").
#' @param diploid_cov Numeric. The expected coverage for a diploid region (e.g., 100).
#'
#' @return Numeric. The estimated original coverage for the segment.
#'
#' @details
#' For males and sex chromosomes (\code{"X"} or \code{"Y"}), the diploid coverage is halved.
#' For all other cases, the segment mean is scaled by the full diploid coverage.
#'
#' @examples
#' SegmentMeanToOriCov(SM = 0, gender = "male", chromosome = "X", diploid_cov = 100)
#' SegmentMeanToOriCov(SM = 0.5, gender = "female", chromosome = "3", diploid_cov = 100)
#'
#' @export
SegmentMeanToOriCov <- function( SM, gender, chromosome, diploid_cov ){
  # According to the segment mean, calculate the original coverage of the segment.
  # Assume the diploid coverage is diploid_cov (e.g., 100).
  if( gender == "male" && chromosome %in% c("X", "Y")){
    SegCov <- 2^( SM )* (diploid_cov/2)

  }else{

    SegCov <- 2^( SM )* diploid_cov
  }

  return(SegCov)
}


#' Estimate Standard Deviation of a Numeric Vector
#'
#' Calculates the standard deviation (sigma) from the sample variance of a numeric vector.
#'
#' @param x Numeric vector.
#'
#' @return Numeric. The estimated standard deviation of \code{x}.
#'
#' @details
#' This function computes the sample variance using \code{var(x)} and returns its square root.
#'
#' @examples
#' EstimateVariance(c(1, 2, 3, 4, 5))
#'
#' @export
EstimateVariance <- function(x){
  # Estimate variance and standard deviation
  variance_x <- var(x)  # Variance
  sigma_x <- sqrt(variance_x)          # Standard deviation
  return(sigma_x)
}


#' Calculate Cancer Cell Fraction (CCF)
#'
#' @description
#' Computes the cancer cell fraction (CCF) based on total copy number,
#' mutation multiplicity, tumor purity, and integer copy number.
#'
#' @param CN Numeric. The total copy number in the segment.
#' @param mu Numeric. Diploid coverage scale factor (e.g. A value of 1 indicates that the segment mean from GATK does not require adjustment. Additionally, the pseudo-diploid coverage is set to 100.).
#' @param rho Numeric. Tumor purity (fraction between 0 and 1).
#' @param C_i Numeric. Integer copy number of the mutation.
#'
#' @details
#' The formula is: \cr
#' \code{ccf = ((C_i * 2 / (mu * 100)) - ((1 - rho) * 2) - (rho * 2)) / (rho * (CN - 2))} \cr
#' If CCF is \code{-Inf}, returns 0.01; if \code{Inf}, returns 1; otherwise returns the computed value.
#'
#' @return Numeric. The estimated cancer cell fraction (CCF), with Inf and -Inf values capped at 1 and 0.01, respectively.
#'
#' @examples
#' Calccf(CN = 4, mu = 100, rho = 0.6, C_i = 2)
#'
#' @export
Calccf <- function(CN, mu, rho, C_i){

  ccf <- ((C_i * 2 / (mu * 100)) - ((1 - rho) * 2) - (rho * 2))/(rho* ((CN ) -2))
  if( !is.na(ccf) && ccf == -Inf){
    ccf <- 0.01
  }else if( !is.na(ccf) && ccf == Inf ){
    ccf <- 1
  }else if( is.na(ccf) ){
    ccf <- 1
  }

  return(ccf)

}


#' Calculate Cancer Cell Fraction (CCF) from Allelic Imbalance (LOH)
#'
#' Estimates the cancer cell fraction (CCF) based on minor and major allele copy numbers, tumor purity, and observed minor-allele frequency (MAF).
#'
#' @param minor Numeric. Minor allele copy number.
#' @param major Numeric. Major allele copy number.
#' @param rho Numeric. Tumor purity (fraction between 0 and 1).
#' @param MAF Numeric. Observed B-allele frequency.
#'
#' @return Numeric. The estimated cancer cell fraction (CCF), with \code{Inf} replaced by 1, \code{-Inf} and \code{NA} replaced by 0.
#'
#' @details
#' The expected MAF is calculated as \code{MAF_ex = minor / (minor + major)}.
#' The formula for CCF is: \cr
#' \code{ccf = (MAF - 0.5) / ((MAF_ex - 0.5) * rho)} \cr
#' If CCF is \code{Inf}, returns 1; if \code{-Inf} or \code{NA}, returns 0.
#'
#' @examples
#' CcfLOH(minor = 1, major = 2, rho = 0.7, MAF = 0.6)
#'
#' @export
CcfLOH <- function(minor, major, rho, MAF  ){
  # ccf from AI information
  # if MAF_ex is 0.5 then ccf is 0
  # if MAF is 0.5 then ccf is 0
  MAF_ex <- minor/(minor + major)

  ccf <- (MAF - 0.5 )/( ( MAF_ex - 0.5) * rho )
  ccf <- ifelse(is.infinite(ccf), ifelse(ccf > 0, 1, 0), ccf)
  ccf <- ifelse( is.na(ccf), 0 , ccf)
  return(ccf)
}



#' Generate All Major/Minor Allele Combinations for a Given Total Copy Number
#'
#' Generates all valid combinations of major and minor allele copy numbers that sum to the specified total copy number, with the constraint that major >= minor.
#'
#' @param CN Integer. The total copy number.
#'
#' @return A data frame with columns:
#'   \item{major}{Major allele copy number.}
#'   \item{minor}{Minor allele copy number.}
#'   \item{CN}{Sum of major and minor (equals input \code{CN}).}
#'
#' @importFrom dplyr mutate
#' @examples
#' GenerateCombinations(3)
#'
#' @export
GenerateCombinations <- function(CN) {
  tmp <- expand.grid(
    major = 0:CN,
    minor = 0:CN
  ) %>%
    subset(major >= minor & major + minor == CN) %>%
    dplyr::mutate( CN = major + minor )
  return(tmp) }

#' Assign Priors to Major/Minor Allele Combinations Based on Biological Difficulty
#'
#' Assigns prior probabilities to major/minor allele combinations using an exponential decay model that incorporates biological plausibility.
#'
#' @param combinations Data frame. Output from [GenerateCombinations()], with columns `major`, `minor`, and `CN`.
#' @param lambda Numeric. Decay rate parameter for the exponential prior.
#'
#' @return A data frame with the input columns plus:
#' \describe{
#'   \item{Bio\_diff}{Assigned biological difficulty score}
#'   \item{prior}{Normalized prior probability for each combination}
#' }
#' @details
#' The function assigns a biological difficulty score to each combination. For unlisted combinations, the score is set according to rules based on \code{CN}, \code{major}, and \code{minor}.
#' The prior is then calculated as \eqn{\exp(-\lambda \times \text{Bio\_diff})}{exp(-lambda * Bio_diff)} and normalized to sum to 1.
#'
#' @importFrom dplyr left_join mutate rowwise ungroup select
#'
#' @examples
#' combos <- GenerateCombinations(3)
#' AssignPriors(combos, lambda = 0.5)
#'
#' @export

AssignPriors <- function(combinations, lambda ) {
  # prior assigned by exponential decay, lambda controls decay rate
  # Complex model will have lower weight

  bio_diff <- data.frame(
    major = c(0, 1, 2, 1, 3, 2, 4, 3, 2, 5, 4, 3),
    minor = c(0, 0, 0, 1, 0, 1, 0, 1, 2, 0, 1, 2),
    Bio_diff = c(3, 2, 3, 1, 4, 2, 5, 3, 4, 6, 4, 4)
  )

  combinations <- combinations %>%
    dplyr::left_join( bio_diff, by = c("major","minor")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate( Bio_diff = ifelse( is.na(Bio_diff), ifelse(minor == 0, CN + 1, ifelse( major == minor, CN, CN - 1) ), Bio_diff ) ) %>%
    dplyr::mutate( Bio_diff_prior = exp(-lambda * Bio_diff) ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate( prior = Bio_diff_prior/sum(Bio_diff_prior,na.rm = T)) %>%
    dplyr::select( -Bio_diff_prior)

  return(combinations)

}



#' Estimate Tumor Purity from Coverage and Segmentation Data
#'
#' Estimates tumor purity by fitting coverage and MAF segmentation data across a range of purity values, selecting the purity that best matches integer copy number models.
#'
#' @param seg A data frame or tibble of segment data, with columns: Chromosome, Segment_Mean, MAF, MAF_gmm_G, MAF_Probes, MAF_gmm_weight, size, etc.
#' @param gender Gender male or femal.
#'
#' @return A list with elements:
#'   \item{dis_df}{Data frame of model fit statistics for each purity.}
#'   \item{Final_model}{List with \code{Final_mu} and \code{Final_rho} (the selected purity).}
#'   \item{min_dis}{Minimum distance to integer copy numbers.}
#'   \item{seg_call}{Segment calls for the selected purity.}
#'   \item{Model_source}{String, always "Coverage".}
#'
#' @details
#' This function fits models for a sequence of purity values, applies segment correction and calling, and selects the purity with the best fit to integer copy number.
#'
#' @importFrom dplyr rowwise mutate group_by summarise arrange filter ungroup
#' @importFrom tidyr unnest_wider
#' @importFrom tibble as_tibble
#' @export
EstimatePurityCov <- function( seg, gender ){
  seg_df <- seg
  purity_array <- seq(0.02, 1, by = 0.02)
  models_cov <- lapply(purity_array,function(p){
    tmp <- seg_df %>%
      dplyr::rowwise() %>%
      dplyr::mutate( correct = list(CorrectPurity(chromosome = Chromosome,
                                           cov_segmentmean = Segment_Mean,
                                           MAF_observe = MAF,
                                           gender = gender,
                                           purity = p ))) %>%
      tidyr::unnest_wider(col = "correct") %>%
      dplyr::rowwise() %>%
      dplyr::mutate( Call = CallwoModel(chromosome = Chromosome,
                                 CNF_correct = CNF_correct,
                                 MAF_correct = MAF_correct,
                                 MAF_gmm_G = MAF_gmm_G,
                                 MAF_Probes = MAF_Probes,
                                 MAF_gmm_weight = MAF_gmm_weight,
                                 gender = gender
                                  )) %>%
      dplyr::mutate( CN = RoundCN(Chrom = Chromosome, Call = Call, CNF = CNF_correct)) %>%
      dplyr::mutate( rho = p )

    tmp_dis <- tmp %>%
      dplyr::group_by( rho ) %>%
      dplyr::filter( abs( Segment_Mean )>=0.1 & size >= 10000000 ) %>%
      dplyr::summarise( dis_integer_CN = sum( abs(CN - CNF_correct)),
                 ploidy = mean(CN) ,.groups = "drop")
    models_cov <- list( results_cov = tmp, results_dis = tmp_dis)

    return(models_cov)
  })
  models_cov <- do.call(rbind,models_cov)
  dis <- do.call(rbind, models_cov[,2])
  results_cov <- tibble::as_tibble(do.call(rbind, models_cov[,1]))
  # select the one with lowest total dis to CN

  dis <- dis %>%
    dplyr::arrange( abs( ploidy - 2 ), dis_integer_CN )
  purity <- as.numeric(dis[1,"rho"])

  # search values in the from the purity range
  purity_range <- round(seq(0.1, 1, length.out = 30),digits = 3)
  purity_range_diff <- abs(purity_range - purity)
  purity <- purity_range[order(purity_range_diff)][1]
  purities <- unique(results_cov$rho)
  purity <- purities[which(purities <= (purity + 0.015) & purities >= (purity - 0.015) )][1]
  CovModel <- list(dis_df = dis,
                   Final_model = list( Final_mu = 1, Final_rho = purity),
                   min_dis = as.numeric(dis[1,"dis_integer_CN"]),
                   seg_call = results_cov %>% dplyr::filter(rho == purity),
                   Model_source = "Coverage")

  return(CovModel)
}




#' Determine Model Source Based on MAF and Coverage Information
#'
#' Determines whether the model source for a segment should be "Coverage", "Diploid", or "Coverage + MAF", based on MAF and segment mean values.
#'
#' @param seg_df A data frame or tibble containing segment information, including columns: \code{size}, \code{FILTER}, \code{MAF_Probes}, \code{MAF}, and \code{Segment_Mean}.
#' @param modelminprobes Numeric. minimum number of probes for a segment to include for model decision
#'
#' @return Character. One of "Coverage", "Diploid", or "Coverage + MAF".
#'
#' @details
#' - If all segments with sufficient probes have minimum MAF > 0.4 and the maximum absolute segment mean is >= 0.2, returns "Coverage".
#' - If all segments with sufficient probes have minimum MAF > 0.4 and the maximum absolute segment mean is < 0.2, returns "Diploid".
#' - Otherwise, returns "Coverage + MAF".
#'
#' @importFrom dplyr filter arrange
#' @examples
#' # seg_df must have columns: size, FILTER, MAF_Probes, MAF, Segment_Mean
#' # ModelSource(seg_df, modelminprobes)
#'
#' @export
ModelSource <- function( seg_df, modelminprobes  ) {
  # Define model source according to MAF
  # If all 0.5 >= MAF > 0.4 and abs(Segment_Mean) >= 0.1, then Coverage otherwise Coverage + MAF
  seg_df <- seg_df %>%
    dplyr::filter( size >= 10000000 )


  min_MAF <- seg_df %>%
    dplyr::filter( FILTER != "FAILED" & MAF_Probes >= modelminprobes ) %>%
    dplyr::arrange(MAF)
  min_MAF <- min(min_MAF$MAF,na.rm = T)
  max_cov <- max(abs(seg_df$Segment_Mean), na.rm = T)


  if( min_MAF > 0.4 ){
    if( max_cov >= 0.2 ){
      model_source <- "Coverage"
    }else{
      model_source <- "Diploid"
    }

  }else{model_source <- "Coverage + MAF"}

  return(model_source)
}



#' Calculate Marginal Likelihood for a Single Segment
#'
#' Computes the marginal likelihood for a segment given copy number, MAF, purity, and model parameters.
#'
#' @param C_i Numeric. Integer copy number for the segment.
#' @param B_i Numeric. Observed B-allele frequency (MAF) for the segment.
#' @param mu Numeric. Diploid coverage scale factor (e.g. A value of 1 indicates that the segment mean from GATK does not require adjustment. Additionally, the pseudo-diploid coverage is set to 100.).
#' @param rho Numeric. Tumor purity (fraction between 0 and 1).
#' @param sigma_C Numeric. Not directly used in this function, but included for compatibility.
#' @param k Numeric. Beta distribution concentration parameter.
#' @param lambda Numeric. Exponential decay parameter for the prior.
#' @param gamma Numeric. Weight for the prior in the likelihood calculation.
#' @param epsilon Numeric. Small value to avoid log(0) and zero parameters in beta.
#'
#' @return A data frame with columns: major, minor, CN, ccf, Bio_diff, prior, expected_MAF, MAF_ll, weighted_prior, exp_MAF_ll, exp_prior, MAF_likelihood, Segcov, MAF, mu, rho.
#'
#' @details
#' The returned data frame contains the following columns:
#' \itemize{
#'   \item major: Major allele count
#'   \item minor: Minor allele count
#'   \item CN: Total copy number
#'   \item ccf: Cancer cell fraction
#'   \item Bio_diff: Biological difference
#'   \item prior: Prior probability
#'   \item expected_MAF: Expected MAF
#'   \item MAF_ll: MAF log-likelihood
#'   \item weighted_prior: Weighted prior
#'   \item exp_MAF_ll: Expected MAF log-likelihood
#'   \item exp_prior: Expected prior
#'   \item MAF_likelihood: MAF likelihood
#'   \item Segcov: Segment coverage
#'   \item MAF: Observed MAF
#'   \item mu: Mutation multiplicity
#'   \item rho: Tumor purity
#' }
#'
#' @export
CalSegmentLikelihood <- function(C_i, B_i,  mu, rho, sigma_C, k, lambda, gamma, epsilon ) {

  #ccf range 0.1-1, so the CNF_tumor range is
  # Calculate C_tumor when ccf = 0.05
  ccf <- 0.2             # Value for ccf
  CN_tumor_max <- (C_i * 2 / (mu * 100) - (1 - rho) * 2 - rho * (1 - ccf) * 2) / (rho * ccf)
  # Calculate C_tumor when ccf = 1
  ccf <- 1                 # Value for ccf
  CN_tumor_min <- (C_i * 2 / (mu * 100) - (1 - rho) * 2 - rho * (1 - ccf) * 2) / (rho * ccf)

  CN_tumor <- c(floor(min(CN_tumor_min,CN_tumor_max)): ceiling(max(CN_tumor_max, CN_tumor_min)))
  CN_tumor <- CN_tumor[which(CN_tumor<= 10 & CN_tumor >= 0 )]
  if( length(CN_tumor) > 0){
    combinations <- GenerateCombinations( CN = CN_tumor[1] )
    combinations<- combinations %>%
      dplyr::rowwise() %>%
      dplyr::mutate( ccf = Calccf( CN = CN, mu = mu, rho = rho, C_i = C_i) )
  }else{combinations <- NA}

  if(length(CN_tumor) > 1){
    for ( tumor_n in CN_tumor[-1]) {
      combinations2 <- GenerateCombinations( CN = tumor_n )
      combinations2 <- combinations2 %>%
        dplyr::rowwise() %>%
        dplyr::mutate( ccf = Calccf( CN = CN, mu = mu, rho = rho, C_i = C_i) )
      combinations <- rbind(combinations,combinations2)
      combinations <- combinations %>%
        dplyr::distinct_all() %>%
        dplyr::filter( major >= 0 & minor >= 0 )
    }
  }
  if( !is.na(combinations)[1] ){


    combinations <- AssignPriors(combinations = combinations, lambda = lambda)
    combinations$MAF_likelihood <- apply(combinations, 1, function(row) {
      m <- as.numeric(row["major"])
      n <- as.numeric(row["minor"])
      prior <- as.numeric(row["prior"])
      ccf <- as.numeric(row["ccf"])

      # Calculate likelihood according to the possible combinations determined from CNF
      # MAF = 0 leads to alpha = 0 , the gamma funcfion will be seperately defined to avoid an infinite log-gamma term [ lgamma(alpha) ]
      # when the total copy number is 2, the ccf value will be NA, so here we assume ccf is 1 for all diploid segments

      if( m + n > 0 ){
        if( ccf < 0.1 & m + n == 2 & n ==0){
          ccf <- 1
        }
        maf <- (1 - rho) * 0.5 + rho * ( ccf*( n / (m + n) ) + (1-ccf) * 0.5 )
        alpha <- k * maf + epsilon
        beta <- k * (1 - maf) + epsilon
        maf_ll <- lgamma(alpha + beta) - lgamma(alpha) - lgamma(beta) +
          (alpha - 1) * log(max(B_i,epsilon,na.rm = T)) + (beta - 1) * log(max(1 - B_i, epsilon, na.rm = T))

        # Introduce a weighted score that combines both the likelihood and the prior
        combine_maf_ll <- maf_ll + gamma * log(prior)
        expected_maf = maf
      }else{
        expected_maf = NA
        maf_ll <- NA
        combine_maf_ll <- NA}


      l <- list( expected_maf = expected_maf,
                 maf_ll = maf_ll,
                 weighted_prior = gamma*log(prior),
                 exp_maf_ll = exp(maf_ll),
                 exp_prior = exp( gamma * log(prior) ),
                 MAF_likelihood = exp( combine_maf_ll ) )

      return(l)
    })

    combinations <- combinations %>%
      tidyr::unnest_wider(col = MAF_likelihood)
    combinations <- combinations %>%
      dplyr::arrange(-MAF_likelihood) %>%
      dplyr::mutate( Segcov = C_i, MAF = B_i, mu= mu, rho = rho)
  }else{
    combinations <- data.frame(
      major = NA,
      minor = NA,
      CN = NA,
      ccf = NA,
      Bio_diff = NA,
      prior = NA,
      expected_maf = NA,
      maf_ll = NA,
      weighted_prior = NA,
      exp_maf_ll = NA,
      exp_prior = NA,
      MAF_likelihood = NA,
      Segcov = C_i,
      MAF = B_i,
      mu = mu,
      rho = rho
    )

  }

  # Marginalize over all (m, n) combinations
  return(combinations)
}



#' Calculate Likelihood for Each \eqn{\mu} and \eqn{\rho} Combination
#'
#' For each segment in the data, computes the MAF likelihood across all major/minor allele combinations for given mutation multiplicity (\eqn{\mu}) and tumor purity (\eqn{\rho}).
#'
#' @param mu Numeric. Diploid coverage scale factor (e.g. A value of 1 indicates that the segment mean from GATK does not require adjustment. Additionally, the pseudo-diploid coverage is set to 100.).
#' @param rho Numeric. Tumor purity (fraction between 0 and 1).
#' @param data Data frame or tibble. Must include columns: \code{Segcov}, \code{MAF}, \code{index}, \code{Tag}, and \code{k}.
#' @param sigma_C Numeric. Parameter for segment likelihood (passed to \code{CalSegmentLikelihood}).
#' @param lambda Numeric. Exponential decay parameter for the prior.
#' @param gamma Numeric. Weight for the prior in the likelihood calculation.
#' @param epsilon Numeric. Small value to avoid log(0) and zero parameters in beta.
#'
#' @return A data frame with the segment likelihoods and all columns from \code{\link{CalSegmentLikelihood}}, plus \code{index} and \code{Tag}.
#'
#' @details
#' For each segment, calls \code{\link{CalSegmentLikelihood}} and filters out combinations with \code{ccf > 1.2}.
#' If there is only one combination and both major and minor are zero, sets \code{MAF_likelihood} to 1.
#'
#' @importFrom dplyr filter mutate
#' @export
Callikelihood <- function(mu, rho, data, sigma_C, lambda, gamma, epsilon) {
  MAF_likelihood <- lapply( 1:nrow(data), function(i) {
    C_i <- data$Segcov[i]
    B_i <- data$MAF[i]
    index_i <- data$index[i]
    Tag_i <- data$Tag[i]
    k_i <- data$k[i]
    segment_likelihood <- CalSegmentLikelihood(C_i = C_i,
                                               B_i = B_i,
                                               mu = mu,
                                               rho = rho,
                                               sigma_C =  sigma_C,
                                               k = k_i,
                                               lambda = lambda,
                                               gamma = gamma,
                                               epsilon= epsilon )
    ## if the ccf is higher than 1.2 will be removed
    segment_likelihood <- segment_likelihood %>%
      dplyr::filter(is.na( ccf ) |ccf <= 1.3 ) %>%
      dplyr::mutate( index = index_i,
              Tag = Tag_i )

    if( nrow(segment_likelihood) == 1 ){
      segment_likelihood <- segment_likelihood %>%
        dplyr::mutate( MAF_likelihood = ifelse( minor == 0 & major ==0 , 1, MAF_likelihood ) )
    }

    return( segment_likelihood)
  })
  MAF_likelihood <- do.call(rbind,MAF_likelihood)

  return(MAF_likelihood)
}



#' Run Likelihood Calculation Across All Purity and Scale Factor Combinations
#'
#' For each combination of purity and scale factor in \code{purity_sf}, computes the segment likelihoods for all segments in \code{data} using the specified model parameters.
#'
#' @param purity_sf Data frame or tibble. Grid of purity and scale factor values to evaluate; must include columns \code{mu} and \code{rho}.
#' @param data Data frame or tibble. Segment data (see \code{\link{Callikelihood}} for required columns).
#' @param sigma_C Numeric. Parameter for segment likelihood (passed to \code{Callikelihood}).
#'
#' @return A data frame with all likelihoods for each segment and each (mu, rho) combination, with NA minor alleles removed and NA ccf values set to 1.
#'
#' @importFrom dplyr filter mutate
#' @examples
#' # RunCallikelihood(purity_sf, data, sigma_C = 0.1, lambda = 1, gamma = 1, epsilon = 0.01)
#'
#' @export
RunCallikelihood <- function( purity_sf, data, sigma_C, lambda = 1, gamma = 1, epsilon = 0.01){

  results_ll <- lapply(1:nrow(purity_sf ), function(n) {
    mu <- purity_sf[n,"mu"]
    rho <- purity_sf[n,"rho"]
    tmp_result <- Callikelihood(mu = mu,
                                rho = rho,
                                data = data,
                                sigma_C = var_sf,
                                lambda = lambda,
                                gamma = gamma,
                                epsilon = epsilon
    )
    return(tmp_result)
  })

  results <- do.call(rbind, results_ll)
  results_refine <- results %>%
    dplyr::filter( !is.na( minor)) %>%
    dplyr::mutate( ccf = ifelse( is.na(ccf), 1, ccf))

  return(results_refine)

}



#' Select the Most Likely Call per Segment
#'
#' For each segment, selects the allele combination (major/minor) with the highest MAF likelihood,
#' with special handling for subclonal events and cases where minor allele is zero.
#'
#' @param results A data frame or tibble containing all possible calls for each segment, with columns including:
#' \code{Segcov}, \code{MAF}, \code{mu}, \code{rho}, \code{index}, \code{major}, \code{minor}, \code{ccf},
#' \code{MAF_likelihood}, \code{Bio_diff}, \code{maf_ll}, \code{cov_diff}, etc.
#'
#' @return A data frame or tibble with the most likely call per segment, including a \code{log_MAF_likelihood} column.
#'
#' @details
#' - Assigns a minimum likelihood for cases where \code{MAF_likelihood} is zero.
#' - If possible subclonal events are present, compares top two likelihoods and may select the second as a subclonal event.
#' - For segments where \code{minor == 0}, the call is refined using coverage difference.
#'
#' @importFrom dplyr filter mutate group_by group_modify arrange slice_max ungroup select rows_update slice_head between row_number n
#' @importFrom tidyr unnest_wider
#' @export
SelectCallpersegment <- function( results ){
  # Assign extremely small likelihood for likelihood 0 term
  likelihood_min <- results %>% dplyr::filter( MAF_likelihood != 0 )
  likelihood_min <- min(likelihood_min$MAF_likelihood, na.rm = T)

  results <- results %>%
    dplyr::mutate(
      expected_cov = mu * 100* (rho * (major + minor) + (1 - rho) * 2 ) / 2,
      cov_diff = abs(Segcov - expected_cov) )
  #1. choose top likelihood allele combinations for each segment
  top_likelihood_rows <- results %>%
    dplyr::filter( ccf <= 1.3) %>%
    dplyr::group_by(Segcov, MAF, mu, rho, index) %>%
    dplyr::group_modify(~ {
      .x %>%
        dplyr::arrange(desc(MAF_likelihood), Bio_diff ) %>%
        dplyr::slice_max(order_by = MAF_likelihood, n = 2, with_ties = TRUE) %>%
        dplyr::mutate(
          MAF_likelihood = dplyr::case_when(
            dplyr::n() == 1 & is.na(MAF_likelihood) ~ 1,
            is.na(MAF_likelihood) ~ likelihood_min,
            TRUE ~ MAF_likelihood
          ),
          choose_second = dplyr::if_else(
            ## if the first call is ref, but the second call is a subclonal has ccf > 0.3 and the maf_ll of the second call
            ## is higher than the first call then choose second one,
            ## when the first call ccf_MAF is higher than 1.3 then choose the second one, this is to avoid choos something is really off from the MAF.
            ( ( minor[1] == 1 & major[1] == 1 &  between(ccf[2], 0.3, 1) & maf_ll[2] >= maf_ll[1] ) ) |
              ( ccf_MAF[1] >= 1.3) ,
            2, 1
          )
        ) %>%
        dplyr::filter(row_number() == choose_second) %>%
        dplyr::select(-choose_second)
    }) %>%
    dplyr::ungroup()


  # refine the ones with minor ==0 with coverage difference
  top_likelihood_rows_minor0 <- top_likelihood_rows %>% dplyr::filter( minor == 0 )

  if(nrow(top_likelihood_rows_minor0) > 0){

    cov_top <- results %>%
      dplyr::filter( minor == 0) %>%
      dplyr::group_by(Segcov, MAF, mu, rho, index) %>%
      dplyr::arrange( cov_diff) %>%
      dplyr::slice_head(n = 1)

    id_cols <- c("index", "Segcov", "MAF", "mu", "rho")
    top_likelihood_rows_minor0 <- top_likelihood_rows_minor0  %>%
      dplyr::rows_update(cov_top, by = id_cols, unmatched = "ignore")

    top_likelihood_rows_minor_non0 <- top_likelihood_rows %>%
      dplyr::filter( minor != 0 )
    top_likelihood_rows <- rbind(top_likelihood_rows_minor0, top_likelihood_rows_minor_non0)

  }
  top_likelihood_rows <- top_likelihood_rows %>%
    dplyr::mutate( log_MAF_likelihood = log(MAF_likelihood) )

  return( top_likelihood_rows )

}



#' Find the "Elbow" Point in a Ranked Vector (Tier 1 Model Selection)
#'
#' Identifies the first major drop ("elbow") in a sorted vector of values, using a smoothed difference approach.
#'
#' @param values Numeric vector. Typically sorted (e.g., decreasing likelihoods or scores).
#' @param fold Numeric. Multiplicative factor for the median drop to define a significant elbow (e.g., 1.5).
#'
#' @return Integer. The index of the first major drop ("elbow") in the \code{values} vector.
#'
#' @details
#' Uses a rolling mean of negative differences to smooth noise and find the first drop exceeding \code{fold} times the median drop.
#' If no such drop is found, uses the geometric "elbow" method as a fallback.
#'
#' @importFrom zoo rollmean
#' @importFrom purrr map_dbl
#' @examples
#' set.seed(1)
#' vals <- sort(runif(20, 0, 10), decreasing = TRUE)
#' FindTier1Models(vals, fold = 1.5)
#'
#' @export
FindTier1Models<- function(values, fold) {
  n <- length(values)
  diffs <- -diff(values)  # Compute negative differences (drops)
  smoothed_diffs <- zoo::rollmean(diffs, k = 3, fill = NA)  # Smooth noise

  # Find the first drop exceeding 1.5x the median drop
  median_drop <- median(smoothed_diffs, na.rm = TRUE)
  nonzero_smallest_drop <- quantile(smoothed_diffs[which(smoothed_diffs != 0 & smoothed_diffs != Inf & smoothed_diffs != -Inf)],na.rm = T)[2]
  if( !is.na(median_drop) & median_drop > 0 ){
    threshold <- fold * median_drop
  }else{
    threshold <- nonzero_smallest_drop
  }

  first_elbow <- which(smoothed_diffs > threshold)[1]

  # Fallback: Use original elbow method if threshold fails
  if (is.na(first_elbow)) {
    distances <- purrr::map_dbl(1:n,
                                ~ abs((values[n] - values[1]) * .x - (n - 1) * values[.x] + n*values[1] - values[n]) /
                           sqrt((values[n] - values[1])^2 + (n - 1)^2))
    first_elbow <- which.max(distances)
  }

  first_elbow

}



#' Select Model by Distance to Integer Copy Number
#'
#' For each candidate model in \code{tier1}, evaluates the mean distance to integer copy number for diploid and non-diploid segments, and returns a summary table sorted by total distance.
#'
#' @param tier1 Data frame or tibble. Each row is a candidate model with columns \code{mu}, \code{rho}, \code{total_log_likelihood}, and \code{segments_n}.
#' @param df Data frame or tibble. All segment call likelihoods (see \code{\link{ExtractCall}}).
#' @param seg Data frame or tibble. Segment information.
#' @param gender gender options to be passed to \code{\link{RefineCalls}}.
#'
#' @return A data frame with one row per model, including columns:
#'   \item{mu}{Mutation multiplicity.}
#'   \item{rho}{Tumor purity.}
#'   \item{total_log_likelihood}{Total log-likelihood for the model.}
#'   \item{segments_n}{Number of segments.}
#'   \item{diploid_n}{Number of diploid segments.}
#'   \item{diploid_distance_to_integer}{Mean distance to integer CN for diploid segments.}
#'   \item{nondiploid_n}{Number of non-diploid segments.}
#'   \item{nondiploid_distance_to_integer}{Mean distance to integer CN for non-diploid segments.}
#'   \item{total_distance_to_integer}{Sum of diploid and non-diploid mean distances.}
#'
#' @importFrom dplyr filter mutate arrange
#' @export
SelectModelByDis <- function(tier1, df, seg, gender){


  tier1_dis <- lapply(1:nrow(tier1), function(index){

    tmp_model <- tier1[index,]
    tmp_mu <- tmp_model$mu
    tmp_rho <- tmp_model$rho
    tmp_call <- ExtractCall(df = top_likelihood_rows, max_L_mu = tmp_mu, max_L_rho = tmp_rho, seg = seg)
    tmp_call <- RefineCalls(df = tmp_call, max_L_mu = tmp_mu, max_L_rho = tmp_rho, gender = gender )
    tmp_call <- tmp_call %>%
      dplyr::mutate(dis_integer_CN = abs(CNF_correct - CN) ) %>%
      dplyr::mutate( dis_integer_CN = ifelse(MAF_Probes <= 50, 0, dis_integer_CN ))

    tmp_diploid_call <- tmp_call %>%
      dplyr::filter( CN == 2 )
    tmp_nondiploid_call <- tmp_call %>%
      dplyr::filter( CN != 2)

    diploid_n = nrow(tmp_diploid_call)
    diploid_distance_to_integer = ifelse( diploid_n > 0, mean(tmp_diploid_call$dis_integer_CN, na.rm = T), 0 )
    nondiploid_n = nrow(tmp_nondiploid_call)
    nondiploid_distance_to_integer =  ifelse( nondiploid_n > 0 , mean(tmp_nondiploid_call$dis_integer_CN, na.rm = T), 0 )
    total_distance_to_integer = diploid_distance_to_integer + nondiploid_distance_to_integer

    tmp_model <-   data.frame ( mu = tmp_mu,
                                rho = tmp_rho,
                                total_log_likelihood = tmp_model$total_log_likelihood,
                                segments_n = tmp_model$segments_n,
                                diploid_n = diploid_n,
                                diploid_distance_to_integer = diploid_distance_to_integer,
                                nondiploid_n = nondiploid_n,
                                nondiploid_distance_to_integer = nondiploid_distance_to_integer,
                                total_distance_to_integer = total_distance_to_integer
    )

    return(tmp_model)
  })

  tier1_dis <- do.call(rbind,tier1_dis)
  tier1_dis <- tier1_dis %>%
    dplyr::arrange( total_distance_to_integer, diploid_distance_to_integer, nondiploid_distance_to_integer)
  return(tier1_dis)

}



#' Estimate Minimum Tumor Purity from MAF Data
#'
#' Estimates the minimum tumor purity based on the lowest MAF among segments tagged as "Include".
#'
#' @param df A data frame or tibble containing at least the columns \code{MAF} and \code{Tag}.
#'
#' @return Numeric. The estimated minimum tumor purity, or 0 if no suitable segments are found.
#'
#' @details
#' Only rows with \code{Tag == "Include"} are considered. If the minimum MAF is less than 0.4, the minimum purity is estimated as \code{max((0.5 - min\_MAF) / 0.5 - 0.05, 0)}. Otherwise, returns 0.
#'
#' @importFrom dplyr filter arrange
#' @examples
#' df <- data.frame(MAF = c(0.3, 0.45, 0.5), Tag = c("Include", "Exclude", "Include"))
#' EstimateMinPurity(df)
#'
#' @export
EstimateMinPurity <- function( df ){

  df <- df %>%
    dplyr::filter( Tag == "Include" ) %>%
    dplyr::arrange( MAF )
  if( nrow(df ) > 0 ){

    min_MAF <- df[1,"MAF"] %>% as.numeric()
    if( min_MAF < 0.4 ){
      min_purity <- (0.5- min_MAF)/0.5
      min_purity <- max(min_purity - 0.05, 0)
    }else{min_purity <- 0}

  }else{min_purity <- 0}
  return(min_purity)
}



#' Refine Tier 1 Models and Their Segment Calls
#'
#' For each candidate (mu, rho) model in \code{tier1}, extracts and refines segment calls, computes fit statistics, and returns a list of calls and model summaries.
#'
#' @param tier1 Data frame or tibble. Each row is a candidate model with columns \code{mu} and \code{rho}.
#' @param results Data frame or tibble. The full set of segment likelihoods/results.
#' @param top_likelihood_rows Data frame or tibble. Top likelihood calls for each segment.
#' @param likelihood_min Numeric. The minimum log-likelihood value to use as a penalty.
#' @param seg Data frame or tibble. Segment information.
#' @param modelminprobes Integer. Minimum number of probes for model inclusion.
#' @param gender Character. Sample gender ("male" or "female").
#' @param callcov Numeric. Subclonal events calling cutoff based on CNF.
#' @param modelminAIsize Numeric. Minimum segment size for model inclusion.
#'
#' @return A list of length equal to \code{nrow(tier1)}. Each element is a list with:
#'   \item{calls}{Refined segment calls for the model.}
#'   \item{model}{Summary statistics for the model (data.frame).}
#'
#' @details
#' For each (mu, rho) combination, extracts and refines calls, computes distances to integer copy number, ploidy, and likelihoods, and summarizes the results.
#'
#' @importFrom dplyr filter mutate
#' @examples
#' # RefineTier1Models(tier1, results, top_likelihood_rows, likelihood_min, seg, modelminprobes = 20, gender = "female", callcov = 0.3, modelminAIsize = 500000)
#'
#' @export
RefineTier1Models <- function( tier1, results, top_likelihood_rows, likelihood_min, seg, modelminprobes, gender, callcov,modelminAIsize  ){

  tier1_refine_tmp <- lapply( 1:nrow(tier1),function(i){
    tmp_mu <- as.numeric(tier1[i,"mu"])
    tmp_rho <- as.numeric(tier1[i,"rho"])

    ## extract and refine indivisual calls
    tmp_call <- ExtractCall(df = top_likelihood_rows,max_L_mu = tmp_mu,max_L_rho = tmp_rho, seg = seg)
    tmp_call <- RefineCalls(df = tmp_call,max_L_mu = tmp_mu, max_L_rho = tmp_rho, gender = gender)
    tmp_call <- RefineCallsSecond(df = tmp_call,results = results,final_mu = tmp_mu, final_rho = tmp_rho, gender = gender, callcov = callcov)
    tmp_model_call <- tmp_call %>%
      dplyr::filter( ! Chromosome %in% c("X","Y")) %>%
      dplyr::filter( size >= modelminAIsize & FILTER != "FAILED"  & MAF_Probes >= modelminprobes ) %>%
      dplyr::filter(!(major == 0 & minor == 0)) %>%
      dplyr::mutate( dis_integer_CN = abs(CNF_correct - CN) ) %>%
      dplyr::mutate( dis_integer_CN = ifelse(MAF_Probes <= modelminprobes, 0, dis_integer_CN )) %>%
      dplyr::mutate( log_MAF_lilelihood = log(MAF_likelihood) )

    ## ploidy

    ploidy <- mean(tmp_model_call$CN,na.rm = T)
    tmp_diploid_call <- tmp_model_call %>%
      dplyr::filter( CN == 2 )
    tmp_nondiploid_call <- tmp_model_call %>%
      dplyr::filter( CN != 2)
    diploid_n = nrow(tmp_diploid_call)
    diploid_distance_to_integer = ifelse( diploid_n > 0, mean(tmp_diploid_call$dis_integer_CN, na.rm = T), 0 )
    nondiploid_n = nrow(tmp_nondiploid_call)
    nondiploid_distance_to_integer =  ifelse( nondiploid_n > 0 , mean(tmp_nondiploid_call$dis_integer_CN, na.rm = T), 0 )
    total_distance_to_integer = diploid_distance_to_integer + nondiploid_distance_to_integer
    segments_n <- nrow(tmp_model_call)
    penalty_n <- tmp_model_call %>% filter( minor == 0 & major == 0 ) %>% nrow()
    total_likelihood <- sum(tmp_model_call$log_MAF_lilelihood,na.rm = T) + likelihood_min * penalty_n
    tmp_model <-   data.frame ( mu =  tmp_mu,
                                rho = tmp_rho,
                                total_log_likelihood = total_likelihood ,
                                segments_n = segments_n,
                                diploid_n = diploid_n,
                                diploid_distance_to_integer = diploid_distance_to_integer,
                                nondiploid_n = nondiploid_n,
                                nondiploid_distance_to_integer = nondiploid_distance_to_integer,
                                total_distance_to_integer = total_distance_to_integer,
                                ploidy = ploidy
    )
    tmp_call$mu <- tmp_mu
    tmp_call$rho <- tmp_rho
    re <- list( calls = tmp_call, model = tmp_model)
    return(re)
  })

}


#' Group Values by Change Points in the Mean
#'
#' Uses the \code{changepoint} package to detect change points in the mean of a numeric vector and assigns group IDs accordingly.
#'
#' @param x Numeric vector.
#'
#' @return An integer vector of the same length as \code{x}, where each element indicates the group assignment based on detected change points.
#'
#' @details
#' Uses \code{cpt.mean} from the \pkg{changepoint} package with the PELT method and MBIC penalty to detect change points in the mean.
#'
#' @importFrom changepoint cpt.mean cpts
#' @examples
#' # Example with a simulated change in mean
#' set.seed(123)
#' x <- c(rnorm(50, 0), rnorm(50, 3))
#' group_id <- Groupvalues(x)
#' table(group_id)
#'
#' @export
Groupvalues <- function(x) {
  # Load required package

  # Detect change points in the mean
  cpt <- cpt.mean(x, method = "PELT", penalty = "MBIC" )
  change_points <- cpts(cpt)

  # Assign group IDs
  group_id <- rep(1:(length(change_points) + 1),
                  times = diff(c(0, change_points, length(x))))

  return(group_id)
}




#' Cluster and Select the Best Models Based on Likelihood and Copy Number Distance
#'
#' Clusters models by total likelihood, diploid and nondiploid distance to integer copy number, and selects the best (mu, rho) model.
#'
#' @param models_dis A data frame or tibble with model metrics, including columns: \code{total_log_likelihood}, \code{diploid_distance_to_integer}, \code{nondiploid_distance_to_integer}, \code{ploidy}, \code{mu}, \code{rho}, \code{diploid_n}, etc.
#'
#' @return A list with:
#'   \item{Final_mu}{Best model's mutation multiplicity.}
#'   \item{Final_rho}{Best model's tumor purity.}
#'   \item{total_distance_to_integer}{Total distance to integer copy number for the best model.}
#'   \item{total_likelihood}{Total log-likelihood for the best model.}
#'   \item{models}{The full, clustered, and sorted model data frame.}
#'
#' @details
#' Models are clustered using change point detection on likelihood and distance metrics. The best model is selected by a series of sort criteria.
#'
#' @importFrom dplyr arrange mutate group_by ungroup filter rowwise select desc
#' @importFrom changepoint cpt.mean cpts
#' @examples
#' # models_dis <- ... # output from SelectModelByDis
#' # ClusterModels(models_dis)
#'
#' @export
ClusterModels <- function( models_dis ){

  ## Group models based on total likelihood, diploid distance to integer and nondiploid distance to integer when there is
  ## andy nondiploid regions
  models_dis <- models_dis %>%
    dplyr::arrange( desc(total_log_likelihood)  )%>%
    dplyr::mutate( total_likelihood_cluster = Groupvalues(total_log_likelihood)) %>%
    dplyr::arrange( diploid_distance_to_integer ) %>%
    dplyr::mutate( diploid_distance_cluster = Groupvalues( scale(diploid_distance_to_integer)[, 1])) %>%
    dplyr::arrange( nondiploid_distance_to_integer) %>%
    dplyr::mutate( nondiploid_distance_cluster = Groupvalues( scale(nondiploid_distance_to_integer)[, 1])) %>%
    dplyr::group_by(total_likelihood_cluster) %>%
    dplyr::mutate(total_likelihood_cluster_mean = mean(total_log_likelihood)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(diploid_distance_cluster) %>%
    dplyr::mutate(diploid_distance_cluster_mean = mean(diploid_distance_to_integer)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(nondiploid_distance_cluster) %>%
    dplyr::mutate(nondiploid_distance_cluster_mean = mean(nondiploid_distance_to_integer)) %>%
    dplyr::ungroup()

  ## if the distance is 0 due to no segments then group the second group into the first one

  diploid_distance_cluster1_check <- models_dis %>%
    dplyr::filter( diploid_distance_cluster == 1  & diploid_n > 0 ) %>%
    nrow()


  nondiploid_distance_cluster1_check <- models_dis %>%
    dplyr::filter( nondiploid_distance_cluster == 1 & diploid_n > 0  ) %>%
    nrow()

  if(diploid_distance_cluster1_check <= 3  ){
    models_dis <- models_dis %>%
      dplyr::mutate( diploid_distance_cluster =
                ifelse( diploid_distance_cluster > 1, diploid_distance_cluster - 1, diploid_distance_cluster ) )

  }
  if(nondiploid_distance_cluster1_check <= 3  ){
    models_dis <- models_dis %>%
      dplyr::mutate( nondiploid_distance_cluster =
                ifelse( nondiploid_distance_cluster > 1, nondiploid_distance_cluster - 1, nondiploid_distance_cluster ))
  }


  models_dis <- models_dis %>%
    dplyr::arrange( total_likelihood_cluster,
             diploid_distance_cluster ,
             nondiploid_distance_cluster,
             abs( ploidy - 2),
             desc(total_log_likelihood), desc( rho ) )

  final_diploid_distance_cluster <- models_dis[1,"diploid_distance_cluster"] %>% as.numeric()
  final_nondiploid_distance_cluster <- models_dis[1,"nondiploid_distance_cluster"] %>% as.numeric()

  if( final_diploid_distance_cluster !=1 || final_nondiploid_distance_cluster !=1){
    models_dis <- models_dis %>%
      dplyr::rowwise() %>%
      dplyr::mutate( total_likelihood_cluster_update = ifelse( diploid_distance_cluster ==1 &
                                                          nondiploid_distance_cluster ==1 &
                                                          total_likelihood_cluster < 5 &
                                                          total_log_likelihood > 0 ,
                                                        1,total_likelihood_cluster ))
    models_dis <- models_dis %>%
      dplyr::arrange( total_likelihood_cluster_update,
               diploid_distance_cluster ,
               nondiploid_distance_cluster,
               abs( ploidy - 2),
               desc(total_log_likelihood), desc( rho ) ) %>%
      dplyr::select( -total_likelihood_cluster_update)
  }

  final_mu <- as.numeric( models_dis[1,"mu"] )
  final_rho <- as.numeric(models_dis[1,"rho"] )
  Final_model <- list(Final_mu = final_mu,
                      Final_rho = final_rho,
                      total_distance_to_integer = as.numeric( models_dis[1,"total_distance_to_integer"] ),
                      total_likelihood = as.numeric(models_dis[1,"total_log_likelihood"]),
                      models = models_dis)

  return( Final_model )
}


#' Select and Refine the Final CNV Model and Segment Calls
#'
#' Selects the best model (mu, rho) combination based on total log-likelihood and distance to integer copy number, refines segment calls, and writes results to output files.
#'
#' @param results Data frame or tiblle. Likelihood results calculated by RunCallikelihood function.
#' @param top_likelihood_rows Data frame or tibble. Top likelihood calls for each segment, including columns \code{mu}, \code{rho}, \code{log_MAF_likelihood}, \code{major}, \code{minor}, \code{Tag}, etc.
#' @param groupinfo Data frame or tibble. Information for estimating minimum purity (see \code{\link{EstimateMinPurity}}).
#' @param seg Data frame or tibble. Segment information.
#' @param prefix Character. Output file prefix.
#' @param gender Character. Sample gender ("male" or "female").
#' @param out_dir Character. Output directory for saving results.
#' @param model_source Character. Description of model source (e.g., "Coverage", "Coverage + MAF").
#' @param callcov Numeric. Subclonal events calling cutoff based on coverage.
#' @param modelminAIsize Numeric. Minimum segment size for model inclusion.
#' @param modelminprobes Integer. Minimum number of probes for model inclusion.
#' @param minsf Numeric. Minimum scale factor will be estimated. default is 0.4. range 0~0.4.
#'
#'
#' @return A list with:
#'   \item{Final_model}{The final selected model (see \code{\link{ClusterModels}}).}
#'   \item{models}{The full model table (before and after refinement).}
#'   \item{refined_calls}{The final calls of tier1 models.}
#'
#' @details
#' - Refines calls and models through several selection and clustering steps.
#' - Writes the model table and final calls to output files in the specified output directory.
#'
#' @importFrom dplyr filter group_by summarise left_join mutate arrange desc n
#' @importFrom tidyr replace_na
#' @importFrom data.table setDT
#' @importFrom utils write.table
#' @examples
#' # SelectFinalModel(top_likelihood_rows, groupinfo, prefix = "Sample1", gender = "female", out_dir = "results/",model_source = "Coverage + MAF", callcov = 0.3,modelminAIsize=5000000, modelminprobes= 20)
#'
#' @export
SelectFinalModel <- function( results, top_likelihood_rows, groupinfo, prefix, gender, out_dir,model_source,seg, callcov,modelminAIsize, modelminprobes,minsf ){
  # Remove all possible homozygous deletion segments

  min_purity <- EstimateMinPurity(groupinfo)
  tmp <- top_likelihood_rows

  likelihood_min <- min(tmp[which(tmp$MAF_likelihood !=0),"log_MAF_likelihood"], na.rm = T)
  # Penalize segments with >= modelminAIsize and non failed tag but could not estimate likelihood
  penalty_counts <- tmp %>%
    dplyr::filter(major == 0 & minor == 0 & Tag == "Include") %>%
    dplyr::group_by(mu, rho) %>%
    dplyr::summarise(Likelihood_penalty_rows = dplyr::n(), .groups = "drop")


  models <- tmp %>%
    dplyr::filter(!(major == 0 & minor == 0)) %>%
    dplyr::filter(Tag == "Include") %>%
    dplyr::group_by(mu, rho) %>%
    dplyr::summarise(
      total_log_likelihood = sum(log_MAF_likelihood, na.rm = TRUE),
      segments_n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::left_join(penalty_counts, by = c("mu", "rho")) %>%
    dplyr::mutate(
      Likelihood_penalty_rows = replace_na(Likelihood_penalty_rows, 0),
      total_log_likelihood = total_log_likelihood + likelihood_min * Likelihood_penalty_rows
    ) %>%
    dplyr::arrange(desc(total_log_likelihood))


  # Second round of model selection
  # 1. Find Tier1 models
  # 2. Estimate the distance between CNF_correct to integer value then calculate total distance of all segments.
  # 3. Prioritize models with lowest distance


  tier1_index <- 1
  fold <- 3
  while( tier1_index <= 50 && fold <= 15 ){
    tier1_index <- FindTier1Models(values = models$total_log_likelihood, fold = fold)
    fold <- fold + 1
  }
  if(tier1_index <= 50){
    tier1_index <- 50
  }

  tier1 <- models[c(1:tier1_index),]
  ## Refine Tier1 model calls by using RefineCalls and RefineCallsecond function

  tier1_calls <- RefineTier1Models(tier1 = tier1,
                                   results = results,
                                   top_likelihood_rows =  top_likelihood_rows,
                                   likelihood_min = likelihood_min,
                                   seg = seg,
                                   modelminprobes = modelminprobes,
                                   gender = gender,
                                   callcov = callcov,
                                   modelminAIsize = modelminAIsize
                                   )

  refined_calls <- do.call(rbind, lapply(tier1_calls, function(x) x$calls))
  refined_models <- do.call(rbind, lapply(tier1_calls, function(x) x$model))

  PlotModelDot( models = models , out_dir = out_dir ,prefix = prefix, tier1_index = tier1_index)

  ## remove models with lowever purity determined by EstimateMinPurity function
  refined_models$Tier1 <- "Tier1_Models"
  tier1_dis_fil <- refined_models %>%
    dplyr::filter( rho >= min_purity) %>%
    dplyr::filter( mu >= minsf )

  Final_model <- ClusterModels(models_dis = tier1_dis_fil )


  ## udpate the models table

  models <- models %>% dplyr::left_join(Final_model$models, by = c("mu","rho","segments_n"))
  colnames(models) <- gsub("\\.x","_before_refine",colnames(models))
  colnames(models) <- gsub("\\.y","_after_refine",colnames(models))
  models[which(models$mu == as.numeric(Final_model$Final_mu) &
                 models$rho == as.numeric(Final_model$Final_rho)),"Tier1"] <- "Final_model_MAF"

  setDT(refined_calls)
  Final_call <- refined_calls[mu == Final_model$Final_mu & rho == Final_model$Final_rho]
  Final_call$Model_source <- model_source
  OutModelFile <- paste0(out_dir, "/", prefix, "_Models_likelihood.tsv")
  write.table(models, OutModelFile, row.names = F,quote = F, sep = "\t" )
  outms <- paste0( "The model table is saved at: ", OutModelFile)
  print(outms)
  OutFile_calls <- paste0(out_dir, "/", prefix, "_final_calls.tsv")
  write.table(Final_call, OutFile_calls, row.names = F,quote = F, sep = "\t" )
  outms <- paste0( "The final call is saved at: ", OutFile_calls)
  print(outms)
  OutFile_calls <- paste0(out_dir, "/", prefix, "_refined_tier1model_calls.tsv")
  write.table(refined_calls, OutFile_calls, row.names = F,quote = F, sep = "\t" )
  outms <- paste0( "The Tier1 model's call is saved at: ", OutFile_calls)
  print(outms)
  return_models <- list( Final_model = Final_model, models = models, refined_calls = refined_calls )
}

#' Extract Final Calls for a Given (mu, rho) Model
#'
#' For a given (mu, rho) combination, joins the segment call results to the segment metadata and returns the final calls.
#'
#' @param df Data frame or tibble. Contains segment call results, including columns \code{mu}, \code{rho}, and \code{index}.
#' @param max_L_mu Numeric. The selected value of mutation multiplicity (\code{mu}).
#' @param max_L_rho Numeric. The selected value of tumor purity (\code{rho}).
#' @param seg Data frame or tibble. Segment metadata, must include column \code{index}.
#'
#' @return A data frame or tibble with the final calls for the given (mu, rho) model, joined to segment metadata.
#'
#' @importFrom dplyr filter right_join select
#' @examples
#' # ExtractCall(df, max_L_mu = 1, max_L_rho = 0.7, seg)
#'
#' @export
ExtractCall <- function( df, max_L_mu, max_L_rho, seg  ){
  seg$index <- as.character(seg$index)
  df$index <- as.character(df$index)
  final_call <- df %>%
    dplyr::filter( mu == max_L_mu & rho == max_L_rho  ) %>%
    dplyr::right_join( seg, by = "index" )

  if( "Segcov.x" %in% colnames(final_call) ){
    final_call <- final_call %>% dplyr::select(-Segcov.x)
  }
  if( "MAF.x" %in% colnames(final_call) ){
    final_call <- final_call %>% dplyr::select(-MAF.x)
  }

  colnames(final_call) <- gsub("\\..*","",colnames(final_call))

  return(final_call)
}



#' Refine Segment Calls According to the Final Model
#'
#' Refines values and calls for each segment based on the selected (mu, rho) model, adjusting segment mean, purity, and copy number calls.
#'
#' @param df Data frame or tibble. Segment call data, must include columns used in refinement.
#' @param max_L_mu Numeric. The selected value of mutation multiplicity (\code{mu}).
#' @param max_L_rho Numeric. The selected value of tumor purity (\code{rho}).
#' @param gender Character. Gender female or male.
#'
#' @return A data frame or tibble with refined segment calls, ordered by chromosome and start position.
#'
#' @importFrom dplyr rowwise mutate select arrange
#' @importFrom tidyr unnest_wider
#' @examples
#' # RefineCalls(df, max_L_mu = 1, max_L_rho = 0.7, gender = "female")
#'
#' @export
RefineCalls<- function( df , max_L_mu, max_L_rho, gender){
  ## Refine the values and calls according to the final model
  chrom_levels <- c(c(1:22,"X","Y"))
  col_name <- c("Chromosome", "Start", "End","size", "Num_Probes", "Call", "ccf", "ccf_MAF",
                "Segment_Mean", "CNF_correct", "major", "minor", "CN",
                "MAF", "MAF_correct", "expected_maf", "expected_cov", "MAF_Probes",
                "MAF_gmm_G", "MAF_gmm_weight","balance_tag", "BreakpointSource", "FILTER",
                "maf_ll", "MAF_likelihood", "mu", "rho",  "index",
                "gatk_SM_raw", "gatk_count", "gatk_baselinecov",
                "gatk_gender", "pipeline_gender"
  )
  final_call <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate( Segment_Mean = CalSM( max_L_mu = max_L_mu,
                                  Segcov = Segcov,
                                  gender = gender ,
                                  Chromosome = Chromosome) ) %>%
    dplyr::mutate( Correct_purity = list(CorrectPurity(cov_segmentmean = Segment_Mean,
                                                MAF_observe = MAF,
                                                gender = gender,
                                                purity = max_L_rho,
                                                chromosome = Chromosome))) %>%
    tidyr::unnest_wider(col = Correct_purity) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(Call = CallWTModel(major = major,
                              minor = minor,
                              Chromosome = Chromosome,
                              gender = gender,
                              CNF_correct = CNF_correct) ) %>%
    dplyr::mutate( CN = ifelse( !(is.na(minor) | is.na(major)),
                         major + minor ,
                         RoundCN(Chrom = Chromosome, Call = Call, CNF = CNF_correct, gender = gender ))) %>%
    dplyr::select( all_of( col_name )  ) %>%
    dplyr::mutate( Chromosome = factor(Chromosome,levels = chrom_levels)) %>%
    dplyr::arrange( Chromosome, Start ) %>%
    dplyr::mutate( ccf = ifelse( CN == 2 & ccf < 0.1 , 1, ccf )) %>%
    dplyr::mutate( FILTER = ifelse( Chromosome == "Y", "FAILED", FILTER))

  return(final_call)
}



#' Second Refinement of Segment Calls for Coverage/AI Profile Mismatches
#'
#' Further refines segment calls where the coverage profile and allelic imbalance (AI) profile do not match, by considering alternative copy number states and updating calls and expected values.
#'
#' @param df Data frame or tibble. Refined segment calls (see \code{\link{RefineCalls}}), must include columns as in \code{col_name}.
#' @param results Data frame or tibble. All likelihood results for possible copy number states, must include columns \code{rho}, \code{mu}, \code{index}, \code{CN}, etc.
#' @param final_mu Numeric. The selected value of scale factor (\code{mu}).
#' @param final_rho Numeric. The selected value of tumor purity (\code{rho}).
#' @param gender Character. Sample gender ("male" or "female").
#' @param callcov Numeric. Subclonal events calling cutoff based on coverage.
#'
#' @return A data frame or tibble with further refined segment calls, including updated columns such as \code{ccf_COV}, \code{ccf_final}, and \code{CN_mix}.
#'
#' @details
#' Segments with discordant coverage and AI profiles are re-evaluated by searching for alternative copy number states in the likelihood results, and calls are updated accordingly.
#'
#' @importFrom dplyr filter select arrange rowwise mutate relocate
#' @examples
#' # RefineCallsSecond(df, results, final_mu = 1, final_rho = 0.7, gender = "female", callcov = 0.3)
#'
#' @export
RefineCallsSecond <- function( df, results, final_mu, final_rho, gender, callcov ){
  # refine individual call that not matching cov profile and AI profile
  col_name <- c("Chromosome", "Start", "End","size", "Num_Probes", "Call", "ccf", "ccf_MAF",
                "Segment_Mean", "CNF_correct", "major", "minor", "CN",
                "MAF", "MAF_correct", "expected_maf", "expected_cov", "MAF_Probes",
                "MAF_gmm_G", "MAF_gmm_weight","balance_tag",
                "BreakpointSource", "FILTER",
                "maf_ll", "MAF_likelihood", "mu", "rho",
                "index", "gatk_SM_raw", "gatk_count", "gatk_baselinecov",
                "gatk_gender", "pipeline_gender"  )
  df <- df %>%
    dplyr::mutate( cov_diff = abs( CNF_correct - CN) )
  bad <- df %>%
    dplyr::filter( ( cov_diff >= 0.6 & ccf <= 0.8 & ! is.na(MAF) ) |
              ( cov_diff >= 0.6 & ccf == 1 & ! is.na(MAF) ) |
              ( cov_diff >= callcov & CN == 2 & !is.na(MAF)))
  if( nrow(bad) > 0 ){
    good <- df %>%
      dplyr::filter( ! index %in% bad$index)
    for (i in c(1:nrow(bad))) {
      cnf_correct <- bad[i,"CNF_correct"] %>% as.numeric()
      cn_total <- bad[i,"CN"] %>% as.numeric()
      cn_expect <- setdiff( c(floor(cnf_correct), ceiling(cnf_correct)), cn_total ) %>%
        unlist() %>%
        as.numeric()
      cn_expect <- ifelse(length(cn_expect) >1 & cnf_correct >= (1+callcov) | cnf_correct <= (2-callcov) , cn_expect[which(cn_expect != 2)], cn_expect )
      chrom <- bad[i,"Chromosome"] %>% as.character()
      tmp_index <- as.character( bad[i,"index"] )
      tmp <- results %>%
        dplyr::filter( rho == final_rho & mu == final_mu )
      tmp$index <- as.character(tmp$index)
      tmp <- tmp %>%
        dplyr::filter(index == tmp_index & CN == cn_expect)
      if(nrow(tmp)> 0 ){
        tmp <-  tmp %>%
          dplyr::mutate( Chromosome = chrom) %>%
          dplyr::mutate( CNF_correct = cnf_correct ) %>%
          dplyr::rowwise() %>%
          dplyr::mutate( Call = CallWTModel(major = major,
                                     minor = minor,
                                     Chromosome = Chromosome,
                                     gender = gender,
                                     CNF_correct = CNF_correct) )


        if( nrow(tmp) > 0){
          col_inter <- intersect(colnames(tmp), colnames(bad))
          col_inter <- unique(c(col_inter,"MAF","CN"))
          bad[i,col_inter] <- tmp[1,col_inter ]
          bad[i,"expected_cov"] <- final_mu * 100* (final_rho * tmp[1,"CN"] + (1 - final_rho) * 2 ) / 2
        }
      }


    }
    df <- rbind(good, bad)
  }
  df <- df %>%
    dplyr::select(all_of(col_name)) %>%
    dplyr::arrange( Chromosome, Start ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate( ccf = ifelse( ccf > 1, 1 , ccf)) %>%
    dplyr::mutate( ccf_MAF = ifelse( ccf_MAF > 1, 1, ccf_MAF )) %>%
    dplyr::mutate( CN_mix = ifelse( CN %% 2 != 0 & abs(ccf - ccf_MAF) >= 0.3 & ccf_MAF != 0  , "CN_Mix", "No")) %>%
    dplyr::mutate( ccf_final = ifelse( minor == 0 , ccf_MAF, ccf ) )

  colnames(df)[7] <- "ccf_COV"

  df <- df %>%
    dplyr::relocate( ccf_final, .after = "ccf_MAF"  )
  return(df)
}

