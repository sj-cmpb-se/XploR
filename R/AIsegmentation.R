#' Check Whether Two AI Segments Should Be Merged
#'
#' This function determines whether two adjacent AI (allelic imbalance) segments should be merged,
#' based on the difference in their GMM means and quality thresholds.
#'
#' @param cur_row A data frame row (list or tibble row) representing the current segment.
#' @param next_row A data frame row (list or tibble row) representing the next segment.
#' @param opt A list of options, must include \code{mergeai} (numeric threshold for mean difference) and \code{snpmin} (minimum SNP count).
#'
#' @return Logical value: \code{TRUE} if the two segments should be merged, \code{FALSE} otherwise.
#'
#' @export

MergeAICheck <- function(cur_row,next_row,opt){
  # Merge if
  # ai segment diff <= mergeai or one of the ai segment quality is fail
  if(abs( as.numeric(cur_row$gmm_mean) - as.numeric(next_row$gmm_mean)) <= opt$mergeai ){
    re <- TRUE
  }else if( !is.na(cur_row$gmm_weight) & !is.na(next_row$gmm_weight)){
    if(min(cur_row$gmm_weight , next_row$gmm_weight ) <= 0.35 ||
       min(cur_row$nonzero_count, next_row$nonzero_count) < 2*opt$snpmin ){ re <- TRUE }else{re <- FALSE}
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



#' Estimate BAF Cluster Mean and Weight Using Gaussian Mixture Modeling
#'
#' Fits a Gaussian mixture model (GMM) to B-allele frequency (BAF) values and returns the most prominent cluster's mean and weight.
#' Clusters with similar means (within a threshold) are merged. If there are insufficient values or nearly no variance, a simple mean is returned.
#'
#' @param baf_values Numeric vector of BAF values.
#'
#' @return Named numeric vector with elements:
#'   \item{gmm_mean}{The mean of the most prominent cluster (capped at 0.5).}
#'   \item{gmm_weight}{The mixture weight of the most prominent cluster.}
#'   \item{gmm_G}{Number of clusters after merging adjacent means.}
#'
#' @details
#' Uses \code{\link[mclust]{Mclust}} for Gaussian mixture modeling and \code{\link{ClusterAdjacent}} to merge similar clusters.
#'
#' @importFrom mclust Mclust
#' @export
EstimateBAFbyGMM <- function(baf_values){

  n_bafs <- length(baf_values)
  sd <- sd( baf_values, na.rm = T )
  if( n_bafs >= 3 && sd >= 1e-6 ){
    gmm_model <- mclust::Mclust(baf_values, G = 1:10)
    means <- gmm_model$parameters$mean
    weights <- gmm_model$parameters$pro
    if(gmm_model$G >= 2){
      re <- ClusterAdjacent(values = means, weights = weights, threshold = 0.01)
      means <- re$means
      weights <- re$weights
      #variances <- gmm_model$parameters$variance$sigmasq
      num_components <- length(re$means)
    }else{
      num_components <- gmm_model$G
    }

    ## cluster the means if the difference is small


    sorted_indices <- order(-weights)
    gmm_mean <- means[sorted_indices][1] %>% as.numeric()
    gmm_mean <- ifelse( gmm_mean > 0.48, 0.5, gmm_mean )
    gmm_mean <- min( gmm_mean, 0.5 )
    # gmm_variance <- variances[sorted_indices][1]
    gmm_weight <- weights[sorted_indices][1]

    Estimated_baf <- c(
      "gmm_mean" = gmm_mean,
      #"gmm_variance" = gmm_variance,
      "gmm_weight" = gmm_weight,
      "gmm_G" = num_components
    )}else{
      gmm_mean <- mean(baf_values,na.rm = T)
      gmm_mean <- ifelse( gmm_mean > 0.48, 0.5, gmm_mean )
      gmm_mean <- min( gmm_mean, 0.5 )
      Estimated_baf <- c(
        "gmm_mean" =  gmm_mean,
        #"gmm_variance" = sd,
        "gmm_weight" = 1,
        "gmm_G" = 0)

    }

  return(Estimated_baf)

}


#' Merge Adjacent AI Segments Based on Similarity Criteria
#'
#' Iteratively merges adjacent rows in a data frame of AI segments if they meet criteria defined by \code{MergeAICheck}.
#' For merged segments, BAF values are re-estimated and segment information is updated.
#'
#' @param df A data frame or tibble of AI segments, with columns such as Chromosome, bin, Start, End, nonzero_count, each_baf, etc.
#' @param opt A list of options, must include \code{mergeai} (threshold for merging) and \code{snpmin} (minimum SNP count).
#' @param tmp_baf A data frame or tibble containing BAF values and genomic coordinates (must include Chromosome, Start, End, baf).
#'
#' @return A data frame or tibble with merged AI segments, updated BAF estimates, and segment information.
#'
#' @details
#' This function uses \code{\link{MergeAICheck}} to determine if two adjacent segments should be merged.
#' For merged segments, BAF values are re-estimated using \code{\link{EstimateBAFbyGMM}}.
#'
#' @importFrom dplyr filter arrange
#' @importFrom tibble tibble
#' @importFrom tidyr unnest_wider
#' @export
MergeAIRow <- function(df, opt, tmp_baf  ) {

  if(nrow(df) > 1 ){
    i <- 1
    while ( i < ( nrow(df) ) ) {

      cur_row <- df[i,]
      next_row <- df[ i+1,]

      if (  MergeAICheck(cur_row = cur_row, next_row = next_row, opt) ) {
        #cur_row_bafs = str_split(cur_row$each_baf, pattern = ";", simplify = TRUE) %>% as.numeric() %>% na.omit()
        #next_row_bafs = str_split(next_row$each_baf, pattern = ";", simplify = TRUE ) %>% as.numeric() %>% na.omit()
        cur_next_bafs <- tmp_baf %>%
          dplyr::filter( Chromosome == cur_row$Chromosome) %>%
          dplyr::filter( Start >= cur_row$Start & End <= next_row$End)

        new_df <- tibble::tibble(
          Chromosome = cur_row$Chromosome,
          bin = paste(c(cur_row$bin, next_row$bin),collapse = ";"),
          Start = cur_row$Start,
          End = next_row$End,
          #Estimated_baf= list(EstimateBAFbyGMM(baf_values = c(cur_row_bafs, next_row_bafs))) ,
          Estimated_baf = list(EstimateBAFbyGMM(baf_values = cur_next_bafs$baf )),
          nonzero_count = cur_row$nonzero_count + next_row$nonzero_count,
          each_baf = paste( cur_next_bafs$baf , collapse = ";")
        ) %>%
          tidyr::unnest_wider(col = Estimated_baf)

        df <- rbind(new_df, df[-c(i,i+1),])
        df <- df %>% arrange(by=Chromosome,Start,na.last = T)
        i <- 1
      }else{i <- i+1}
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
#' @param opt A list of options to be passed to the merge functions.
#' @param tmp_baf A data frame for use with \code{MergeAIRow}. Not used if \code{AIorSeg == "Seg"}.
#'
#' @return A data frame or tibble with merged segments or AI rows for each chromosome.
#'
#' @details
#' Uses \code{\link{MergeAIRow}} if \code{AIorSeg == "AI"}, otherwise \code{\link{MergeSegRow}}.
#'
#' @importFrom dplyr arrange
#' @export
CallMerge <- function(data,AIorSeg,opt, tmp_baf){
  data<- data %>%
    dplyr::arrange(by=Chromosome, Start, na.last = T)
  bychr <- split(data,f=data$Chromosome)
  if(AIorSeg == "AI"){
    bychr <- lapply(bychr,function(x) MergeAIRow(df = x,opt, tmp_baf = tmp_baf))}
  if(AIorSeg == "Seg"){
    bychr <- lapply(bychr,function(x) MergeSegRow(df = x,opt))
  }
  processed_df <- do.call(rbind,bychr)
  return(processed_df)
}


#' Bin BAF Data and Estimate BAF Metrics per Bin
#'
#' Bins input BAF data into fixed-size windows and calculates summary metrics (e.g., estimated BAF, nonzero counts) for each bin.
#'
#' @param data A data frame or tibble containing at least \code{Chromosome}, \code{End}, and \code{baf} columns.
#' @param bin_size Numeric. The bin size to use (e.g., 1e6 for 1Mb bins).
#'
#' @return A tibble with one row per bin and columns:
#'   \item{Chromosome}{Chromosome identifier.}
#'   \item{bin}{Bin start coordinate.}
#'   \item{Start}{Minimum \code{End} value in the bin.}
#'   \item{End}{Maximum \code{End} value in the bin.}
#'   \item{Estimated_baf}{List column with GMM-based BAF estimation results.}
#'   \item{nonzero_count}{Number of nonzero BAF values in the bin.}
#'   \item{each_baf}{Semicolon-separated string of nonzero BAF values in the bin.}
#'   \item{gmm_mean, gmm_weight, gmm_G}{Unnested BAF GMM metrics.}
#'
#' @importFrom dplyr mutate group_by summarise ungroup filter
#' @importFrom tidyr unnest_wider
#' @export
binbaf <- function(data, bin_size){

  # Bin the data and calculate the desired metrics
  binned_data <- data %>%
    dplyr::mutate(bin = floor(End / bin_size) * bin_size) %>%
    dplyr::group_by(Chromosome, bin) %>%
    dplyr::summarise(
      Start = min(End),
      End = max(End),
      Estimated_baf = list(EstimateBAFbyGMM(baf_values = baf[baf != 0] )),
      nonzero_count = sum(baf != 0),
      each_baf = paste(baf[baf != 0 ], collapse = ";")
    ) %>% tidyr::unnest_wider(col = "Estimated_baf")%>%
    dplyr::ungroup() %>%
    dplyr::filter( !is.na( gmm_mean))

  return(binned_data)

}


#' Estimate BAF Shift Scaling Factor
#'
#' Estimates a scaling factor to correct B-allele frequency (BAF) shift using Gaussian mixture modeling.
#' The function selects BAF values greater than 0.4, fits a GMM (1-3 components), and computes the scaling factor as 0.5 divided by the maximum GMM mean.
#' If there are fewer than 500 SNPs above the threshold, returns \code{NA}.
#'
#' @param baf A data frame or tibble with at least a \code{baf} column (numeric BAF values).
#'
#' @return Numeric. The estimated BAF scaling factor, or \code{NA} if there are not enough SNPs.
#'
#' @details
#' Uses \code{\link[mclust]{Mclust}} for Gaussian mixture modeling.
#'
#' @importFrom dplyr filter
#' @importFrom mclust Mclust
#' @export
EstimateBAFshift <- function( baf){

  tmp_baf <- baf %>%
    dplyr::filter( baf > 0.4)
  if( nrow( tmp_baf ) >= 500 ){
    #baf_shift <- mean(tmp_baf$baf)
    baf_shift <- mclust::Mclust(data = tmp_baf$baf, G = c(1:3))
    baf_shift <- max(baf_shift$parameters$mean,na.rm = T)
    baf_scale <- 0.5/baf_shift
  }else{
    print( " Not enough SNP sites to estimate the BAF variance, skipping")
    baf_scale <- NA
  }
  return(baf_scale)

}



#' Correct BAF Values Using a Scaling Factor
#'
#' Applies a scaling factor to B-allele frequency (BAF) values greater than 0.4, and caps all BAF values at 0.5.
#'
#' @param baf A data frame or tibble with a numeric column named \code{baf}.
#' @param baf_scale Numeric. The scaling factor to be applied to BAF values greater than 0.4.
#'
#' @return A data frame or tibble with the corrected \code{baf} column.
#'
#' @importFrom dplyr mutate rowwise
#' @export
CorrectBAF <- function( baf , baf_scale){
  baf <- baf %>%
    dplyr::mutate( baf = ifelse( baf > 0.4, baf * baf_scale, baf ) )%>%
    dplyr::rowwise() %>%
    dplyr::mutate( baf = pmin( baf, 0.5 ) )
  return(baf)
}



#' Search for Breakpoints within a Segment Using BAF Data
#'
#' Refines breakpoints within a segment using B-allele frequency (BAF) data and AI binning/merging procedures.
#' If enough informative BAF sites are present, the segment may be split into finer regions based on BAF clustering and merging.
#'
#' @param seg_row A data frame row (list or tibble row) representing a single segment. Must have columns: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean, Segment_Mean_raw, Count, Baseline_cov, gatk_gender, pipeline_gender, size.
#' @param baf A data frame or tibble containing BAF data. Must include columns: Chromosome, Start, End, baf.
#' @param opt A list of options. Must include \code{aibinsize} (bin size for AI), \code{snpmin} (minimum SNP count), and \code{minaisize} (minimum AI size).
#'
#' @return A data frame with the refined segment(s), including updated breakpoints, BAF metrics, and a \code{BreakpointSource} column indicating whether breakpoints were post-processed or from GATK.
#'
#' @importFrom dplyr filter select mutate
#' @export
SearchBreakpoint <- function(seg_row, baf, opt ){
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
    MAF = NA,
    MAF_Probes = 0,
    MAF_gmm_G = NA,
    MAF_gmm_weight = NA,
    #MAF_gmm_variance = NA,
    size = seg_row$size)

  tmp_baf <- baf %>%
    dplyr::filter( Chromosome == seg_row$Chromosome & Start >= seg_row$Start & End <= seg_row$End ) %>%
    dplyr::filter( baf != 0 )
  n_baf <- nrow(tmp_baf)
  if( is.null(n_baf) ){ n_baf <- 0}
  if(n_baf > 0 ){
    binned_data <- binbaf(data = tmp_baf, bin_size = opt$aibinsize)
    max_nonzero_n <- max(binned_data$nonzero_count,na.rm = T)
    for ( min_prob in c(0: min(opt$snpmin, max_nonzero_n - 1 ) ) ){
      binned_data <- binned_data %>% filter(nonzero_count >= min_prob )
      binned_data <- CallMerge(data = binned_data,AIorSeg = "AI",opt = opt, tmp_baf = tmp_baf)
    }
    binned_data <- binned_data %>% filter( End - Start >= opt$minaisize)
    if( nrow(binned_data) > 0 ){
      merge_ai <- CallMerge(data = binned_data, AIorSeg = "AI", opt = opt, tmp_baf = tmp_baf)
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


  if(nrow(tmp_seg) > 1){

    tmp_seg$BreakpointSource <- "Postprocess"}else{tmp_seg$BreakpointSource <- "GATK"}

  return(tmp_seg)

}


#' Assign Quality Tag to a Segment Based on BAF GMM Metrics
#'
#' Assigns a quality tag ("PASS" or "FAILED") to a segment based on the BAF GMM weight, probe count, and number of clusters.
#'
#' @param MAF_gmm_weight Numeric. The mixture weight of the most prominent BAF GMM cluster.
#' @param MAF_Probes Integer. The number of nonzero BAF probes in the segment.
#' @param MAF_gmm_G Integer. Number of BAF GMM clusters.
#' @param snpmin Integer. The minimum SNP count (should be passed from \code{opt$snpmin}).
#'
#' @return Character. "PASS" if the segment passes quality checks, "FAILED" otherwise.
#'
#' @export
AddQualTag <- function(MAF_gmm_weight, MAF_Probes, MAF_gmm_G){
  if(!is.na(MAF_gmm_G)){
    if( MAF_Probes >= 2*opt$snpmin ){
      if( MAF_gmm_weight >= 0.35){Qualtag <- "PASS" }
      else{ Qualtag <- "FAILED" }


    }else{Qualtag <- "FAILED"}
  }else{Qualtag <- "FAILED"}
  return(Qualtag)
}


