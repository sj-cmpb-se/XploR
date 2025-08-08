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
#'
#' @return A data frame with columns \code{Chromosome}, \code{Pos}, and \code{maf}, containing only autosomal positions (\code{Chromosome} not X or Y), with minor allele frequency (\code{maf}) between 0 and 0.5.
#'
#' @details
#' For \code{aitype = "gatk"}, the function reads the file, filters for loci with at least 10 total counts, calculates the B-allele frequency (\code{baf}), and computes the minor allele frequency (\code{maf}).
#'
#' For \code{aitype = "dragen"}, the function extracts the first, second, and fourth columns, and renames them as \code{Chromosome}, \code{Pos}, and \code{baf}.
#'
#' For \code{aitype = "other"}, the function expects columns \code{Chromosome}, \code{Pos}, and \code{baf}.
#'
#' @examples
#' \dontrun{
#' # Read AI data from a count file
#' maf_df <- ReadAI("gatk", "sample.allele_counts")
#'
#' # Read AI data from a DRAGEN file
#' maf_df <- ReadAI("dragen", "sample_dragen.tumor.ballele.counts.gz")
#'
#' # Read AI data from a BAF file
#' maf_df <- ReadAI("other", "sample.tsv")
#' }
#'
#' @importFrom dplyr filter mutate select
#' @export
ReadAI <- function(aitype, ballele){
  if( aitype == "gatk"){
    filtered_lines <- grep("^[^@]", readLines(ballele), value = TRUE)
    baf <- read.table(text = paste(filtered_lines, collapse = "\n"), sep = "\t", header = TRUE)
    baf <- baf %>%
      dplyr::filter( REF_COUNT + ALT_COUNT >= 10) %>%
      dplyr::mutate( baf = ALT_COUNT/(ALT_COUNT + REF_COUNT)) %>%
      dplyr::select( -ALT_COUNT, -REF_COUNT, -REF_NUCLEOTIDE, -ALT_NUCLEOTIDE )
  }else if( aitype == "dragen" ){
    baf <- read.table(ballele, sep = "\t", header = TRUE,stringsAsFactors = F)
    baf <- baf %>%
      dplyr::filter( allele1Count + allele2Count >= 10 ) %>%
      dplyr::mutate( baf = allele2Count/(allele2Count + allele1Count))
    baf <- baf[,c("contig","start","baf")]
  }else{ baf <- read.table(ballele, sep = "\t", header = TRUE,stringsAsFactors = F)}
  colnames(baf) <- c("Chromosome", "Pos", "baf")
  maf <- baf %>%
    dplyr::filter(!Chromosome %in% c("X", "Y")) %>%
    dplyr::filter(baf != 0 & baf != 1) %>%
    dplyr::mutate(maf = ifelse(baf > 0.5, 1 - baf, baf)) %>%
    dplyr::select( -baf )
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

  n_mafs <- length(maf_values)
  sd <- sd( maf_values, na.rm = T )
  if( n_mafs >= 3 && sd >= 1e-6 ){
    gmm_model <- mclust::Mclust(maf_values, G = 1:10)
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

    Estimated_maf <- c(
      "gmm_mean" = gmm_mean,
      #"gmm_variance" = gmm_variance,
      "gmm_weight" = gmm_weight,
      "gmm_G" = num_components
    )}else{
      gmm_mean <- mean(maf_values,na.rm = T)
      gmm_mean <- ifelse( gmm_mean > 0.48, 0.5, gmm_mean )
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
#' Bins input MAF data into fixed-size windows and calculates summary metrics (e.g., estimated MAF, nonzero counts) for each bin.
#'
#' @param data A data frame or tibble containing at least \code{Chromosome}, \code{Start}, and \code{maf} columns.
#' @param bin_size Numeric. The bin size to use (e.g., 1e6 for 1Mb bins).
#'
#' @return A tibble with one row per bin and columns:
#'   \item{Chromosome}{Chromosome identifier.}
#'   \item{bin}{Bin start coordinate.}
#'   \item{Start}{Minimum \code{End} value in the bin.}
#'   \item{End}{Maximum \code{End} value in the bin.}
#'   \item{Estimated_maf}{List column with GMM-based maf estimation results.}
#'   \item{nonzero_count}{Number of nonzero maf values in the bin.}
#'   \item{each_maf}{Semicolon-separated string of nonzero maf values in the bin.}
#'   \item{gmm_mean, gmm_weight, gmm_G}{Unnested maf GMM metrics.}
#'
#' @importFrom dplyr mutate group_by summarise ungroup filter
#' @importFrom tidyr unnest_wider
#' @export
BinMaf <- function(data, bin_size){

  # Bin the data and calculate the desired metrics
  binned_data <- data %>%
    dplyr::mutate(bin = floor(Pos / bin_size) * bin_size) %>%
    dplyr::group_by(Chromosome, bin) %>%
    dplyr::summarise(
      Start = min(Pos),
      End = max(Pos),
      Estimated_maf = list(EstimateMAFbyGMM(maf_values = maf[maf != 0] )),
      nonzero_count = sum(maf != 0),
      each_maf = paste(maf[maf != 0 ], collapse = ";"),
      .groups = "drop"
    ) %>% tidyr::unnest_wider(col = "Estimated_maf")%>%
    dplyr::ungroup() %>%
    dplyr::filter( !is.na( gmm_mean))

  return(binned_data)

}


#' Estimate MAF Shift Scaling Factor
#'
#' Estimates a scaling factor to correct MAF shift.
#'
#' @param maf A data frame or tibble with at least a \code{maf} column (numeric MAF values).
#'
#' @return Numeric. The estimated MAF scaling factor, or \code{NA} if there are not enough SNPs.
#'
#' @details
#' Uses \code{\link[mclust]{Mclust}} for Gaussian mixture modeling.
#'
#' @importFrom dplyr filter
#' @importFrom mclust Mclust
#' @export
EstimateMAFshift <- function( maf){

  tmp_maf <- maf %>%
    dplyr::filter( maf > 0.4)
  if( nrow( tmp_maf ) >= 500 ){
    #maf_shift <- mean(tmp_maf$maf)
    maf_shift <- mclust::Mclust(data = tmp_maf$maf, G = c(1:3))
    maf_shift <- max(maf_shift$parameters$mean,na.rm = T)
    maf_scale <- 0.5/maf_shift
  }else{
    print( " Not enough SNP sites to estimate the maf variance, skipping")
    maf_scale <- NA
  }
  return(maf_scale)

}



#' Correct MAF Values Using a Scaling Factor
#'
#' Applies a scaling factor to MAF values.
#'
#' @param maf A data frame or tibble with a numeric column named \code{maf}.
#' @param maf_scale Numeric. The scaling factor to be applied to MAF values.
#'
#' @return A data frame or tibble with the corrected \code{maf} column.
#'
#' @importFrom dplyr mutate rowwise
#' @export
CorrectMAF <- function( maf , maf_scale){
  maf <- maf %>%
    dplyr::mutate( maf = ifelse( maf > 0.4, maf * maf_scale, maf ) )%>%
    dplyr::rowwise() %>%
    dplyr::mutate( maf = pmin( maf, 0.5 ) )
  return(maf)
}



#' Search for Breakpoints within a Segment Using MAF Data
#'
#' Refines breakpoints within a segment using MAF data.
#' If enough informative MAF sites are present, the segment may be split into finer regions based on MAF clustering and merging.
#'
#' @param seg_row A data frame row (list or tibble row) representing a single segment. Must have columns: Sample, Chromosome, Start, End, Num_Probes, Segment_Mean, Segment_Mean_raw, Count, Baseline_cov, gatk_gender, pipeline_gender, size.
#' @param maf A data frame or tibble containing MAF data. Must include columns: Chromosome, Start, End, maf.
#' @param aibinsize Numeric, AI bin size (default: 500000).
#' @param snpmin Numeric. Minimum SNP count required for a segment to be considered as a separate segment.
#' @param minaisize Numeric, smallest AI segment size (default: 1000000).
#' @param mergeai Numeric. Threshold for the difference in MAF (gmm_mean) between adjacent segments to allow merging.
#'
#'
#' @return A data frame with the refined segment(s), including updated breakpoints, MAF metrics, and a \code{BreakpointSource} column indicating whether breakpoints were post-processed or from GATK.
#'
#' @importFrom dplyr filter select mutate
#' @export
SearchBreakpoint <- function(seg_row, maf, mergeai, snpmin, minaisize, aibinsize ){
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

  tmp_maf <- maf %>%
    dplyr::filter( Chromosome == seg_row$Chromosome & Pos >= seg_row$Start & Pos <= seg_row$End ) %>%
    dplyr::filter( maf != 0 )
  n_maf <- nrow(tmp_maf)
  if( is.null(n_maf) ){ n_maf <- 0}
  if(n_maf > 0 ){
    binned_data <- BinMaf(data = tmp_maf, bin_size = aibinsize)
    max_nonzero_n <- max(binned_data$nonzero_count,na.rm = T)
    for ( min_prob in c(0: min( snpmin, max_nonzero_n - 1 ) ) ){
      binned_data <- binned_data %>% filter(nonzero_count >= min_prob )
      binned_data <- CallMerge(data = binned_data, AIorSeg = "AI", mergeai = mergeai, snpmin = snpmin, tmp_maf = tmp_maf)
    }
    binned_data <- binned_data %>% filter( End - Start >= minaisize)
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


#' Assign Quality Tag to a Segment Based on maf GMM Metrics
#'
#' Assigns a quality tag ("PASS" or "FAILED") to a segment based on the maf GMM weight and probe count thresholds.
#'
#'A segment is marked as "PASS" if the maf GMM weight is at least 0.35 and the probe count is at least twice the minimum SNP count (\code{snpmin}). Otherwise, it is marked as "FAILED".
#'
#' @param MAF_gmm_weight Numeric. The mixture weight of the most prominent maf GMM cluster.
#' @param MAF_Probes Integer. The number of nonzero maf probes in the segment.
#' @param MAF_gmm_G Integer. Number of maf GMM clusters.
#' @param snpmin Integer. The minimum nonzero SNP count.
#'
#' @return Character. "PASS" if the segment passes quality checks, "FAILED" otherwise.
#'
#' @examples
#' AddQualTag(MAF_gmm_weight = 0.4, MAF_Probes = 20, MAF_gmm_G = 2, snpmin = 7)
#' AddQualTag(MAF_gmm_weight = 0.2, MAF_Probes = 10, MAF_gmm_G = 1, snpmin = 7)
#'
#' @export
AddQualTag <- function(MAF_gmm_weight, MAF_Probes, MAF_gmm_G, snpmin){
  if(!is.na(MAF_gmm_G)){
    if( MAF_Probes >= 2*snpmin ){
      if( MAF_gmm_weight >= 0.35){Qualtag <- "PASS" }
      else{ Qualtag <- "FAILED" }


    }else{Qualtag <- "FAILED"}
  }else{Qualtag <- "FAILED"}
  return(Qualtag)
}


