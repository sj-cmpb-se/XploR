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
      #variances <- gmm_model$parameters$variance$sigmasq
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
#' @param maxgap Maximum gap size inside a bin. If exceed then start another bin.( default: 1000000)
#' @param snpnum SNP number in each bin.( default: 30 )
#' @param maxbinsize Maximum bin size.( default: 2000000 )
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
#'   \item{gmm_mean, gmm_weight, gmm_G}{Unnested maf GMM metrics.}
#'
#' @importFrom dplyr mutate group_by summarise ungroup filter n
#' @importFrom tidyr unnest_wider
#' @export
BinMaf <- function(data, datatype,
                   maxgap = 1000000, snpnum = 30,
                   maxbinsize = 2000000, minbinsize = 500000,
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
          bin_count <- i - start_index + 1
          current_span <- df$Pos[i] - df$Pos[start_index]

          # Check bin break conditions
          if (i < n) {
            if ( (bin_count == snpnum  || df$gap[i] > maxgap || current_span > maxbinsize) && current_span >= minbinsize ) {
              current_bin_id <- current_bin_id + 1
              start_index <- i + 1
            }
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
SearchBreakpoint <- function(seg_row, maf, pon_ref,
                             mergeai = 0.15,
                             snpmin = 7,
                             maxgap = 2000000, snpnum = 30,
                             maxbinsize = 2000000,minbinsize = 500000,
                             minsnpcov = 20,
                             segmethod = "cbs",
                             cbssmooth = "no"){
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
if( ! seg_row$Chromosome %in% c("X","Y")){
  tmp_maf <- maf %>%
    dplyr::filter( Chromosome == seg_row$Chromosome & Pos >= seg_row$Start & Pos <= seg_row$End ) %>%
    dplyr::filter( maf != 0 )
  n_maf <- nrow(tmp_maf)
  if( is.null(n_maf) ){ n_maf <- 0}
  if(n_maf > 0 ){
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
        binned_data <- binned_data %>% filter(nonzero_count >= min_prob )
        binned_data <- CallMerge(data = binned_data, AIorSeg = "AI", mergeai = mergeai, snpmin = snpmin, tmp_maf = tmp_maf)
      }
      binned_data <- binned_data %>% filter( End - Start >= 2*minbinsize)
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
    if(segmethod == "cbs"){
      # Prepare data for CBS (use gmm_mean as the segmentation track)
      maf_seg_data <- data.frame(
        chrom = binned_data$Chromosome,
        maploc = binned_data$Start,    # or bin center if you prefer
        maf = binned_data$gmm_mean
      )

      # Create CNA object
      maf_CNA <- DNAcopy::CNA(
        genomdat = maf_seg_data$maf,
        chrom = maf_seg_data$chrom,
        maploc = maf_seg_data$maploc,
        data.type = "logratio"
      )

      if( cbssmooth == "yes"){
        maf_CNA_smoothed <- DNAcopy::smooth.CNA(maf_CNA)
        maf_cbs <- DNAcopy::segment(maf_CNA_smoothed, verbose = 1)
      }else{
        maf_cbs <- DNAcopy::segment(maf_CNA, verbose = 1)
      }
      maf_cbs$segRows <- maf_cbs$segRows %>%
        dplyr::mutate( n_bins = endRow - startRow + 1,
                seg_id = dplyr::row_number())
      binned_data$merge_index <- rep( maf_cbs$segRows$seg_id, times = maf_cbs$segRows$n_bins  )
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
                MAF = gmm_mean,
                MAF_Probes = nonzero_count,
                MAF_gmm_G = gmm_G,
                MAF_gmm_weight = gmm_weight
        ) %>%
        dplyr::select( -gmm_mean, -nonzero_count, -gmm_G, -gmm_weight, -merge_index, -each_maf) %>%
        dplyr::relocate( Sample, Chromosome, Start, End, Num_Probes, Segment_Mean, gatk_SM_raw,
                  gatk_count, gatk_baselinecov, gatk_gender, pipeline_gender,
                  MAF, MAF_Probes, MAF_gmm_G, MAF_gmm_weight, size )

    }

  }
  tmp_seg<- CorrectBias( tmp_seg, pon_ref, tmp_maf)
  if(nrow(tmp_seg) > 1){
    tmp_seg$BreakpointSource <- "Postprocess"}else{tmp_seg$BreakpointSource <- "GATK"}
}else{
 tmp_seg["BreakpointSource"] <- "GATK"
}
  return(tmp_seg)
}


#' Assign Quality Tag to a Segment Based on maf GMM Metrics
#'
#' Assigns a quality tag ("PASS" or "FAILED") to a segment based on the maf GMM weight and probe count thresholds.
#'
#'A segment is marked as "PASS" if the maf GMM weight is at least 0.35 and the probe count is at least twice the minimum SNP count (\code{snpmin}). Otherwise, it is marked as "FAILED".
#' @param Chromosome Character.
#' @param MAF_gmm_weight Numeric. The mixture weight of the most prominent maf GMM cluster.
#' @param MAF_Probes Integer. The number of nonzero maf probes in the segment.
#' @param MAF_gmm_G Integer. Number of maf GMM clusters.
#' @param snpmin Integer. The minimum nonzero SNP count.
#'
#' @return Character. "PASS" if the segment passes quality checks, "FAILED" otherwise.
#'
#' @examples
#' AddQualTag(Chromosome = "1", MAF_gmm_weight = 0.4, MAF_Probes = 20, MAF_gmm_G = 2, snpmin = 7)
#' AddQualTag(Chromosome = "X", MAF_gmm_weight = 0.2, MAF_Probes = 10, MAF_gmm_G = 1, snpmin = 7)
#'
#' @export
AddQualTag <- function(Chromosome, MAF_gmm_weight, MAF_Probes, MAF_gmm_G, snpmin){

  if(!is.na(MAF_gmm_G)){
    if( MAF_Probes >= 2*snpmin ){
      if( MAF_gmm_weight >= 0.35){Qualtag <- "PASS" }
      else{ Qualtag <- "FAILED" }


    }else{Qualtag <- "FAILED"}
  }else{Qualtag <- "FAILED"}
 if(Chromosome %in% c("X","Y")){
   Qualtag <- "Exclude"
 }
  return(Qualtag)
}

#' Correct MAF Bias in Segments Using Panel of Normal Reference
#'
#' Adjusts the minor allele frequency (MAF) values in each segment for systematic bias using a panel of normal (PoN) reference. For each segment, compares the segment MAF to the PoN MAF distribution, and applies a logit-based correction if the segment MAF is not significantly different from the PoN. If the segment MAF is significantly different, retains the original value.
#'
#' @param tmp_seg Data frame or data.table. Segmented data to be corrected, must include columns: Chromosome, Start, End, Sample, Num_Probes, Segment_Mean, gatk_SM_raw, gatk_count, gatk_baselinecov, gatk_gender, pipeline_gender, MAF, MAF_Probes, MAF_gmm_G, MAF_gmm_weight, size.
#' @param pon_ref Data frame or data.table. Panel of normal reference, must include columns: Chromosome, Start, End, pon_mafs (comma-separated string of PoN MAF values).
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
CorrectBias <- function( tmp_seg, pon_ref, tmp_maf ){
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
              pipeline_gender, MAF, MAF_Probes, MAF_gmm_G, MAF_gmm_weight, size ) %>%
    dplyr::summarise( pon_mafs = paste( pon_mafs, collapse = ","))
  logit  <- function(p) qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  ilogit <- function(x) plogis(x)
  correct <- function( pon_mafs, each_mafs, gmm_mean){
    gmm_mean <- gmm_mean[1]
    if( !is.na(gmm_mean)){
    pon_mafs <- as.numeric(unlist(strsplit(pon_mafs, ",")))
    each_mafs <- as.numeric(unlist(strsplit(each_mafs, ";")))
    pon_median <- median(pon_mafs,na.rm = T)
    gmm_mean_corr <- NA
    if( length(pon_mafs) >= 10 && length(each_mafs) >= 5 ){
      w <- wilcox.test(pon_mafs, each_mafs )
      if( gmm_mean >= pon_median ||  (gmm_mean - pon_median <= 0.05 && w$p.value >= 0.01 ) ){
        y_centered = logit(gmm_mean) - logit( pon_median )
        gmm_mean_corr = pmin( ilogit(y_centered), 0.5 )
      }
    }else{
      if( !is.na(pon_median) ){
        if( gmm_mean - pon_median <= 0.1 ){
          y_centered = logit(gmm_mean) - logit( pon_median )
          gmm_mean_corr = pmin( ilogit(y_centered), 0.5 )
        }
      }
    }
    gmm_mean_corr <- ifelse( is.na(gmm_mean_corr), gmm_mean, gmm_mean_corr) }else{ gmm_mean_corr <- NA}
    return( gmm_mean_corr )
  }
  extract_maf <- function( tmp_maf, Start, End){

    each_mafs <- tmp_maf %>%
      dplyr::filter( Start <= Pos & End >= Pos ) %>%
      dplyr::filter( maf != 0 )

    each_mafs <- paste0( each_mafs$maf, collapse = ";")
    return(each_mafs)
  }

  tmp_seg_corr <- seg_corr %>%
    dplyr::mutate( each_mafs = extract_maf( tmp_maf, Start, End)) %>%
    dplyr::mutate( gmm_mean_corr = correct(pon_mafs, each_mafs, gmm_mean=MAF)) %>%
    dplyr::mutate( MAF = gmm_mean_corr) %>%
    dplyr::select( -gmm_mean_corr, -each_mafs, -pon_mafs ) %>%
    dplyr::relocate( MAF, .before = MAF_Probes)

  return(tmp_seg_corr)
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
