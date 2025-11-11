#' Correct Copy Ratio File if GATK and Pipeline Gender Disagree
#'
#' Adjusts the X chromosome log2 copy ratio in the coverage data if the gender assigned by GATK and the pipeline do not match.
#'
#' @param cov Data frame or tibble. Coverage data, must include columns \code{CONTIG} and \code{LOG2_COPY_RATIO}.
#' @param seg Data frame or tibble. Segmentation data, must include columns \code{gatk_gender} and \code{pipeline_gender}.
#'
#' @return The corrected coverage data frame or tibble.
#'
#' @details
#' - If GATK gender is "male" and pipeline gender is "female", subtracts 1 from X chromosome log2 copy ratio.
#' - If GATK gender is "female" and pipeline gender is "male", adds 1 to X chromosome log2 copy ratio.
#' - Otherwise, returns the input coverage unchanged.
#'
#' @importFrom dplyr mutate
#' @examples
#' # CorrectGender(cov, seg)
#'
#' @export
CorrectGender <- function(cov, seg){
  gatkgender <- unique(seg$gatk_gender)[1]
  pipelinegender <- unique(seg$pipeline_gender)[1]
  print( paste0("GATK gender is:", gatkgender))
  print( paste0("Real gender is:", pipelinegender))
  if(gatkgender == "male" && pipelinegender == "female"){
    cov <- cov %>%
      dplyr::mutate( LOG2_COPY_RATIO = ifelse( CONTIG == "X" , LOG2_COPY_RATIO - 1, LOG2_COPY_RATIO) )
    print( "Gender didn't match. GATK is male, Real gender is female, Correcting X chromosome CR.")
  }else if(gatkgender == "female" && pipelinegender == "male"){
    cov <- cov %>%
      dplyr::mutate( LOG2_COPY_RATIO = ifelse( CONTIG == "X" , LOG2_COPY_RATIO + 1, LOG2_COPY_RATIO) )
    print( "Gender didn't match. GATK is female, Real gender is male, Correcting X chromosome CR.")
  }else{

    print( "Gender is correct, no need to correct the CR file.")
  }

  return(cov)
}



#' Calculate CNF based on purity and diploid coverage size factor
#'
#' Adjusts the log2 copy ratio and computes the copy number fraction (CNF) for each region, given sample purity and scale factor.
#'
#' @param cr Data frame or tibble. Coverage data, must include columns \code{CONTIG} and \code{LOG2_COPY_RATIO}.
#' @param gender Character. Sample gender ("male" or "female").
#' @param purity Numeric. Tumor purity (fraction between 0 and 1).
#' @param sf Numeric. Scale factor for diploid coverage.
#'
#' @return The coverage data frame or tibble with additional columns \code{Segcov}, updated \code{LOG2_COPY_RATIO}, and \code{CNF}.
#'
#' @details
#' - Adjusts the log2 copy ratio using the scale factor.
#' - Calculates CNF based on purity and scale factor.
#'
#' @importFrom dplyr rowwise mutate
#' @examples
#' # CRCorrectPurity(cr, gender = "female", purity = 0.7, sf = 1)
#'
#' @export
CRCorrectPurity <- function(cr, gender, purity, sf ){

  # Calculate the CNF according to purity and cr from denoised file.
  cr <- cr %>%
    dplyr::rowwise() %>%
    dplyr::mutate( Segcov = SegmentMeanToOriCov( SM = LOG2_COPY_RATIO, gender = gender, chromosome = CONTIG, diploid_cov = 100 )) %>% ## change to segcov
    dplyr::mutate( LOG2_COPY_RATIO = CalSM( max_L_mu = sf,
                                     Segcov = Segcov,
                                     gender = gender,
                                     Chromosome = CONTIG) ) %>% ## correct the sf
    dplyr::mutate( CNF = ifelse( gender == "male" && CONTIG %in% c("X","Y") ,
                          (2 ^LOG2_COPY_RATIO - (1 - purity) ) / purity,
                          (2 * 2^LOG2_COPY_RATIO - 2 * (1 - purity)) / purity ) ) # correct the purity
  return(cr)
}


#' Bin and Smooth Copy Number Data from Coverage File
#'
#' Splits the coverage data by chromosome, bins into 300bp windows, and computes smoothed copy number values using a larger window size.
#'
#' @param cov Data frame or tibble. Must include columns \code{CONTIG}, \code{START}, \code{END}, and \code{CNF}.
#' @param cov_binsize Numeric. Bin size (in bases) for smoothing (e.g., 10000 for 10kb).
#'
#' @return A data frame with columns: \code{CONTIG}, \code{bin_start}, \code{bin_norm_id}, and smoothed CNF values.
#'
#' @details
#' - Bins are generated every 300bp within each region.
#' - Smoothed copy number is computed as the median CNF within each larger (\code{cov_binsize}) window.
#'
#' @importFrom dplyr mutate group_by summarize select if_else
#' @importFrom tidyr unnest
#' @importFrom purrr map2
#' @examples
#' # CovBin(cov, cov_binsize = 10000)
#'
#' @export
CovBin <- function(cov,cov_binsize){
  # from cov file generate bin file with window 300bp
  # smooth the bin file to cov_binsize window.
  cov_bin_chr <- split(cov, f = "CONTIG")
  sep_bin <- lapply(cov_bin_chr, function(chr_cov){

    chr_cov <- chr_cov %>%
      dplyr::mutate(bin = map2( START, END, ~ seq(.x, .y, by = 300))) %>%
      tidyr::unnest(bin) %>%
      dplyr::mutate(bin_start = floor((bin - 1) / 300) * 300 + 1) %>%
      dplyr::mutate(bin_end = bin_start + 299) %>%
      dplyr::mutate(bin_end = if_else(bin_end > END, END, bin_end),
             bin_start = if_else(bin_start < START + 1 , START + 1 , bin_start)) %>%
      dplyr::mutate(bin_norm_id = ceiling(bin_start/cov_binsize)*cov_binsize ) %>%
      dplyr::group_by( CONTIG , bin_norm_id ) %>%
      dplyr::summarize(smoothed_cnf = median(CNF, na.rm = TRUE)) %>%
      dplyr::mutate( bin_start = bin_norm_id - cov_binsize) %>%
      dplyr::select( CONTIG, bin_start, bin_norm_id, everything() ) %>%
      dplyr::mutate( smoothed_bin_cnf = round(smoothed_cnf * 10 )/ 10 )
    colnames(chr_cov)[3] <- "bin_end"

    return(chr_cov)
  } )
  cov_bin <- do.call(rbind,sep_bin)
  return(cov_bin)

}


#' Round BAF Values into 20 Bins
#'
#' For each chromosome, rounds the B-allele frequency (BAF, or \code{af}) values into 20 bins (increments of 0.05).
#'
#' @param ai A list of data frames or tibbles, one per chromosome. Each must include columns \code{af}, \code{contig}, \code{start}, \code{stop}, \code{allele1Count}, \code{allele2Count}.
#'
#' @return A data frame with columns: \code{contig}, \code{start}, \code{stop}, \code{allele1Count}, \code{allele2Count}, \code{norm_af}.
#'
#' @details
#' - BAF values (\code{af}) are rounded to the nearest 0.05.
#' - Output is a single data frame combining all chromosomes.
#'
#' @importFrom dplyr mutate select
#' @examples
#' # ai_list <- list(chr1 = data.frame(af = runif(10), contig = "1", start = 1:10, stop = 11:20, allele1Count = 1:10, allele2Count = 11:20))
#' # RoundAI(ai_list)
#'
#' @export
RoundAI <- function(ai){
  ## round the baf value into 20bin.
  round_chr_ai <- lapply(ai,function(chr_ai){

    chr_ai <- chr_ai %>%
      dplyr::mutate( norm_af = round(af*20)/20) %>%
      dplyr::select( CONTIG, POSITION, REF_COUNT, ALT_COUNT, norm_af )
    return(chr_ai)
  })

  round_ai <- do.call(rbind, round_chr_ai)
  round_ai <- round_ai
  return(round_ai)

}


#' Keep Only Coverage Bins in Whitelisted Regions
#'
#' Filters smoothed coverage bins to include only those overlapping whitelisted regions.
#'
#' @param smooth_cov Data frame or tibble. Smoothed coverage data, must include columns \code{CONTIG}, \code{bin_start}, and \code{bin_end}.
#' @param whitelist Data frame or tibble. Whitelist regions, must include columns \code{chrom}, \code{start}, and \code{end}.
#' @param gender Character. femle or male.
#'
#' @return A data frame of smoothed coverage bins that overlap whitelisted regions.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' # KeepWhitelistCov(smooth_cov, whitelist)
#'
#' @export
KeepWhitelistCov <- function(smooth_cov, whitelist, gender){

  # keep bins in whitelist only
  smooth_cov_range <- GenomicRanges::makeGRangesFromDataFrame(
    smooth_cov,
    seqnames.field = "CONTIG",
    start.field = "bin_start",
    end.field = "bin_end",
    keep.extra.columns = TRUE )

  whitelist_range <- GenomicRanges::makeGRangesFromDataFrame(
    whitelist,
    seqnames.field = "chrom",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = FALSE

  )

  subsetted_smooth_cov_range <- as.data.frame(subsetByOverlaps(smooth_cov_range, whitelist_range))


  return(subsetted_smooth_cov_range)
}

#' Keep Only AI Bins in Whitelisted Regions
#'
#' Filters AI (allelic imbalance) bins to include only those overlapping whitelisted regions.
#'
#' @param rounded_ai Data frame or tibble. Rounded AI data, must include columns \code{contig}, \code{start}, and \code{stop}.
#' @param whitelist Data frame or tibble. Whitelist regions, must include columns \code{seqnames}, \code{start}, and \code{end}.
#'
#' @return A data frame of AI bins that overlap whitelisted regions.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom IRanges subsetByOverlaps
#' @examples
#' # KeepWhitelistAI(rounded_ai, whitelist)
#'
#' @export
KeepWhitelistAI <- function(rounded_ai, whitelist){
  # keep bins in whitelist only
  rounded_ai_range <- GenomicRanges::makeGRangesFromDataFrame(
    rounded_ai,
    seqnames.field = "contig",
    start.field = "start",
    end.field = "stop",
    keep.extra.columns = TRUE

  )

  whitelist_range <- GenomicRanges::makeGRangesFromDataFrame(
    whitelist,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    keep.extra.columns = FALSE

  )

  subsetted_rounded_ai_range <- as.data.frame(GenomicRanges::subsetByOverlaps(rounded_ai_range, whitelist_range))

  return(subsetted_rounded_ai_range)
}


#' Clean Homozygous AI Bins Based on LOH Segments
#'
#' Removes homozygous allelic imbalance (AI) bins that do not overlap with LOH (loss of heterozygosity) segments.
#'
#' @param rounded_ai Data frame or tibble. AI data, must include columns \code{contig}, \code{start}, \code{stop}, \code{allele1Count}, \code{allele2Count}.
#' @param call_seg Data frame or tibble. Segment calls, must include columns \code{chrom}, \code{loc.start}, \code{loc.end}, and \code{MAF}.
#'
#' @return A data frame of AI bins, with homozygous bins outside LOH regions removed.
#'
#' @importFrom dplyr filter select
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @importFrom S4Vectors subjectHits queryHits
#' @examples
#' # CleanHomalt(rounded_ai, call_seg)
#'
#' @export
CleanHomalt <- function(rounded_ai, call_seg){
  # Clean homozygous AI bins
  call_seg <- call_seg %>% dplyr::filter(MAF > 0.7 | MAF < 0.3)
  if( nrow(call_seg) > 0 ){
    loh_range <- GenomicRanges::makeGRangesFromDataFrame(
      call_seg,
      seqnames.field = "chrom",
      start.field = "loc.start",
      end.field = "loc.end",
      keep.extra.columns = FALSE

    )

    ai_range <- GenomicRanges::makeGRangesFromDataFrame(
      rounded_ai,
      seqnames.field = "CONTIG",
      start.field = "POSITION",
      end.field = "POSITION",
      keep.extra.columns = TRUE

    )

    overlaps <- GenomicRanges::findOverlaps(ai_range, loh_range)
    overlapping_loh <- loh_range[S4Vectors::subjectHits(overlaps)]
    overlapping_loh_df <- as.data.frame(overlapping_loh)
    rounded_ai[S4Vectors::queryHits(overlaps),c( "matched_seqnames", "matched_start", "matched_end" )] <-  overlapping_loh_df[,c(1:3)]
    final_ai <- rounded_ai %>%
      #filter( (allele1Count > 1 & allele2Count > 1 ) | seqnames %in% c("X","Y") ) %>%
      filter( ( ! ( is.na(matched_seqnames) & ( REF_COUNT < 2 | ALT_COUNT < 2  ) ) ) ) %>%
      select( -matched_seqnames, -matched_start, -matched_end )
  }else{ final_ai <- rounded_ai
  }
  colnames(final_ai)[1] <- "seqnames"

  return(final_ai)
}


#' Smooth AI Bin Values by Downsampling
#'
#' If the length of the input vector is greater than 30, returns a vector of up to 30 values sampled according to the observed frequency of each unique value. Otherwise, returns the vector unchanged.
#'
#' @param baf Numeric or character vector. Typically a vector of rounded AI (e.g., BAF) values.
#'
#' @return A numeric vector, either unchanged (if length <= 30) or of length <= 30, with values sampled according to observed frequencies.
#'
#' @examples
#' set.seed(1)
#' x <- sample(c(0.1, 0.2, 0.3, 0.4), 100, replace = TRUE)
#' AIbinSmooth(baf = x)
#'
#' @export
AIbinSmooth <- function(baf) {
  # Return a vector of same frequency if the length is higher than 30
  baf <- baf[!is.na(baf)]
  if( length(baf) > 0){
    counts <- data.frame(table(baf))
    mode_value <- counts[order(counts$Freq,decreasing = T),]
    colnames(mode_value) <- c("value","Freq" )
    total_dots <- length(baf)
    mode_value$prop <- mode_value$Freq/total_dots
    if( total_dots > 30 ){

      re <- rep(mode_value$value,round(mode_value$prop * 30))

    }else{re <- baf }

    re <- as.numeric(as.character(re))
  }else{re <- NA}

  return(re)
}



#' Smooth AI Bins by Window Size and Downsampling
#'
#' Smooths allelic imbalance (AI) bins by a specified window size, keeping up to 30 data points per bin using frequency-based downsampling.
#'
#' @param df Data frame or tibble. Must include columns \code{seqnames}, \code{start}, and \code{norm_af}.
#' @param ai_binsize Numeric. Bin size (in bases) for smoothing.
#' @param gender Character. Sample gender ("male" or "female"). (Currently not used in this function, but included for compatibility.)
#'
#' @return A data frame with columns: \code{seqnames}, \code{bin_start}, \code{bin_end}, and smoothed AI values.
#'
#' @details
#' - Splits the data by chromosome.
#' - For each bin, keeps up to 30 values using \code{\link{ai_bin_smooth}}.
#'
#' @importFrom dplyr mutate group_by summarize ungroup select
#' @importFrom tidyr unnest
#' @examples
#' # SmoothAI(df, ai_binsize = 10000, gender = "female")
#'
#' @export
SmoothAI <- function(df, ai_binsize , gender){
  # Smooth AI bins by ai_binsize, for each bin keep 30 datapoints.
  ai_bin_chr <- split(df, f = df$seqnames)
  ai_bin_chr <- Filter(function(x) nrow(x) > 0, ai_bin_chr )
  sep_bin <- lapply(c(1:length(ai_bin_chr)), function(aa){

    chr_ai_smoothed <- ai_bin_chr[[aa]] %>%
      dplyr::mutate(bin_norm_id = ceiling(POSITION/ai_binsize)*ai_binsize ) %>%
      dplyr::group_by( seqnames , bin_norm_id ) %>%
      dplyr::summarize(smoothed_ai = list( AIbinSmooth( norm_af ) ))   %>%
      tidyr::unnest( smoothed_ai ) %>%
      dplyr::mutate( bin_start = bin_norm_id - ai_binsize) %>%
      dplyr::select( seqnames, bin_start, bin_norm_id, everything() ) %>%
      dplyr::ungroup()

    colnames( chr_ai_smoothed )[3] <- "bin_end"
    return(chr_ai_smoothed)
  })

  ai_bin_smoothed <- do.call(rbind,sep_bin)
  return(ai_bin_smoothed)
}


#' Multi-panel Plot of Copy Number, BAF, Clonality, and Quality Across Chromosomes
#'
#' Produces a summary plot with four panels: copy number, BAF, clonality, and segment quality, faceted by chromosome.
#'
#' @param df_cov Data frame or tibble. Smoothed coverage data, must include columns \code{seqnames}, \code{start}, \code{end}, and \code{smoothed_bin_cnf}.
#' @param df_ai Data frame or tibble. Smoothed AI data, must include columns \code{seqnames}, \code{bin_start}, and \code{smoothed_ai}.
#' @param call_seg Data frame or tibble. Segment calls, must include columns \code{seqnames}, \code{loc.start}, \code{loc.end}, \code{CN}, \code{MAF}, \code{ccf_final}, \code{FILTER}.
#' @param whitelist Data frame or tibble. Whitelist regions, must include columns \code{seqnames}, \code{start}, and \code{end}.
#' @param gender Character. Sample gender ("male" or "female").
#' @param prefix Character. Plot title or output prefix.
#' @param purity Numeric. Tumor purity (not directly used in the plot, but included for compatibility).
#'
#' @return A \code{grob} object (from \code{gridExtra::grid.arrange}) representing the combined plot.
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_segment geom_point scale_fill_manual scale_color_manual facet_grid theme_minimal scale_y_continuous labs theme element_blank element_text element_line margin vars
#' @importFrom dplyr mutate filter select
#' @importFrom scales label_number
#' @importFrom gridExtra grid.arrange
#' @examples
#' # PlotCov(df_cov, df_ai, call_seg, whitelist, gender = "female", prefix = "Sample1", purity = 0.7)
#'
#' @export
PlotCov <- function(df_cov, df_ai, call_seg , whitelist , gender, prefix, purity ){
  colnames(whitelist)[1] <- "seqnames"
  ## Color code
  if(gender == "male"){
    color <- c("darkblue","darkgreen","darkred","darkorchid4","darkgoldenrod" )
    color <- rep(color,length.out =24)
    chrom_levels <- c(seq(1:22),"X","Y")
    names(color) <- chrom_levels
  }else{
    color <- c("darkblue","darkgreen","darkred","darkorchid4","darkgoldenrod" )
    color <- rep(color,length.out =23)
    chrom_levels <- c(seq(1:22),"X")
    names(color) <- chrom_levels
  }
  df_cov$seqnames <- factor(df_cov$seqnames,levels = chrom_levels)
  df_ai$seqnames <- factor(df_ai$seqnames, levels = chrom_levels)
  call_seg$seqnames <- factor(call_seg$seqnames, levels = chrom_levels)
  whitelist$seqnames <- factor(whitelist$seqnames, levels = chrom_levels)
  y_lim <- call_seg %>%
    dplyr::mutate(size = loc.end - loc.start) %>%
    dplyr::filter(size > 10000000)

  if(nrow(y_lim) > 10 ){
    y_lim <- max(y_lim$CN) + 2
    y_lim <- ifelse(y_lim >= 16 , 16, y_lim)
  }else{
    y_lim <- call_seg %>% dplyr::mutate(size = loc.end - loc.start) %>% dplyr::filter(size > 1000000)
    if(nrow(y_lim) > 10){
      y_lim <- max(y_lim$CN) + 2
      y_lim <- ifelse(y_lim >= 16 , 16, y_lim)
    }else{ y_lim <- call_seg %>% dplyr::mutate( size = loc.end - loc.start) %>% dplyr::filter( size > 500000)
    y_lim <- median(y_lim$CN) + 4
    y_lim <- ifelse( y_lim > 16, 16, y_lim) }
  }

  y_lim <- ifelse(y_lim < 4 , 4, y_lim)



  if( y_lim >=10){
    p_margin <- ggplot2::margin(t = 1, r = 1, b = 0.5, l = 2, unit = "pt")
    q_margin <- ggplot2::margin(t = 0, r = 1, b = 0.5, l = 7, unit = "pt")
    line_pos <- c(4,6,8)
  }else{
    p_margin <- ggplot2::margin(t = 1, r = 1, b = 0.5, l = 2, unit = "pt")
    q_margin <- ggplot2::margin(t = 0, r = 1, b = 0.5, l = 2, unit = "pt")
    line_pos <- c(1,4)
  }

  df_cov <- df_cov %>% dplyr::mutate(smoothed_bin_cnf = ifelse( smoothed_bin_cnf >= y_lim, y_lim, smoothed_bin_cnf ))


  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 2,color="black",linewidth=1) +
    ggplot2::geom_hline(yintercept = line_pos,color="grey",linewidth=0.5,linetype = "dashed") +
    ggplot2::geom_segment(data = whitelist, ggplot2::aes(x = start-100000, xend = end+100000,y = y_lim,yend = y_lim,color = seqnames ), linewidth = 2) +
    ggplot2::geom_segment(data = df_cov,
                 ggplot2::aes(x = start, xend = end, y = smoothed_bin_cnf - 0.1, yend = smoothed_bin_cnf + 0.1, color = seqnames),
                 alpha = 0.2, linewidth = 1 ) +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::geom_segment( data = call_seg, ggplot2::aes(x = loc.start, xend = loc.end,y = CN, yend = CN),color="cyan3",linewidth = 1 ) +
    ggplot2::facet_grid(cols = vars(seqnames), scales = "free_x", space = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(limits = c(0, y_lim), breaks = seq(0,y_lim), labels = label_number(accuracy = 0.01) )  +
    ggplot2::labs( title = prefix, y = "Copy Number") +
    ggplot2::theme(legend.position = "none",
          panel.spacing = grid::unit(0, "lines"),
          strip.background = ggplot2::element_blank(),
          strip.text = ggplot2::element_text(face = "bold",size = 10),
          plot.title = ggplot2::element_text(hjust = 0.5),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.line.y.left = ggplot2::element_line(color = "black"),
          axis.ticks.y = ggplot2::element_line(color = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          plot.margin = p_margin
    )


  # ggsave(filename = "test.png",plot = p,device = "png",path = "/Volumes/mjin47/CNV/DRAGEN_RUN/annotation",width = 15,height = 5,dpi = 900, bg = 'white')



  q <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = 0.5,color="black",linewidth=0.7) +
    ggplot2::geom_hline(yintercept = c(0,1),color="grey",linewidth=0.5,linetype="dashed") +
    ggplot2::geom_segment(data=df_ai, ggplot2::aes(x = bin_start,
                                 xend = bin_end,
                                 y= smoothed_ai - 0.05,
                                 yend = smoothed_ai + 0.05 ,
                                 color = seqnames ),alpha = 0.2,linewidth = 1) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::geom_segment(data = call_seg,
                 ggplot2::aes(x = loc.start, xend = loc.end, y= MAF, yend = MAF),
                 linewidth = 1,color="cyan3") +
    ggplot2::geom_segment( data = call_seg,
                  ggplot2::aes(x = loc.start, xend = loc.end, y = 1-MAF, yend = 1-MAF),
                  linewidth = 1, color = "cyan3") +
    ggplot2::facet_grid(cols = vars(seqnames), scales = "free_x", space = "free_x")+
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25)) +
    ggplot2::theme_minimal() +
    ggplot2::labs( y = "BAF" ) +
    ggplot2::theme(legend.position = "none",
          panel.spacing = grid::unit(0, "lines"),
          strip.text = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_line(color = "black"),
          axis.line.y.left = ggplot2::element_line(color = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          plot.margin =  q_margin
    )
  call_seg <- call_seg %>%
    dplyr::mutate( ccf_final = ifelse(is.na(ccf_final),1,ccf_final))
  clonality <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = c(0,1),color="black",linewidth=0.5, linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0.8,color="grey40",linewidth=0.5, linetype = "dashed") +
    ggplot2::geom_segment(data = call_seg ,
                 ggplot2::aes(x = loc.start, xend = loc.end,y = ccf_final,yend = ccf_final,color = "grey40"),
                 linewidth = 2) +
    ggplot2::geom_segment(data = whitelist,
                 ggplot2::aes(x= start - 100000, xend = end + 100000, y = 0,yend = 0, color = seqnames, fill = seqnames),
                 linewidth = 2) +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::scale_color_manual(values = color) +
    ggplot2::facet_grid(cols = vars(seqnames), scales = "free_x", space = "free_x") +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0,1), labels = label_number(accuracy = 0.01) )  +
    ggplot2::labs( title = prefix, y = "Clonality") +
    ggplot2::theme(legend.position = "none",
          panel.spacing = grid::unit(0, "lines"),
          strip.text = ggplot2::element_blank(),
          strip.background = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.y = ggplot2::element_line(color = "black"),
          axis.line.y.left = ggplot2::element_line(color = "black"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          plot.margin =  q_margin,
          plot.title = ggplot2::element_blank()
    )

  quality <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = call_seg, ggplot2::aes(x = loc.start, xend = loc.end, y= 0.1, yend = 0.1, color = FILTER ), linewidth = 4) +
    ggplot2::scale_color_manual(values = c("PASS" = "grey15", "FAILED" = "grey" , "REVIEW" = "grey35")) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0,1), labels = c("QC","") )+
    ggplot2::facet_grid(cols = vars(seqnames), scales = "free_x", space = "free_x")+
    ggplot2::theme_minimal() +
    ggplot2::labs( y = "" ) +
    ggplot2::theme(legend.position = "none",
          strip.text = ggplot2::element_blank(),
          panel.spacing = grid::unit(0, "lines"),
          strip.background = ggplot2::element_blank(),
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          plot.margin = grid::unit(c(0, 1, 0.5, 2.2), "pt")
    )


  summary_plot <- gridExtra::grid.arrange(p,q,clonality, quality,ncol=1, heights= c(5,3,1.5,0.5))
  return(summary_plot)

}
