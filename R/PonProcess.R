#' Generate Whitelist, Blacklist, and Cytoband-Detectable Edge Files from GATK PoN HDF5
#'
#' Processes a GATK Panel of Normals (PoN) HDF5 file to extract final intervals (whitelist) and excluded regions (blacklist), and outputs them as BED files. Also creates a cytoband-detectable edge file by overlapping whitelist intervals with cytoband annotations to define p and q arm coverage.
#'
#' This function calls \code{\link{FindCytoband}}, which must be available in the package namespace.
#'
#' @param pon_file Character. Path to the GATK PoN HDF5 file.
#' @param blacklist_bed Character. Output file path for the blacklist BED file (intervals excluded from the panel).
#' @param whitelist_bed Character. Output file path for the whitelist BED file (final intervals retained in the panel).
#' @param detectable_edge Character. Output file path for the cytoband-detectable edge file, used in downstream annotation.
#' @param gender Character. female or male. Gender of the GATK PON file. Mixed gender are not supported.
#' @param cytoband Character. Path to the cytoband annotation file with columns: chrom, chromStart, chromEnd, name, gieStain.
#'
#' @return Invisibly returns a list containing data frames: \code{blacklist}, \code{whitelist}, and \code{cytoband_edge}.
#'
#' @details
#' This function performs the following:
#' \itemize{
#'   \item Reads original and final intervals from the PoN file and maps contig indices to chromosome names.
#'   \item Writes a whitelist BED file with intervals retained for CNV calling.
#'   \item Writes a blacklist BED file with intervals filtered out by GATK PoN processing.
#'   \item Generates a cytoband-detectable edge file by intersecting whitelist intervals with cytoband annotations, summarizing detectable p and q arm boundaries per chromosome.
#' }
#'
#' @seealso \code{\link{FindCytoband}}
#'
#' @importFrom hdf5r H5File
#' @importFrom dplyr anti_join filter group_by summarize arrange mutate bind_rows rowwise ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom stringr str_replace
#' @importFrom utils write.table
#'
#' @examples
#' PonProcess("gatk.pon.hdf5", "pon_blacklist.bed", "pon_whitelist.bed", "detectable_edge.bed", "cytoband.txt")
#'
#' @export
PonProcess <- function(pon_file, blacklist_bed, whitelist_bed, cytoband, detectable_edge, gender ) {
  f <- hdf5r::H5File$new(pon_file, mode = "r")
  ori_contig_names <- f[["/original_data/intervals/indexed_contig_names"]]$read()
  ori_intervals_pos <- as.data.frame(f[["/original_data/intervals/transposed_index_start_end"]]$read())
  panel_contig_names <- f[["/panel/intervals/indexed_contig_names"]]$read()
  panel_intervals_pos <- as.data.frame(f[["/panel/intervals/transposed_index_start_end"]]$read())
  f$close_all()

  ori_intervals_pos[,1] <- ori_contig_names[as.numeric(ori_intervals_pos$V1) + 1]
  panel_intervals_pos[,1] <- panel_contig_names[as.numeric(panel_intervals_pos$V1) + 1]
  colnames(ori_intervals_pos) <- c("chrom", "start", "end")
  colnames(panel_intervals_pos) <- c("chrom", "start", "end")

  blacklist <- dplyr::anti_join(ori_intervals_pos, panel_intervals_pos, by = c("chrom", "start", "end"))



  cytoband <- read.table(cytoband,header=T)
  cytoband$chrom <- gsub("chr","",cytoband$chrom)
  cytoband$chromStart <- cytoband$chromStart + 1
  detectable_edge_df <- panel_intervals_pos %>%
    dplyr::rowwise() %>%
    dplyr::mutate( start_cytoband = FindCytoband( cytoband, chrom = chrom, pos = start ),
            end_cytoband = FindCytoband( cytoband, chrom = chrom, pos = end)) %>%
    dplyr::mutate( start_arm = gsub("\\d+.*","",start_cytoband),
            end_arm = gsub("\\d+.*","",start_cytoband))

  if( ! identical( detectable_edge_df$start_arm, detectable_edge_df$end_arm) ){
    arm_diff <- detectable_edge_df %>%
      dplyr::filter( start_arm != end_arm )
    arm_diff_cytoband <- lapply(1:nrow(arm_diff), function(x){
      chromosome <- arm_diff[x,"chrom"] %>% as.character()
      start <- arm_diff[x,"start"] %>% as.numeric()
      end <- arm_diff[x,"end"] %>% as.numeric()
      median_cov <- arm_diff[x,"median_cov"] %>% as.numeric()
      cytoband_tmp <- cytoband %>%
        dplyr::filter( chrom == chromosome ) %>%
        dplyr::filter( chromEnd >= start  ) %>%
        dplyr::filter( chromStart <= end ) %>%
        dplyr::mutate( median_cov = median_cov,
                start_arm = gsub("\\d+.*","",name),
                end_arm = start_arm)
      colnames(cytoband_tmp) <- c("chrom","start","end","start_cytoband", "gieStain", "median_cov", "start_arm", "end_arm")
      cytoband_tmp[1,"start"] <- start
      cytoband_tmp[nrow(cytoband_tmp),"end"] <- end
      cytoband_tmp$end_cytoband <- cytoband_tmp$start_cytoband
      cytoband_tmp <- cytoband_tmp[,colnames(arm_diff)]
      return(cytoband_tmp)
    })
    arm_diff_cytoband <- do.call(rbind, arm_diff_cytoband)
    detectable_edge_df <- detectable_edge_df %>%
      dplyr::filter( start_arm != end_arm) %>%
      dplyr::bind_rows(arm_diff_cytoband) %>%
      dplyr::arrange(chrom, start )
  }

  detectable_edge_df <- detectable_edge_df %>%
    dplyr::group_by( chrom, start_arm ) %>%
    dplyr::summarize( chromStart = min(start),
               chromEnd = max(end),
               first_name = start_cytoband[which.min(start)],
               last_name = end_cytoband[which.max(end)] ) %>%
    tidyr::pivot_wider(
      names_from = start_arm,
      values_from = c(chromStart, chromEnd, first_name, last_name) )

  detectable_edge_df <- detectable_edge_df %>%
    dplyr::mutate( chrom = factor(chrom, levels = c(as.character(1:22), "X", "Y"))) %>%
    dplyr::arrange(chrom ) %>%
    dplyr::rename_with(~stringr::str_replace(., "^chrom(Start|End)_([pq])$", "\\2_chrom\\1")) %>%
    dplyr::rename_with(~stringr::str_replace(., "^first_name_([pq])$", "\\1_first_name")) %>%
    dplyr::rename_with(~stringr::str_replace(., "^last_name_([pq])$", "\\1_last_name"))

  colnames(detectable_edge_df)[1] <- "Chrom"
  detectable_edge_df <- detectable_edge_df[,c("Chrom",	"p_chromStart",	"p_chromEnd",	"p_first_name",	"p_last_name",
                                                "q_chromStart",	"q_chromEnd",	"q_first_name",	"q_last_name")]
  if(gender == "female"){
    blacklist <- blacklist %>% dplyr::filter(chrom !="Y")
    panel_intervals_pos <- panel_intervals_pos %>% dplyr::filter( chrom !="Y")
    detectable_edge_df <- detectable_edge_df %>% dplyr::filter( Chrom !="Y")
  }
  utils::write.table(blacklist, file = blacklist_bed, quote = FALSE, row.names = FALSE, sep = "\t")
  utils::write.table(panel_intervals_pos, file = whitelist_bed, quote = FALSE, row.names = FALSE, sep = "\t")
  utils::write.table(detectable_edge_df, file = detectable_edge, quote = FALSE, row.names = FALSE, sep = "\t")


}



