#' Read and Combine Panel of Normal AI Files
#'
#' Reads multiple allelic imbalance (AI) files from a panel of normal (PoN), processes and combines the results into a single data frame.
#'
#' @param ai_pon_file Character. Path to a text file listing the paths to all PoN AI files (one file path per line).
#' @param aitype Character. Type of AI input file. Passed to \code{ReadAI()} (e.g., \code{"gatk"}, \code{"dragen"}, \code{"other"}).
#' @param minsnpcov Integer. Minimum SNP coverage to include a site in the AI calculation. Passed to \code{ReadAI()}.
#' @param gender Character. Sample gender, either \code{"female"} or \code{"male"}. This value is passed to \code{ReadAI()}.
#'
#' @return A data frame combining all AI data from the PoN files.
#'
#' @details
#' Each file path listed in \code{ai_pon_file} is read and processed using \code{ReadAI()}, with the resulting data frames combined by row.
#'
#' @importFrom dplyr bind_rows
#' @export

ReadPonAI <- function( ai_pon_file, aitype, minsnpcov, gender ){
  ai_pon <- read.table(ai_pon_file,sep="\t")
  aipondt <- lapply(ai_pon$V1,function(x){
    print(paste0("Loading AI PON files.......  ",  x ))
    baf <- ReadAI(aitype = aitype,ballele = x, minsnpcov = minsnpcov, gender = gender)
    sampleID <- gsub(".*\\/|\\..*","",x)
    baf$sampleID <- sampleID
    return(baf)
  })

  aipondt <- do.call(rbind, aipondt)

  return(aipondt)

}

#' Choose Optimal Number of Depth Strata (Bins) for Panel of Normals
#'
#' Determines the optimal number of bins (strata) for grouping samples by SNP median depth, ensuring each bin is sufficiently populated and diverse for downstream modeling.
#'
#' @param normals_dt Data frame or data.table. Must has depth information of each bin.
#'
#' @return Integer. The chosen number of bins that satisfies the population and diversity constraints.
#'
#' @details
#' This function iteratively tests different bin (strata) numbers to find the optimal value based on SNP counts and sequencing depth.
#'
#' @importFrom data.table as.data.table
#' @importFrom stats quantile
#' @export
ChooseNbins <- function(normals_dt,
                        target_bins = 8, min_rows_per_stratum = 2000,
                        min_unique_bins = 100, max_bins = 10, min_bins = 3) {
  dt <- data.table::as.data.table(normals_dt)
  N <- nrow(dt)
  # Try decreasing number of bins until every stratum is sufficiently populated
  for (k in pmin(target_bins, max_bins):min_bins) {
    db <- cut(dt$psb_snp_median_depth,
              breaks = stats::quantile(dt$psb_snp_median_depth, probs = seq(0, 1, length.out = k + 1),
                                       na.rm = TRUE, names = FALSE),
              include.lowest = TRUE, labels = FALSE)
    tab_rows <- table(db)
    # how many distinct genomic bins in each stratum?
    lib <- dt[, .(bin, db)][, .N, by = .(db, bin)][, .N, by = db]$N
    if (all(tab_rows >= min_rows_per_stratum) && all(lib >= min_unique_bins)) return(k)
  }
  return(min_bins)
}

#' Estimate Beta-Binomial Over-Dispersion Parameter (Theta)
#'
#' Estimates the over-dispersion parameter (\code{theta}) for a beta-binomial model stratified by sequencing depth, using panel of normal (PoN) BAF data.
#'
#' @param normals_dt Data frame or tibble. Panel of normal data.
#' @param pon_ref Data frame or tibble. Reference PoN data.
#' @param n_bins Integer. Number of depth strata to use.
#'
#' @return
#'   Returns a list with elements:
#'   \item{theta_table}{Data frame with columns \code{depth_bin} and \code{theta} (per-stratum estimates).}
#'   \item{breaks}{Numeric vector of depth quantile breakpoints.}
#'
#' @details
#' The function compares observed BAF variance to expected binomial variance across bins and depths, and estimates theta as:
#' \deqn{\theta = \frac{\mathrm{Var}(p)}{\mathrm{E}[p(1-p)/d]} - 1}{theta = (Var(p) / E[p(1-p)/d] - 1) / (median depth - 1)}
#' where \eqn{p} is the observed BAF, and \eqn{d} is the median depth.
#'
#' @importFrom dplyr left_join filter group_by summarise mutate ungroup
#' @export
EstimateTheta <- function(normals_dt, pon_ref, n_bins) {

  dt <- normals_dt %>%
    dplyr::left_join(pon_ref %>%
                       dplyr::select(bin, pon_mean_snp_baf), by = "bin") %>%
    dplyr::filter(is.finite(pon_mean_snp_baf), is.finite(pon_mean_snp_baf))

    # Depth-quantile strata
    probs <- seq(0, 1, length.out = n_bins + 1)
    brks  <- quantile(dt$psb_snp_median_depth, probs = probs, na.rm = TRUE, names = FALSE)
    brks  <- unique(brks)  # guard against duplicated breaks
    if (length(brks) < 2L) brks <- range(depth_vec, na.rm = TRUE)

    dt2 <- dt %>%
      dplyr::mutate(depth_bin = cut(psb_snp_median_depth, breaks = brks, include.lowest = TRUE, labels = FALSE))

    tmp <- dt2 %>%
      dplyr::group_by(bin, depth_bin) %>%
      dplyr::summarise(
        var_baf   = var(psb_snp_baf, na.rm = TRUE),
        mean_binv = mean(pon_mean_snp_baf * (1 - pon_mean_snp_baf) / psb_snp_median_depth, na.rm = TRUE),
        pon_median_snp_median_depth = median(psb_snp_median_depth, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        theta_bin = pmax(
          (var_baf / pmax(mean_binv, 1e-12) - 1) / pmax(pon_median_snp_median_depth - 1, 1), 0
        )
      )

    theta_table <- tmp %>%
      dplyr::group_by(depth_bin) %>%
      dplyr::summarise(theta = median(theta_bin, na.rm = TRUE), .groups = "drop") %>%
      dplyr::mutate(theta = pmin(pmax(theta, 0), 0.2))  # clamp per stratum

    return(list(theta_table = theta_table, breaks = brks))

}

#' Process Panel of Normal AI Files and Estimate Beta-binomial Dispersion based on the depth
#'
#' Reads and bins allelic imbalance (AI) data from a panel of normal (PoN), summarizes bin-level statistics, and estimates beta-binomial dispersion (theta) for downstream modeling.
#'
#' @param ai_pon_file Character. Path to a text file listing PoN AI file paths (one per line).
#' @param aitype Character. Type of AI input file (passed to \code{ReadPonAI()}, e.g., \code{"gatk"}, \code{"dragen"}, \code{"other"}).
#' @param minsnpcov Integer. Minimum SNP coverage to include a site in the AI calculation. (default: 20)
#' @param maxgap Numeric. Maximum allowed gap between SNPs within a bin. (default: 2000000)
#' @param maxbinsize Numeric. Maximum allowed bin size (bp). (default: 5000000)
#' @param minbinsize Numeric. Minimum allowed bin size (bp). (default: 500000)
#' @param snpnum Integer. Target number of SNPs per bin. (default: )
#' @param output Character. Output directory for the processed PoN AI Rdata file.
#' @param prefix Character. Prefix for the output file.
#'
#' @return Invisibly returns \code{NULL}. Saves an \code{Rdata} file containing the processed PoN reference (\code{pon_ref}) and the estimated dispersion parameters (\code{theta_fit}).
#'
#' @details
#' This function reads all PoN AI files, bins the AI data using \code{BinMaf()}, summarizes bin-level BAF/MAF and depth statistics, and estimates beta-binomial over-dispersion (\code{theta}) stratified by depth. The reference table and theta estimates are saved as an \code{Rdata} file for use in downstream analysis.
#'
#' @importFrom dplyr group_by summarise mutate left_join
#' @importFrom stats median
#' @export
PONAIprocess <- function( ai_pon_file, aitype, minsnpcov = 20, output,
                          prefix, maxgap = 2000000, maxbinsize = 5000000,
                          minbinsize = 500000, snpnum = 30, gender
                          ){

  aipondt <- ReadPonAI(ai_pon_file = ai_pon_file,
                       aitype = aitype,
                       minsnpcov = minsnpcov,
                       gender = gender )

  normals_dt <- BinMaf(data = aipondt,
                       datatype = "pon",
                       maxgap = maxgap,
                       maxbinsize = maxbinsize,
                       minbinsize = minbinsize,
                       snpnum = snpnum,
                       minsnpcov = minsnpcov)


  pon_ref <- normals_dt %>%
    dplyr::group_by(bin, Chromosome, Start, End, size) %>%
    dplyr::summarise(
      pon_mean_snp_baf  = median(psb_snp_baf, na.rm = TRUE),
      pon_mean_snp_median_baf = median(psb_snp_median_baf, na.rm = TRUE),
      pon_mean_snp_median_maf = median(psb_snp_median_maf, na.rm = TRUE),
      pon_mafs = paste( psb_snp_mafs , collapse = ","),
      pon_depth_median  = median(psb_snp_median_depth, na.rm = TRUE),
      n_normals     = sum(is.finite(psb_snp_baf)),
      .groups = "drop"
    )

  ### Estimate beta-binomial dispersion theta based on depth
  n_bins <- ChooseNbins(normals_dt)
  theta_fit <- EstimateTheta(normals_dt, pon_ref, n_bins = n_bins)

  outFile <- paste0( output, "/", prefix, "_PON_AI.Rdata")
  save(pon_ref,theta_fit,file=outFile)


}

