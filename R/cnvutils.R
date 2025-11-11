#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`


#' Calculate Segment Mean Adjusted by Scale Factor
#'
#' @param Segcov Numeric. The observed segment coverage.
#' @param max_L_mu Numeric. The scale factor (mu).
#' @param gender Character. "male" or "female".
#' @param Chromosome Character. Chromosome name (e.g., "1", "X", "Y").
#'
#' @return Numeric. The adjusted segment mean (log2 ratio).
#'
#' @examples
#' CalSM(Segcov = 150, max_L_mu = 1, gender = "female", Chromosome = "1")
#'
#' @export
CalSM <- function(Segcov, max_L_mu, gender, Chromosome) {
  diploid_cov <- 100 * max_L_mu
  uniploid_cov <- 50 * max_L_mu
  if (!is.na(Segcov)) {
    if (gender == "male" & Chromosome %in% c("X", "Y")) {
      Segment_Mean <- log2(Segcov / uniploid_cov)
    } else {
      Segment_Mean <- log2(Segcov / diploid_cov)
    }
  } else {
    Segment_Mean <- NA
  }
  return(Segment_Mean)
}


#' Correct CNF and MAF According to Purity
#'
#' @param chromosome Character. Chromosome name.
#' @param cov_segmentmean Numeric. Segment mean (log2 ratio).
#' @param MAF_observe Numeric. Observed minor allele frequency.
#' @param gender Character. "male" or "female".
#' @param purity Numeric. Tumor purity (fraction between 0 and 1).
#'
#' @return A list with \code{CNF_correct} and \code{MAF_correct}.
#'
#' @export
CorrectPurity <- function(chromosome, cov_segmentmean, MAF_observe, gender, purity) {
  if (gender == "male" & chromosome %in% c("X", "Y")) {
    CNF_observe = 2^cov_segmentmean
    CNF_tumor = (CNF_observe - 1 + purity) / purity
    MAF_tumor = MAF_observe / purity
    MAF_tumor = ifelse(MAF_tumor > 0.5, 0.5, MAF_tumor)
    MAF_tumor = ifelse(MAF_tumor < 0, 0, MAF_tumor)
  } else {
    CNF_observe = 2 * 2^cov_segmentmean
    CNF_tumor = (CNF_observe - 2 + 2 * purity) / purity
    MAF_tumor = (MAF_observe - 0.5 + 0.5 * purity) / purity
    MAF_tumor = ifelse(MAF_tumor > 0.5, 0.5, MAF_tumor)
    MAF_tumor = ifelse(MAF_tumor < 0, 0, MAF_tumor)
  }
  if (CNF_tumor < 0) { CNF_tumor <- 0 }
  correct = list(
    CNF_correct = CNF_tumor,
    MAF_correct = MAF_tumor
  )
  return(correct)
}

#' Assign Copy Number Call Without Model
#'
#' Assigns a copy number call ("GAIN", "LOSS", "REF", "GAINLOH", or "CNLOH") to a segment based on corrected copy number fraction (CNF), minor allele frequency (MAF), and probe characteristics, using user-defined thresholds.
#'
#' @param chromosome Character. Chromosome name (e.g., "1", "X", "Y").
#' @param CNF_correct Numeric. Corrected copy number fraction (CNF).
#' @param MAF_correct Numeric. Corrected minor allele frequency (MAF).
#' @param MAF_gmm_G Integer. Number of GMM clusters (not directly used in logic, but included for compatibility).
#' @param MAF_Probes Integer. Number of probes used for MAF estimation.
#' @param MAF_gmm_weight Numeric. GMM weight (not directly used in logic, but included for compatibility).
#' @param gender Character. Sample gender ("male" or "female").
#' @param callcovcutoff Numeric. Threshold for CNF to call gain or loss (default: 0.3).
#' @param callaicutoff Numeric. Threshold for MAF to call gain with loss of heterozygosity (GAINLOH) (default: 0.3).
#' @param minsnpcallaicutoff Integer. Minimum number of SNP probes to call CNLOH (default: 10).
#' @param gender Character. Gender information "male" or "female".
#'
#' @return Character. Copy number call: one of "GAIN", "LOSS", "REF", "GAINLOH", or "CNLOH".
#'
#' @details
#' Calls are made according to CNF and MAF thresholds:
#' - For sex chromosomes in males, thresholds are adjusted for haploid status.
#' - "GAINLOH" is called if CNF is high and MAF is below a dynamic cutoff.
#' - "CNLOH" is called if MAF is low and probe count meets the minimum threshold.
#'
#' @examples
#' CallwoModel("1", CNF_correct = 2.8, MAF_correct = 0.1, MAF_gmm_G = 2, MAF_Probes = 15, MAF_gmm_weight = 0.5, gender = "female")
#'
#' @export
CallwoModel <- function(chromosome, CNF_correct, MAF_correct, MAF_gmm_G, MAF_Probes, MAF_gmm_weight, callcovcutoff = 0.3, callaicutoff = 0.3, minsnpcallaicutoff = 10, gender ) {
  if ( gender == "male" && chromosome %in% c("X", "Y")) {
    if (CNF_correct >= (1 + callcovcutoff)) {
      Call <- "GAIN"
    } else if (CNF_correct < (1 - callcovcutoff)) {
      Call <- "LOSS"
    } else {
      Call <- "REF"
    }
  } else {
    if (CNF_correct >= (2 + callcovcutoff)) {
      MAF_expected <- 1 / ceiling(CNF_correct)
      gainlohcutoff <- MAF_expected - callaicutoff
      gainlohcutoff <- ifelse(gainlohcutoff < 0.1, 0.1, gainlohcutoff)
      if (MAF_correct <= gainlohcutoff && MAF_Probes >= 10) {
        Call <- "GAINLOH"
      } else {
        Call <- "GAIN"
      }
    } else if (CNF_correct < (2 - callcovcutoff)) {
      Call <- "LOSS"
    } else {
      if (MAF_correct <= 0.3 && MAF_Probes >= minsnpcallaicutoff) {
        Call <- "CNLOH"
      } else {
        Call <- "REF"
      }
    }
  }
  return(Call)
}

#' Round Copy Number to Integer Based on Call and Gender
#'
#' @param gender Character. "male" or "female".
#' @param Chrom Character. Chromosome name.
#' @param Call Character. CNV call.
#' @param CNF Numeric. Corrected CNF.
#'
#' @return Integer. Rounded copy number.
#'
#' @export
RoundCN <- function(gender, Chrom, Call, CNF) {
  diff_cn <- abs(round(CNF) - CNF)
  if (diff_cn < 0.3) {
    CN <- round(CNF)
  } else {
    if (Call == "LOSS") {
      CN <- floor(CNF)
    } else if (Call %in% c("GAIN", "GAINLOH")) {
      CN <- ceiling(CNF)
    } else {
      if (gender == "male" && Chrom %in% c("X", "Y")) {
        CN <- 1
      } else {
        CN <- 2
      }
    }
  }
  return(CN)
}



#' Assign Copy Number Call With Model (Major/Minor)
#'
#' @param minor Integer. Minor allele copy number.
#' @param major Integer. Major allele copy number.
#' @param CNF_correct Numeric. Corrected CNF.
#' @param Chromosome Character. Chromosome name.
#' @param gender Character. "male" or "female".
#'
#' @return Character. Copy number call: "GAIN", "LOSS", "REF", "GAINLOH", or "CNLOH".
#'
#' @export
CallWTModel <- function(minor, major, CNF_correct, Chromosome, gender) {
  if (gender == "male" & Chromosome %in% c("X", "Y")) {
    if (CNF_correct >= 1.4) {
      Call <- "GAIN"
    } else if (CNF_correct < 0.6) {
      Call <- "LOSS"
    } else {
      Call <- "REF"
    }
  } else if (gender == "female" & Chromosome == "X") {
    if (CNF_correct >= 2.4) {
      Call <- "GAIN"
    } else if (CNF_correct <= 1.6) {
      Call <- "LOSS"
    } else {
      Call <- "REF"
    }
  } else {
    if (!is.na(minor) & !is.na(major)) {
      if (minor + major == 2) {
        if (major == minor) {
          Call <- "REF"
        } else if (major > minor) {
          Call <- "CNLOH"
        }
      } else if (minor + major > 2) {
        if (major - minor >= 2) {
          Call <- "GAINLOH"
        } else {
          Call <- "GAIN"
        }
      } else if (minor + major < 2) {
        Call <- "LOSS"
      }
    } else {
      if (CNF_correct >= 2.3) {
        Call <- "GAIN"
      } else if (CNF_correct < 1.7) {
        Call <- "LOSS"
      } else {
        Call <- "REF"
      }
    }
  }
  return(Call)
}
