#' Annotate CNV Segments with ISCN and Gene Information
#'
#' Annotates CNV segments with ISCN (International System for Human Cytogenomic Nomenclature) and gene overlap information.
#'
#' @param input Path to merged segment file (TSV).
#' @param out_dir Output directory.
#' @param prefix Output file prefix.
#' @param cytoband Path to cytoband annotation file.
#' @param whitelist_edge Path to detectable edge file for ISCN.
#' @param gene Path to gene annotation file.
#'
#' @return Invisibly returns the annotated data frame.
#'
#' @importFrom dplyr left_join rowwise mutate group_by summarize arrange filter
#' @importFrom data.table fread
#' @importFrom tidyr unnest_wider
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @export
AnnotateSegments <- function(input, out_dir, prefix, cytoband, whitelist_edge, gene) {
  whitelist_edge <- read.table(whitelist_edge, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cytoband_df <- read.table(cytoband, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cytoband_df$chrom <- gsub("chr", "", cytoband_df$chrom)
  acrocentric_chr <- c("13", "14", "15", "21", "22")
  df <- data.table::fread(input)
  colnames(whitelist_edge)[1] <- "chrom"
  colnames(df)[c(1:3)] <- c("chrom", "loc.start", "loc.end")

  calculate_gaps <- function(Chrom, Start, End, p_chromStart, p_chromEnd, q_chromStart, q_chromEnd) {
    Start <- ifelse(Chrom %in% acrocentric_chr, max(c(q_chromStart, Start)), max(c(p_chromStart, Start)))
    End <- min(c(End, q_chromEnd))
    if (!is.na(p_chromEnd) & Start <= p_chromEnd) {
      p_gap_to_tel <- Start - p_chromStart
      p_gap_to_tel <- ifelse(p_gap_to_tel > 0, p_gap_to_tel, 0)
    } else { p_gap_to_tel <- NA }
    if (((!is.na(p_chromEnd) & Start >= p_chromEnd) | Chrom %in% acrocentric_chr) & Start <= q_chromEnd) {
      q_gap_to_cen <- Start - q_chromStart
      q_gap_to_cen <- ifelse(q_gap_to_cen > 0, q_gap_to_cen, 0)
    } else { q_gap_to_cen <- NA }
    if ((!is.na(p_chromStart) & End >= p_chromStart) & End <= q_chromStart) {
      p_gap_to_cen <- p_chromEnd - End
      p_gap_to_cen <- ifelse(p_gap_to_cen > 0, p_gap_to_cen, 0)
    } else { p_gap_to_cen <- NA }
    if (((!is.na(p_chromEnd) & End >= p_chromEnd) | Chrom %in% acrocentric_chr) & End <= q_chromEnd) {
      q_gap_to_tel <- q_chromEnd - End
      q_gap_to_tel <- ifelse(q_gap_to_tel > 0, q_gap_to_tel, 0)
    } else { q_gap_to_tel <- NA }
    list("p_gap_to_tel" = p_gap_to_tel,
         "p_gap_to_cen" = p_gap_to_cen,
         "q_gap_to_tel" = q_gap_to_tel,
         "q_gap_to_cen" = q_gap_to_cen)
  }

  findcytoband <- function(cytoband, chrom, pos) {
    pos <- as.numeric(pos)
    site_name <- cytoband[cytoband$chrom == chrom & (pos >= cytoband$chromStart & pos <= cytoband$chromEnd), "name"]
    return(site_name)
  }

  generate_ISCN <- function(chrom, Start, End, CNF, CN,
                            p_chromStart, p_chromEnd,
                            p_first_name, p_last_name,
                            q_chromStart, q_chromEnd,
                            q_first_name, q_last_name,
                            p_gap_to_tel, p_gap_to_cen,
                            q_gap_to_tel, q_gap_to_cen,
                            cytoband) {
    start_band <- findcytoband(cytoband = cytoband, chrom = chrom, pos = Start)
    end_band <- findcytoband(cytoband = cytoband, chrom = chrom, pos = End)
    copy_number <- ifelse(is.na(CN), round(CNF, digits = 2), CN)
    if (start_band == end_band) {
      ISCN = paste0(chrom, start_band, "(", Start, "_", End, ")X", copy_number)
    } else {
      ISCN = paste0(chrom, start_band, end_band, "(", Start, "_", End, ")X", copy_number)
    }
    if (!is.na(p_gap_to_cen) & !is.na(p_gap_to_tel)) {
      if ((p_gap_to_tel < 5000000 | start_band == p_first_name) &
          (p_gap_to_cen < 5000000 | end_band == p_last_name)) {
        ISCN = paste0(chrom, "p(", Start, "_", End, ")X", copy_number)
      }
    }
    if (!is.na(q_gap_to_cen) & !is.na(q_gap_to_tel)) {
      if ((q_gap_to_cen < 5000000 | start_band == q_first_name) &
          (q_gap_to_tel < 5000000 | end_band == q_last_name)) {
        if (chrom %in% acrocentric_chr) {
          ISCN = paste0(chrom, "(", Start, "_", End, ")X", copy_number)
        } else {
          ISCN = paste0(chrom, "q(", Start, "_", End, ")X", copy_number)
        }
      }
    }
    if (!is.na(p_gap_to_tel) & !is.na(q_gap_to_tel)) {
      if ((p_gap_to_tel < 5000000 & q_gap_to_tel < 5000000) |
          (start_band == p_first_name & end_band == q_last_name)) {
        ISCN = paste0(chrom, "(", Start, "_", End, ")X", copy_number)
      }
    }
    return(ISCN)
  }

  AnnotateGene <- function(gene, segment) {
    gene = data.table::fread(gene)
    colnames(gene) <- c("chrom", "loc.start", "loc.end", "Gene")
    gene$chrom <- gsub("chr", "", gene$chrom)
    gene_ranges <- GenomicRanges::makeGRangesFromDataFrame(
      gene,
      seqnames.field = "chrom",
      start.field = "loc.start",
      end.field = "loc.end",
      keep.extra.columns = TRUE
    )
    segment_ranges <- GenomicRanges::makeGRangesFromDataFrame(
      segment,
      seqnames.field = "chrom",
      start.field = "loc.start",
      end.field = "loc.end",
      keep.extra.columns = TRUE
    )
    overlapping <- GenomicRanges::findOverlaps(segment_ranges, gene_ranges)
    overlapping <- as.data.frame(overlapping)
    overlapping[,"name"] <- gene[overlapping$subjectHits, "Gene"]
    format_gene <- function(x) {
      if (length(x) > 1) {
        x <- paste(x, collapse = ", ")
      }
      return(x)
    }
    gene_anno <- overlapping %>%
      dplyr::group_by(queryHits) %>%
      dplyr::summarize(Gene = format_gene(name), Gene_number = length(name))
    segment[gene_anno$queryHits, c("Gene", "Gene_count")] <- gene_anno[, -1]
    return(segment)
  }

  df_iscn <- df %>%
    dplyr::left_join(whitelist_edge, by = "chrom") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gaps = list(calculate_gaps(
      Chrom = chrom,
      Start = loc.start,
      End = loc.end,
      p_chromStart = p_chromStart,
      p_chromEnd = p_chromEnd,
      q_chromStart = q_chromStart,
      q_chromEnd = q_chromEnd
    ))) %>%
    tidyr::unnest_wider(gaps) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ISCN = generate_ISCN(
      chrom = chrom, Start = loc.start, End = loc.end,
      CNF = CNF_correct, CN = CN,
      p_chromStart = p_chromStart, p_chromEnd = p_chromEnd,
      p_first_name = p_first_name, p_last_name = p_last_name,
      q_chromStart = q_chromStart, q_chromEnd = q_chromEnd,
      q_first_name = q_first_name, q_last_name = q_last_name,
      p_gap_to_tel = p_gap_to_tel, p_gap_to_cen = p_gap_to_cen,
      q_gap_to_tel = q_gap_to_tel, q_gap_to_cen = q_gap_to_cen,
      cytoband = cytoband_df
    ))

  segment_gene <- AnnotateGene(gene = gene, segment = df_iscn[, c("chrom", "loc.start", "loc.end")])
  result <- df_iscn %>% dplyr::left_join(segment_gene, by = c("chrom", "loc.start", "loc.end"))

  OutFile <- file.path(out_dir, paste0(prefix, "_CNV_annotation.tsv"))
  write.table(result, file = OutFile, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(result)
}
