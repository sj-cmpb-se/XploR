#' Annotate Segments with Overlapping Genes
#'
#' For each segment, identifies overlapping genes and adds gene names and the number of overlapping genes to the segment data.
#'
#' @param gene Character. Path to a gene annotation file (tab-delimited, columns: chrom, loc.start, loc.end, Gene).
#' @param segment Data frame. Segment data with columns \code{chrom}, \code{loc.start}, and \code{loc.end}.
#'
#' @return A data frame with the original segment columns plus \code{Gene} (comma-separated gene names) and \code{Gene_count}.
#'
#' @importFrom data.table fread
#' @importFrom dplyr group_by summarize
#' @importFrom GenomicRanges makeGRangesFromDataFrame findOverlaps
#' @examples
#' # segment <- data.frame(chrom = "1", loc.start = 1000, loc.end = 2000)
#' # AnnotateGene("genes.txt", segment)
#' @export
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


#' Generate Cytoband-Based Annotation String for a CNV Segment
#'
#' Returns a cytoband-based annotation string for a CNV segment, based on the segment's position and its proximity to p and q arm detectable region boundaries.
#'
#' The function determines the annotation as follows:
#' - If the segment start and end are in the same cytoband, returns a single-band range.
#' - If the segment spans multiple cytobands, returns the range from start to end cytoband.
#' - If the segment is close to both p-arm telomere and centromere (or at the p-arm boundaries), returns a p-arm event only.
#' - If the segment is close to both q-arm centromere and telomere (or at the q-arm boundaries), returns a q-arm event only (or chromosome-level event for acrocentric chromosomes).
#' - If the segment is close to both p-arm and q-arm telomeres (or covers the whole chromosome), returns a chromosome-level event.
#' The gap size smaller than 5MB are defined as "close".
#'
#' For acrocentric chromosomes, only q-arm boundaries are considered.
#'
#' @param chrom Character. Chromosome name.
#' @param Start, End Numeric. Segment start and end positions (base pairs).
#' @param CNF Numeric. Fractional copy number.
#' @param CN Integer. Integer copy number.
#' @param p_chromStart, p_chromEnd Numeric. p-arm detectable region start and end positions.
#' @param p_first_name, p_last_name Character. Cytoband names for the p-arm start and end.
#' @param q_chromStart, q_chromEnd Numeric. q-arm detectable region start and end positions.
#' @param q_first_name, q_last_name Character. Cytoband names for the q-arm start and end.
#' @param p_gap_to_tel, p_gap_to_cen, q_gap_to_tel, q_gap_to_cen Numeric. Gaps between segment and detectable region boundaries (from \code{CalculateGaps}).
#' @param cytoband Data frame. Cytoband annotation table.
#'
#' @return Character. A cytoband-based annotation string for the segment.
#'
#' @details
#' The annotation string is determined by the size of the gaps between the segment and p/q arm detectable region boundaries and the cytoband names at the segment start and end.
#' Events are reported as a cytoband range, p-arm only, q-arm only, or chromosome-level event.
#'
#' @examples
#' # GenerateISCN("1", 10000, 50000, 2.1, 2, 0, 25000, "p36.33", "p11.2",
#' #               25001, 100000, "q11.1", "q44", 1000, 1000, 1000, 1000, cytoband)
#'
#' @export
GenerateISCN <- function(chrom, Start, End, CNF, CN,
                          p_chromStart, p_chromEnd,
                          p_first_name, p_last_name,
                          q_chromStart, q_chromEnd,
                          q_first_name, q_last_name,
                          p_gap_to_tel, p_gap_to_cen,
                          q_gap_to_tel, q_gap_to_cen,
                          cytoband) {
  acrocentric_chr <- c("13", "14", "15", "21", "22")
  start_band <- FindCytoband(cytoband = cytoband, chrom = chrom, pos = Start)
  end_band <- FindCytoband(cytoband = cytoband, chrom = chrom, pos = End)
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

#' Find Cytoband Name for Genomic Position
#'
#' Returns the cytoband name for a given chromosome and position by searching a cytoband annotation data frame.
#'
#' @param cytoband Data frame. Cytoband annotation with columns \code{chrom}, \code{chromStart}, \code{chromEnd}, and \code{name}.
#' @param chrom Character. Chromosome name (e.g., "1", "X").
#' @param pos Numeric or character. Genomic position (base pair).
#'
#' @return Character. The cytoband name(s) containing the position.
#'
#' @examples
#' # cytoband <- data.frame(chrom = "1", chromStart = 1, chromEnd = 1000000, name = "p36.33")
#' # FindCytoband(cytoband, chrom = "1", pos = 500000)
#'
#' @export
FindCytoband <- function(cytoband, chrom, pos) {
  pos <- as.numeric(pos)
  site_name <- cytoband[cytoband$chrom == chrom & (pos >= cytoband$chromStart & pos <= cytoband$chromEnd), "name"]
  return(site_name)
}

#' Calculate Gaps Between Segment and Cytoband-Defined Chromosome Landmarks
#'
#' Calculates the distances (gaps) from a segment's start and end positions to the defined start and end sites of the p and q arms of a chromosome.
#' For acrocentric chromosomes, only the q arm start and end sites are considered.
#'
#' @param Chrom Character. Chromosome name (e.g., "1", "X").
#' @param Start Numeric. Segment start position (base pair).
#' @param End Numeric. Segment end position (base pair).
#' @param p_chromStart Numeric. Start position of the p arm detectable region.
#' @param p_chromEnd Numeric. End position of the p arm detectable region.
#' @param q_chromStart Numeric. Start position of the q arm detectable region.
#' @param q_chromEnd Numeric. End position of the q arm detectable region.
#'
#' @return A named list with elements:
#'   \item{p_gap_to_tel}{Gap from segment start to p arm telomere (or NA).}
#'   \item{p_gap_to_cen}{Gap from segment end to p arm centromere (or NA).}
#'   \item{q_gap_to_tel}{Gap from segment end to q arm telomere (or NA).}
#'   \item{q_gap_to_cen}{Gap from segment start to q arm centromere (or NA).}
#'
#' @details
#' For acrocentric chromosomes, p arm calculations are omitted and only q arm start/end sites are considered.
#'
#' @examples
#' # CalculateGaps("1", 10000, 50000, 0, 25000, 25001, 100000)
#'
#' @export
CalculateGaps <- function(Chrom, Start, End, p_chromStart, p_chromEnd, q_chromStart, q_chromEnd) {
  acrocentric_chr <- c("13", "14", "15", "21", "22")
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

#' Annotate CNV Segments with Cytoband and Gene Information
#'
#' Annotates CNV segments with cytoband and gene overlap information.
#'
#' For each segment, the function calculates the gaps between the segment boundaries and the detectable p and q arm region start/end sites (using cytoband and whitelist edge files). For acrocentric chromosomes, only the q arm boundaries are considered. Based on the size of these gaps, the function determines the event type and returns:
#' - a cytoband range (combination of start and stop cytoband),
#' - a p-arm or q-arm only event,
#' - or, if the event spans the whole chromosome, just the chromosome name (chromosome-level event).
#'
#' Gene overlap information is also added for each segment.
#'
#' @param input Path to merged segment file (TSV).
#' @param out_dir Output directory.
#' @param prefix Output file prefix.
#' @param cytoband Path to cytoband annotation file (TSV).
#' @param whitelist_edge Path to detectable edge file for cytoband annotation.
#' @param gene Path to gene annotation file.
#'
#' @return Invisibly returns the annotated data frame with cytoband and gene annotation columns.
#'
#' @importFrom dplyr left_join rowwise mutate
#' @importFrom data.table fread
#' @importFrom tidyr unnest_wider
#' @examples
#' AnnotateSegments(
#'   input = "sample_CNV_final.tsv",
#'   out_dir = "results/",
#'   prefix = "Sample1",
#'   cytoband = "data/cytoBand.txt",
#'   whitelist_edge = "data/whitelist.txt",
#'   gene = "data/gene_anno.txt"
#' )
#' @export
AnnotateSegments <- function(input, out_dir, prefix, cytoband, whitelist_edge, gene) {
  whitelist_edge <- read.table(whitelist_edge, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cytoband_df <- read.table(cytoband, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cytoband_df$chrom <- gsub("chr", "", cytoband_df$chrom)
  acrocentric_chr <- c("13", "14", "15", "21", "22")
  df <- data.table::fread(input)
  colnames(whitelist_edge)[1] <- "chrom"
  colnames(df)[c(1:3)] <- c("chrom", "loc.start", "loc.end")

  df_iscn <- df %>%
    dplyr::left_join(whitelist_edge, by = "chrom") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(gaps = list(CalculateGaps(
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
    dplyr::mutate(ISCN = GenerateISCN(
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
